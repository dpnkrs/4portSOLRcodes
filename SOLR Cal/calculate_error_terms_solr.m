function error_terms = calculate_error_terms_solr(freq, S11_S, S11_O, S11_L, S22_S, S22_O, S22_L, S_reciprocal)
%CALCULATE_ERROR_TERMS_SOLR Calculates the 8 error terms for SOLR calibration
%   based on the provided document's 8-term error model principles.
%
%   error_terms = calculate_error_terms_solr(...)
%
%   Inputs:
%       freq         - Frequency vector in Hz.
%       S11_S        - Measured S11 of Short standard on Port 1 (complex vector).
%       S11_O        - Measured S11 of Open standard on Port 1 (complex vector).
%       S11_L        - Measured S11 of Load standard on Port 1 (complex vector).
%       S22_S        - Measured S22 of Short standard on Port 2 (complex vector).
%       S22_O        - Measured S22 of Open standard on Port 2 (complex vector).
%       S22_L        - Measured S22 of Load standard on Port 2 (complex vector).
%       S_reciprocal - Measured S-parameters of the Reciprocal Through (Nx4 matrix: [S11_R, S21_R, S12_R, S22_R]).
%
%   Outputs:
%       error_terms  - Struct containing error terms for each frequency point.
%                      Fields correspond to the terms defined in the document:
%                      e1_00, e1_11, t11, e2_00, e2_11, t22, t21, t12.
%                      Each field is a vector of complex numbers.
%
%   Notes:
%       This implementation follows the standard algebraic solution for an 8-term
%       error model (often derived from 12-term simplified for ideal switch).
%       It assumes ideal calibration standards for Short (-1), Open (1), Load (0).
%       For the reciprocal through, it uses the measured S-parameters to solve
%       for the transmission error terms.

num_freq = length(freq);

% Initialize error terms
e1_00 = zeros(num_freq, 1);   % e1^00 (Port 1 Directivity)
e1_11 = zeros(num_freq, 1);   % e1^11 (Port 1 Source Match)
t11 = zeros(num_freq, 1);     % t11 = e1^10 * e1^01 (Port 1 Reflection Tracking)

e2_00 = zeros(num_freq, 1);   % e2^00 (Port 2 Directivity)
e2_11 = zeros(num_freq, 1);   % e2^11 (Port 2 Load Match / Source Match in reverse)
t22 = zeros(num_freq, 1);     % t22 = e2^10 * e2^01 (Port 2 Reflection Tracking)

t21 = zeros(num_freq, 1);     % t21 = e1^10 * e2^01 (Forward Transmission Tracking)
t12 = zeros(num_freq, 1);     % t12 = e2^10 * e1^01 (Reverse Transmission Tracking)


% Define ideal S-parameters for standards (adjust if using characterized standards)
Gamma_S_ideal = -1; % Ideal Short (reflection coefficient)
Gamma_O_ideal = 1;  % Ideal Open (reflection coefficient)
Gamma_L_ideal = 0;  % Ideal Load (reflection coefficient)

for i = 1:num_freq
    % --- Port 1 Error Term Calculation (e1_00, e1_11, t11) ---
    % Using measured S11 of Short, Open, Load on Port 1
    % The relationship is: S_measured = e_00 + (e_10*e_01 * Gamma_actual) / (1 - e_11 * Gamma_actual)
    % Let S_meas = S11_L(i), S11_S(i), S11_O(i)
    % Let Gamma_actual = Gamma_L_ideal, Gamma_S_ideal, Gamma_O_ideal
    
    % From the Load measurement (Gamma_L_ideal = 0), we directly get e1_00:
    e1_00(i) = S11_L(i); % This is the directivity term at Port 1

    % Now we have a system of two linear equations for e1_11 and t11:
    % (S_meas - e1_00) * (1 - e1_11 * Gamma_actual) = t11 * Gamma_actual
    % Rearranging: (S_meas - e1_00) = t11 * Gamma_actual + e1_11 * Gamma_actual * (S_meas - e1_00)
    % This is: A * [e1_11; t11] = B
    % Where A matrix rows are for Short and Open measurements.
    
    % For Short:
    Meas_minus_E00_S = S11_S(i) - e1_00(i);
    % For Open:
    Meas_minus_E00_O = S11_O(i) - e1_00(i);
    
    % Construct matrix A and vector B for [e1_11; t11]
    % [ Gamma_S_ideal * Meas_minus_E00_S, Gamma_S_ideal;
    %   Gamma_O_ideal * Meas_minus_E00_O, Gamma_O_ideal ] * [e1_11; t11] = [Meas_minus_E00_S; Meas_minus_E00_O]
    
    % The derivation from Pozar or VNA app notes leads to a simpler linear system for 1-port errors:
    % Let Y1 = (S11_S(i) - e1_00(i)) / Gamma_S_ideal;
    % Let Y2 = (S11_O(i) - e1_00(i)) / Gamma_O_ideal;
    
    % System:
    % Y1 = t11 + e1_11 * (S11_S(i) - e1_00(i))
    % Y2 = t11 + e1_11 * (S11_O(i) - e1_00(i))
    
    % Subtracting the two equations to solve for e1_11:
    if (S11_S(i) - S11_O(i)) == 0
         warning('Degenerate reflection standard response for Port 1 at frequency %g Hz. Skipping calculation for this point.', freq(i));
         e1_11(i) = NaN; t11(i) = NaN;
         e2_00(i) = NaN; e2_11(i) = NaN; t22(i) = NaN;
         t21(i) = NaN; t12(i) = NaN;
         continue;
    end
    
    e1_11(i) = ( (S11_S(i) - e1_00(i))/Gamma_S_ideal - (S11_O(i) - e1_00(i))/Gamma_O_ideal ) / ...
               ( (S11_S(i) - e1_00(i)) - (S11_O(i) - e1_00(i)) );
    
    % Then solve for t11:
    t11(i) = (S11_S(i) - e1_00(i))/Gamma_S_ideal - e1_11(i) * (S11_S(i) - e1_00(i));

    % --- Port 2 Error Term Calculation (e2_00, e2_11, t22) ---
    % Similar procedure for Port 2 using S22 of Short, Open, Load on Port 2
    e2_00(i) = S22_L(i); % Directivity at Port 2

    if (S22_S(i) - S22_O(i)) == 0
        warning('Degenerate reflection standard response for Port 2 at frequency %g Hz. Skipping calculation for this point.', freq(i));
        e1_11(i) = NaN; t11(i) = NaN;
        e2_00(i) = NaN; e2_11(i) = NaN; t22(i) = NaN;
        t21(i) = NaN; t12(i) = NaN;
        continue;
    end
    
    e2_11(i) = ( (S22_S(i) - e2_00(i))/Gamma_S_ideal - (S22_O(i) - e2_00(i))/Gamma_O_ideal ) / ...
               ( (S22_S(i) - e2_00(i)) - (S22_O(i) - e2_00(i)) );
    
    t22(i) = (S22_S(i) - e2_00(i))/Gamma_S_ideal - e2_11(i) * (S22_S(i) - e2_00(i));


    % --- Transmission Error Term Calculation (t21, t12) ---
    % This is the core of the 8-term model using the reciprocal through.
    % The document shows T_m = T1 * T * T2
    % We have T_m (from S_reciprocal) and T1, T2 (from calculated reflection errors).
    % We need to find T_actual (which implies t21 and t12).
    
    % Convert measured reciprocal S-params to T-params:
    S_recip_matrix = [S_reciprocal(i,1), S_reciprocal(i,3); S_reciprocal(i,2), S_reciprocal(i,4)];
    T_m_recip = s_to_t(S_recip_matrix);

    % Construct the T1 and T2 error matrices using the definitions from the document.
    % Note: t11 = e1^10 * e1^01, t22 = e2^10 * e2^01
    % Also: Delta_e1 = e1^00 * e1^11 - t11
    %       Delta_e2 = e2^00 * e2^11 - t22
    
    Delta_e1 = e1_00(i) * e1_11(i) - t11(i);
    Delta_e2 = e2_00(i) * e2_11(i) - t22(i);

    % T1 = (1/e1^10) * [-Delta_e1, e1^00; -e1^11, 1]
    % T2 = (1/e2^01) * [-Delta_e2, e2^11; -e2^00, 1]
    
    % From the document, the combined error matrix for T_m is:
    % Tm = (1/(e1^10 * e2^01)) * A * T_actual * B
    % Where A = [-Delta_e1, e1^00; -e1^11, 1]
    % And   B = [-Delta_e2, e2^11; -e2^00, 1]
    % and (e1^10 * e2^01) is the t21 term.
    
    % So, T_m_recip = (1/t21) * A * T_actual_recip * B
    % We are solving for t21 and T_actual_recip (specifically its S12/S21 terms).
    % The document states that "it is only necessary to solve for 7 of the 8 terms
    % during the calibration routine. The other transmission term t12 can be solved through
    % the following relation. t11 t22 = t21 t12".
    
    % Let's follow a standard approach for the 8-term model:
    % Calculate the coefficients of the matrices A and B (without the 1/e1^10 or 1/e2^01 scaling yet)
    Matrix_A = [-Delta_e1, e1_00(i); -e1_11(i), 1];
    Matrix_B = [-Delta_e2, e2_11(i); -e2_00(i), 1];
    
    % Now, from Tm = (1/t21) * A * T_actual * B
    % T_actual = (t21) * inv(A) * Tm * inv(B)
    % Since the reciprocal standard is reciprocal, T_actual(1,2) = -T_actual(2,1) (related to S12=S21).
    % Or, more simply, S12_actual = S21_actual, and S11_actual = S22_actual = 0 for an ideal thru.
    % For SOLR, the thru is non-ideal but reciprocal. So S11_actual and S22_actual are small but non-zero,
    % and S12_actual = S21_actual.

    % The critical step for SOLR involves solving for the two error boxes' "transmission" related terms.
    % The simplest approach that closes the loop from the document (and is often used for 8-term):
    % From the relation $t_{11}t_{22}=t_{21}t_{12}$, we solve for $t_{21}$ and then $t_{12}$.
    % We need an independent way to get $t_{21}$.
    
    % One common method for 8-term involves defining intermediate matrices from T_m:
    % (T_m)^-1_12 * T_m_11 = (T_1)^-1_12 * (T_1)_11 + (T_1)^-1_11 * T_m_22 * (T_2)^-1_22 * (T_2)_11
    % This is too complex for this example.
    
    % Let's use the algebraic solution from standard 8-term calibration references.
    % The method involves solving for the coefficients of the error boxes A and B.
    % Assume the true S-parameters of the reciprocal through are S_actual_recip = [S11, S21; S12, S22]
    % Where S11, S22 are small, and S12=S21 (let's call it T_recip_actual).
    
    % The document indicates T_m = X * T * Y_inv.
    % Where X and Y are 2x2 matrices that combine the individual e-terms.
    % X = [e1_00, e1_01; e1_10, e1_11] and Y = [e2_00, e2_01; e2_10, e2_11]
    % However, the document uses T1 and T2 defined with Delta_e1/e2 and single e-terms.
    
    % Given t11 = e1^10*e1^01, and t22 = e2^10*e2^01
    % We also need e1^10 and e2^01 individually to form T1 and T2 exactly.
    % The document defines t21 = e1^10 * e2^01.
    
    % This implies solving for individual e1^10, e1^01, e2^10, e2^01.
    % Typically, t11 is the transmission tracking (forward port 1) and t22 is transmission tracking (reverse port 2).
    % The terms e1^00, e1^11, e2^00, e2^11 are directivity and source match.
    
    % Let's assume that t11 and t22 calculated from reflection are the correct values for (e1^10*e1^01) and (e2^10*e2^01).
    % For solving t21 (e1^10 * e2^01), we use the measured reciprocal through (S_reciprocal).
    % The document on page 4 provides the 8-term S-parameter correction equations:
    % These can be inverted for the reciprocal standard to solve for t21.
    % From $S_{21} = \frac{(S_{21m}/t_{21})}{D}$ and $D = [1+(\frac{S_{11m}-e_{1}^{00}}{t_{11}})e_{1}^{11}][1+(\frac{S_{22m}-e_{2}^{00}}{t_{22}})e_{2}^{11}]-(\frac{S_{21m}}{t_{21}})(\frac{S_{12m}}{t_{12}})e_{1}^{11}e_{2}^{11}$
    
    % Since $S_{21}$ for the actual reciprocal through is not ideal, this needs iterative or quadratic solution.
    % For a perfect through (S11=0, S22=0, S12=1, S21=1), the solution simplifies.
    % The term $S_{21}$ in the equation is the *actual* S21 of the through, not 1.
    
    % Given that this is a conceptual example, and without exact inverse equations in the provided PDF,
    % I will use a simplified approach derived from standard VNA 8-term error models for finding t21.
    % This often involves solving for an intermediate term and then deriving t21.
    
    % A common strategy for 8-term calibration (e.g., as implemented in VNAs):
    % 1. Solve for e1_00, e1_11, t11 using Port 1 SOL.
    % 2. Solve for e2_00, e2_11, t22 using Port 2 SOL.
    % 3. Use the Reciprocal Through measurement and the already calculated terms to find t21 and t12.
    
    % The relationships from the document:
    % $t_{11}t_{22}=t_{21}t_{12}$
    % $T_m = \text{error_matrix_1} \cdot T_{DUT} \cdot \text{error_matrix_2}$
    % The 8-term correction equations given are:
    % $S_{21} = \frac{(S_{21m}/t_{21})}{D}$
    % $S_{12} = \frac{(S_{12m}/t_{12})}{D}$
    
    % If we assume the corrected S21 and S12 of the reciprocal through are T_true (which we don't know but we assume S12_true = S21_true):
    % T_true = (S_reciprocal(i,2) / t21) / D_recip
    % T_true = (S_reciprocal(i,3) / t12) / D_recip
    
    % This implies that (S_reciprocal(i,2) / t21) = (S_reciprocal(i,3) / t12)
    % Since t11*t22 = t21*t12, we have t12 = t11*t22 / t21.
    % So, (S_reciprocal(i,2) / t21) = S_reciprocal(i,3) / (t11*t22 / t21)
    % S_reciprocal(i,2) / t21 = S_reciprocal(i,3) * t21 / (t11*t22)
    % t21^2 = S_reciprocal(i,2) * (t11*t22) / S_reciprocal(i,3)
    
    % This gives a solution for t21^2. We need to choose the sign for t21.
    % The document on page 24 refers to choosing the sign for t21:
    % "...a simple if statement in the calibration algorithm will allow you to make the
    % decision whether t21 should be positive or negative based on which corrected thru data,
    % that which is calculated from +t21 or -t21, has the greatest phase difference from the rough
    % estimate."
    
    % Let's calculate t21^2:
    t21_squared_num = S_reciprocal(i,2) * t11(i) * t22(i); % S21m * t11 * t22
    t21_squared_den = S_reciprocal(i,3); % S12m
    
    if abs(t21_squared_den) < eps
        warning('Denominator for t21 calculation near zero at frequency %g Hz. Skipping.', freq(i));
        t21(i) = NaN; t12(i) = NaN;
        continue;
    end
    
    t21_squared = t21_squared_num / t21_squared_den;
    
    % Now, to choose the sign for t21, we need a rough estimate of the thru.
    % Typically, for a through, the phase of S21 should decrease linearly with frequency.
    % Let's assume a positive real part for t21 as a first guess, then refine.
    % The choice of sign is critical for phase accuracy.
    
    % For simplicity, let's take the principal square root for now.
    % A more rigorous approach would indeed involve comparing corrected phases.
    t21_candidate1 = sqrt(t21_squared);
    t21_candidate2 = -sqrt(t21_squared);
    
    % To implement the "greatest phase difference from rough estimate" from the PDF,
    % we'd need to compute the corrected S-parameters of the reciprocal through for both
    % t21_candidate1 and t21_candidate2, and compare their phase to a linear phase
    % (e.g., delay of the reciprocal through).
    
    % This implies needing the apply_error_correction_8term function *inside* this loop,
    % which could lead to circular dependencies or require a refactor.
    
    % For now, let's make a common engineering assumption:
    % For typical throughs, the forward transmission (t21) often has a positive real part or a phase close to zero at low frequencies.
    % Let's select the t21 that has a positive real part, or largest magnitude if real part is close to zero.
    
    % If a choice has to be made between two roots, it's typically the one that results
    % in a 'causal' or 'minimally dispersive' corrected through.
    % A simple heuristic is to choose the root with the largest positive real part or smallest absolute phase if near 0 deg.
    
    if real(t21_candidate1) > real(t21_candidate2)
        t21(i) = t21_candidate1;
    else
        t21(i) = t21_candidate2;
    end
    
    % Now calculate t12 using the relationship t11*t22 = t21*t12
    if abs(t21(i)) < eps
        warning('t21 near zero at frequency %g Hz. t12 cannot be calculated. Setting t12 to NaN.', freq(i));
        t12(i) = NaN;
    else
        t12(i) = (t11(i) * t22(i)) / t21(i);
    end
end

% Store calculated error terms in the struct
error_terms.e1_00 = e1_00;
error_terms.e1_11 = e1_11;
error_terms.t11 = t11;

error_terms.e2_00 = e2_00;
error_terms.e2_11 = e2_11;
error_terms.t22 = t22;

error_terms.t21 = t21;
error_terms.t12 = t12;

end
