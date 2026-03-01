function error_terms = calculate_error_terms_solr_NonIdeal_PhAmbg(freq, S11_S_meas, S11_O_meas, S11_L_meas, S22_S_meas, S22_O_meas, S22_L_meas, S_reciprocal_meas, filepath_S_char, filepath_O_char, filepath_L_char, thru_delay)
%CALCULATE_ERROR_TERMS_SOLR Calculates the 8 error terms for SOLR calibration
%   based on 8-term error model principles.
%   This version reads characterized S-parameter data for the Short, Open,
%   and Load standards directly from provided file paths.
%
%   error_terms = calculate_error_terms_solr(...)
%
%   Inputs:
%       freq              - Frequency vector in Hz.
%       S11_S_meas        - Measured S11 of Short standard on Port 1 (complex vector).
%       S11_O_meas        - Measured S11 of Open standard on Port 1 (complex vector).
%       S11_L_meas        - Measured S11 of Load standard on Port 1 (complex vector).
%       S22_S_meas        - Measured S22 of Short standard on Port 2 (complex vector).
%       S22_O_meas        - Measured S22 of Open standard on Port 2 (complex vector).
%       S22_L_meas        - Measured S22 of Load standard on Port 2 (complex vector).
%       S_reciprocal_meas - Measured S-parameters of the Reciprocal Through (Nx4 matrix: [S11_R, S21_R, S12_R, S22_R]).
%       filepath_S_char   - File path to the characterized S11 of Short standard.
%       filepath_O_char   - File path to the characterized S11 of Open standard.
%       filepath_L_char   - File path to the characterized S11 of Load standard.
%
%   Outputs:
%       error_terms  - Struct containing error terms for each frequency point.
%                      Fields correspond to the terms defined in the document:
%                      e1_00, e1_11, t11, e2_00, e2_11, t22, t21, t12.
%                      Each field is a vector of complex numbers.
%
%   Notes:
%       This implementation follows a standard algebraic solution for an 8-term
%       error model, using the provided measured and characterized standard data.
%       The choice of sign for t21 still relies on a heuristic for simplicity;
%       a more rigorous approach might involve comparing corrected thru phases.

if nargin < 12 || isempty(thru_delay)
    thru_delay = 0; 
end

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

% --- Load Characterized Standard Data Internally ---
% This function now takes file paths and loads the data itself.
try
    [~, S11_S_char] = read_s_parameter_file(filepath_S_char);
    [~, S11_O_char] = read_s_parameter_file(filepath_O_char);
    [~, S11_L_char] = read_s_parameter_file(filepath_L_char);
    
    % Ensure characterized data has the same number of frequency points
    if ~isequal(length(S11_S_char), num_freq) || ...
       ~isequal(length(S11_O_char), num_freq) || ...
       ~isequal(length(S11_L_char), num_freq)
        error('Characterized standard files do not have the same number of frequency points as measured data.');
    end
catch ME
    error('Error loading characterized standard files within calculate_error_terms_solr: %s', ME.message);
end
% ---------------------------------------------------

for i = 1:num_freq
    % --- Port 1 Error Term Calculation (e1_00, e1_11, t11) ---
    % Using measured S11 of Short, Open, Load on Port 1 and their characterized S11 values.
    % The relationship is: S_measured = e_00 + (t11 * Gamma_actual) / (1 - e_11 * Gamma_actual)
    %
    % We are now linearizing using: Gamma_m = 1 * e00 + (Gamma_m * Gamma_a) * e11 + Gamma_a * del_e
    % where del_e = t11 - e00 * e11.
    % After solving for del_e, we will calculate t11 = del_e + e00 * e11.
    
    % Check for potential division by zero or ill-conditioned systems
    if abs(S11_S_char(i) - S11_L_char(i)) < eps || abs(S11_O_char(i) - S11_L_char(i)) < eps || ...
       abs(S11_S_char(i) - S11_O_char(i)) < eps
        warning('Characterized standard S11 values are too close at frequency %g Hz. Skipping Port 1 reflection calculation for this point.', freq(i));
        e1_00(i) = NaN; e1_11(i) = NaN; t11(i) = NaN;
        e2_00(i) = NaN; e2_11(i) = NaN; t22(i) = NaN;
        t21(i) = NaN; t12(i) = NaN;
        continue;
    end

    % Intermediate terms for solving the linear system for [e1_00; e1_11; del_e]
    % A_mat * [e1_00; e1_11; del_e] = B_vec
    A_mat = [1, S11_S_meas(i)*S11_S_char(i), S11_S_char(i);
             1, S11_O_meas(i)*S11_O_char(i), S11_O_char(i);
             1, S11_L_meas(i)*S11_L_char(i), S11_L_char(i)];
    
    B_vec = [S11_S_meas(i); S11_O_meas(i); S11_L_meas(i)];

    if rcond(A_mat) < 1e-10 % Check condition number
        warning('Ill-conditioned matrix for Port 1 reflection at frequency %g Hz. Skipping.', freq(i));
        e1_00(i) = NaN; e1_11(i) = NaN; t11(i) = NaN;
    else
        sol = A_mat \ B_vec;
        e1_00(i) = sol(1); % Directivity
        e1_11(i) = sol(2); % Source Match
        del_e = sol(3);    % del_e = t11 - e1_00 * e1_11
        
        % Derive t11 from del_e
        t11(i) = del_e + e1_00(i) * e1_11(i); % Reflection Tracking (e1^10 * e1^01)
    end
    
    % --- Port 2 Error Term Calculation (e2_00, e2_11, t22) ---
    % Similar procedure for Port 2 using S22 of Short, Open, Load on Port 2
    % (Note: S22_S_meas, S22_O_meas, S22_L_meas are the S11 measurements on Port 2)
    
    if abs(S11_S_char(i) - S11_L_char(i)) < eps || abs(S11_O_char(i) - S11_L_char(i)) < eps || ...
       abs(S11_S_char(i) - S11_O_char(i)) < eps
        warning('Characterized standard S11 values are too close at frequency %g Hz. Skipping Port 2 reflection calculation for this point.', freq(i));
        e2_00(i) = NaN; e2_11(i) = NaN; t22(i) = NaN;
        t21(i) = NaN; t12(i) = NaN;
        continue;
    end

    % Intermediate terms for solving the linear system for [e2_00; e2_11; del_e_p2]
    A_mat_p2 = [1, S22_S_meas(i)*S11_S_char(i), S11_S_char(i); % S11_S_char is Gamma_actual
                1, S22_O_meas(i)*S11_O_char(i), S11_O_char(i);
                1, S22_L_meas(i)*S11_L_char(i), S11_L_char(i)];
    
    B_vec_p2 = [S22_S_meas(i); S22_O_meas(i); S22_L_meas(i)];

    if rcond(A_mat_p2) < 1e-10 % Check condition number
        warning('Ill-conditioned matrix for Port 2 reflection at frequency %g Hz. Skipping.', freq(i));
        e2_00(i) = NaN; e2_11(i) = NaN; t22(i) = NaN;
    else
        sol_p2 = A_mat_p2 \ B_vec_p2;
        e2_00(i) = sol_p2(1); % Directivity Port 2
        e2_11(i) = sol_p2(2); % Source Match Port 2
        del_e_p2 = sol_p2(3); % del_e_p2 = t22 - e2_00 * e2_11
        
        % Derive t22 from del_e_p2
        t22(i) = del_e_p2 + e2_00(i) * e2_11(i); % Reflection Tracking Port 2 (e2^10 * e2^01)
    end

    % --- Transmission Error Term Calculation (t21, t12) ---
    % This part uses the measured reciprocal through and the already calculated
    % reflection error terms (e1_00, e1_11, t11, e2_00, e2_11, t22).
    % The derivation for t21 relies on the reciprocal nature of the through standard.
    
    % t21^2 = S21_measured_reciprocal * (t11 * t22) / S12_measured_reciprocal
    t21_squared_num = S_reciprocal_meas(i,2) * t11(i) * t22(i); % S21m_recip * t11 * t22
    t21_squared_den = S_reciprocal_meas(i,3); % S12m_recip
    
    if abs(t21_squared_den) < eps
        warning('Denominator for t21 calculation near zero at frequency %g Hz. Skipping.', freq(i));
        t21(i) = NaN; t12(i) = NaN;
        continue;
    end
    
    t21_squared = t21_squared_num / t21_squared_den;
    
    % Choose the sign for t21. The document suggests choosing based on phase difference
    % from a rough estimate. For simplicity here, we'll choose the root with positive
    % real part or smallest absolute phase. A more rigorous approach might involve
    % correcting the through with both roots and comparing to a linear phase model.
    t21_candidate1 = sqrt(t21_squared);
    t21_candidate2 = -sqrt(t21_squared); % The other root
    
    % RESOLUTION LOGIC:
    % We expect the phase of the corrected Through to be approximately -2*pi*f*delay.
    % In terms of the error term t21, we compare the candidate phase to a rough estimate.
    % Rough estimate: phase(S21m) - phase(expected_thru).
    expected_thru_phase = -2 * pi * freq(i) * thru_delay;
    measured_thru_phase = angle(S_reciprocal_meas(i,2));
    target_t21_phase = measured_thru_phase - expected_thru_phase;

    % Pick the candidate whose phase is closest to the target
    diff1 = abs(angle(exp(1i*(angle(t21_candidate1) - target_t21_phase))));
    diff2 = abs(angle(exp(1i*(angle(t21_candidate2) - target_t21_phase))));

    if diff1 < diff2
        t21(i) = t21_candidate1;
    else
        t21(i) = t21_candidate2;
    end
    
    % Calculate t12 using the relationship t11*t22 = t21*t12
    if abs(t21(i)) < eps
        warning('Calculated t21 is near zero at frequency %g Hz. t12 cannot be calculated. Setting t12 to NaN.', freq(i));
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
