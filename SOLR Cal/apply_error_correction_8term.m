function S_corrected = apply_error_correction_8term(S_raw, error_terms, freq_idx)
%APPLY_ERROR_CORRECTION_8TERM Applies 8-term error correction to raw S-parameters
%   based on the equations provided in the document "excerptSOLR_MultiportCal.pdf"
%   (Page 4, "8-Term" section).
%
%   S_corrected = apply_error_correction_8term(S_raw, error_terms, freq_idx)
%
%   Inputs:
%       S_raw       - 2x2 complex S-parameter matrix (for a single frequency point).
%                     S_raw = [S11m, S12m; S21m, S22m]
%       error_terms - Struct containing the 8 error terms for ALL frequencies.
%                     Fields: e1_00, e1_11, t11, e2_00, e2_11, t22, t21, t12.
%       freq_idx    - Index of the current frequency point.
%
%   Outputs:
%       S_corrected - 2x2 complex corrected S-parameter matrix.
%                     S_corrected = [S11, S12; S21, S22]

    % Extract error terms for the current frequency from the struct
    e1_00 = error_terms.e1_00(freq_idx);
    e1_11 = error_terms.e1_11(freq_idx);
    t11 = error_terms.t11(freq_idx);
    e2_00 = error_terms.e2_00(freq_idx);
    e2_11 = error_terms.e2_11(freq_idx);
    t22 = error_terms.t22(freq_idx);
    t21 = error_terms.t21(freq_idx);
    t12 = error_terms.t12(freq_idx);

    % Measured S-parameters
    S11m = S_raw(1,1);
    S12m = S_raw(1,2);
    S21m = S_raw(2,1);
    S22m = S_raw(2,2);

    % Calculate the common denominator D as given in the document (Page 4)
    % D = [1+( (S11m-e1^00)/t11 )*e1^11] * [1+( (S22m-e2^00)/t22 )*e2^11] - (S21m/t21)*(S12m/t12)*e1^11*e2^11
    
    % Check for division by zero before calculating D's components
    if abs(t11) < eps || abs(t22) < eps || abs(t21) < eps || abs(t12) < eps
        warning('One or more transmission error terms (t11, t22, t21, t12) are near zero at freq index %d. Correction may be unstable.', freq_idx);
        S_corrected = NaN(2,2);
        return;
    end

    term1 = (S11m - e1_00) / t11;
    term2 = (S22m - e2_00) / t22;
    
    D_factor1 = (1 + term1 * e1_11);
    D_factor2 = (1 + term2 * e2_11);
    D_factor3 = (S21m / t21) * (S12m / t12) * e1_11 * e2_11;
    
    D = D_factor1 * D_factor2 - D_factor3;

    if abs(D) < eps
        warning('Denominator D near zero at frequency index %d. Returning NaNs.', freq_idx);
        S_corrected = NaN(2,2);
        return;
    end
    
    % Calculate corrected S-parameters using the 8-Term equations (Page 4)
    
    % Corrected S11
    S11_num_term1 = ((S11m - e1_00) / t11) * (1 + ((S22m - e2_00) / t22) * e2_11);
    S11_num_term2 = e2_11 * (S21m / t21) * (S12m / t12);
    S11_corr = (S11_num_term1 - S11_num_term2) / D;

    % Corrected S21
    S21_corr = (S21m / t21) / D;

    % Corrected S12
    S12_corr = (S12m / t12) / D;
    
    % Corrected S22 (Equation is slightly ambiguous in PDF, appears as an image.
    % Assuming a symmetric structure to S11:
    % S22 = ( ((S22m-e2^00)/t22) * [1 + ((S11m-e1^00)/t11)*e1^11] - e1^11*(S21m/t21)*(S12m/t12) ) / D
    S22_num_term1 = ((S22m - e2_00) / t22) * (1 + ((S11m - e1_00) / t11) * e1_11);
    S22_num_term2 = e1_11 * (S21m / t21) * (S12m / t12);
    S22_corr = (S22_num_term1 - S22_num_term2) / D;

    S_corrected = [S11_corr, S12_corr; S21_corr, S22_corr];

end
