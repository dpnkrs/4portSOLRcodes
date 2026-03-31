function S_corr = apply_error_correction_8term_local(S_raw, errorTerms, freqIdx, denominatorFloor)
%APPLY_ERROR_CORRECTION_8TERM_LOCAL Apply 8-term correction to one 2x2 matrix.

if nargin < 4 || isempty(denominatorFloor)
    denominatorFloor = 1e-9;
end

e1_00 = errorTerms.e1_00(freqIdx);
e1_11 = errorTerms.e1_11(freqIdx);
t11 = errorTerms.t11(freqIdx);
e2_00 = errorTerms.e2_00(freqIdx);
e2_11 = errorTerms.e2_11(freqIdx);
t22 = errorTerms.t22(freqIdx);
t21 = errorTerms.t21(freqIdx);
t12 = errorTerms.t12(freqIdx);

if any(~isfinite([e1_00, e1_11, t11, e2_00, e2_11, t22, t21, t12])) || ...
        any(abs([t11, t22, t21, t12]) < denominatorFloor)
    S_corr = nan(2, 2);
    return;
end

S11m = S_raw(1, 1);
S12m = S_raw(1, 2);
S21m = S_raw(2, 1);
S22m = S_raw(2, 2);

term1 = (S11m - e1_00) / t11;
term2 = (S22m - e2_00) / t22;
den = (1 + term1 * e1_11) * (1 + term2 * e2_11) - (S21m / t21) * (S12m / t12) * e1_11 * e2_11;
if abs(den) < denominatorFloor
    den = denominatorFloor * exp(1i * angle(den + (den == 0)));
end

S11 = (term1 * (1 + term2 * e2_11) - e2_11 * (S21m / t21) * (S12m / t12)) / den;
S21 = (S21m / t21) / den;
S12 = (S12m / t12) / den;
S22 = (term2 * (1 + term1 * e1_11) - e1_11 * (S21m / t21) * (S12m / t12)) / den;

S_corr = [S11, S12; S21, S22];
end
