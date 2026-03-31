function [S_corr, diagnostics] = apply_switch_correction_2port(S_raw, gammaForward, gammaReverse, denominatorFloor)
%APPLY_SWITCH_CORRECTION_2PORT Apply pair-specific switch-term correction.

if nargin < 4 || isempty(denominatorFloor)
    denominatorFloor = 1e-9;
end

if isempty(gammaForward) || isempty(gammaReverse)
    S_corr = S_raw;
    diagnostics = struct('denominator', ones(1, size(S_raw, 3)), 'denominator_clamped', false(1, size(S_raw, 3)));
    return;
end

nFreq = size(S_raw, 3);
S_corr = zeros(size(S_raw));
denVec = zeros(nFreq, 1);
clamped = false(nFreq, 1);

for idx = 1:nFreq
    gf = gammaForward(idx);
    gr = gammaReverse(idx);

    s11 = S_raw(1, 1, idx);
    s12 = S_raw(1, 2, idx);
    s21 = S_raw(2, 1, idx);
    s22 = S_raw(2, 2, idx);

    den = 1 - s12 * s21 * gf * gr;
    if abs(den) < denominatorFloor
        if den == 0
            den = denominatorFloor;
        else
            den = denominatorFloor * den / abs(den);
        end
        clamped(idx) = true;
    end

    S_corr(1, 1, idx) = (s11 - s12 * s21 * gf) / den;
    S_corr(2, 2, idx) = (s22 - s12 * s21 * gr) / den;
    S_corr(2, 1, idx) = s21 * (1 - s22 * gf) / den;
    S_corr(1, 2, idx) = s12 * (1 - s11 * gr) / den;
    denVec(idx) = den;
end

diagnostics = struct('denominator', denVec, 'denominator_clamped', clamped);
end
