function T = s_to_abcd_local(S, z0)
%S_TO_ABCD_LOCAL Convert a 2-port S matrix to an ABCD matrix.
%   This uses the standard equal-reference-impedance ABCD conversion.

if nargin < 2 || isempty(z0)
    z0 = 50;
end

if abs(S(2, 1)) < eps
    error('S21 is too close to zero for stable S-to-ABCD conversion.');
end

s11 = S(1, 1);
s12 = S(1, 2);
s21 = S(2, 1);
s22 = S(2, 2);

A = ((1 + s11) * (1 - s22) + s12 * s21) / (2 * s21);
B = z0 * ((1 + s11) * (1 + s22) - s12 * s21) / (2 * s21);
C = ((1 - s11) * (1 - s22) - s12 * s21) / (2 * s21 * z0);
D = ((1 - s11) * (1 + s22) + s12 * s21) / (2 * s21);

T = [A, B; C, D];
end
