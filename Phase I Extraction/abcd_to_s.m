function S = abcd_to_s(T, z0)
%ABCD_TO_S Convert a 2-port ABCD matrix to S-parameters.

A = T(1, 1);
B = T(1, 2);
C = T(2, 1);
D = T(2, 2);

den = A + (B / z0) + (C * z0) + D;
if abs(den) < eps
    error('ABCD-to-S conversion denominator is too small.');
end

S11 = (A + (B / z0) - (C * z0) - D) / den;
S21 = 2 / den;
S12 = 2 * (A * D - B * C) / den;
S22 = (-A + (B / z0) - (C * z0) + D) / den;
S = [S11, S12; S21, S22];
end
