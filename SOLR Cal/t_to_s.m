function S = t_to_s(T)
%T_TO_S Converts 2x2 T-parameter matrix to S-parameter matrix.
%   S = T_TO_S(T)
%   Inputs:
%       T - 2x2 T-parameter matrix.
%   Outputs:
%       S - 2x2 S-parameter matrix.

% Ensure T is 2x2
if ~isequal(size(T), [2, 2])
    error('Input T-parameter matrix must be 2x2.');
end

% S-parameter conversion formula
S = (1/T(2,2)) * [T(1,2), T(1,1)*T(2,2) - T(1,2)*T(2,1); 1, -T(2,1)];

end