function T = s_to_t(S)
%S_TO_T Converts 2x2 S-parameter matrix to T-parameter matrix.
%   T = S_TO_T(S)
%   Inputs:
%       S - 2x2 S-parameter matrix.
%   Outputs:
%       T - 2x2 T-parameter matrix.

% Ensure S is 2x2
if ~isequal(size(S), [2, 2])
    error('Input S-parameter matrix must be 2x2.');
end

detS = S(1,1)*S(2,2) - S(1,2)*S(2,1);

% T-parameter conversion formula
T = (1/S(2,1)) * [detS, S(1,1); -S(2,2), 1];

end