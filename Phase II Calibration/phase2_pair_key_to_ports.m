function ports = phase2_pair_key_to_ports(pairKey)
%PHASE2_PAIR_KEY_TO_PORTS Convert P1P4-style key into numeric port indices.

if iscell(pairKey)
    pairKey = pairKey{1};
end
if ~(ischar(pairKey) && isrow(pairKey))
    pairKey = string(pairKey);
    pairKey = pairKey(1);
    pairKey = char(pairKey);
end

tokens = regexp(upper(pairKey), '^P([1-4])P([1-4])$', 'tokens', 'once');
if isempty(tokens)
    error('Invalid pair key: %s', pairKey);
end
ports = [str2double(tokens{1}), str2double(tokens{2})];
end
