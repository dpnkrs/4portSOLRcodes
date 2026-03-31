function S2 = phase2_extract_pair_submatrix(S4, pairKey)
%PHASE2_EXTRACT_PAIR_SUBMATRIX Extract active-pair 2x2 submatrix from a 4-port array.

ports = phase2_pair_key_to_ports(pairKey);
S2 = S4(ports, ports, :);
end
