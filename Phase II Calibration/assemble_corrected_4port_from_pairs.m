function assembled = assemble_corrected_4port_from_pairs(config, pairBlocks)
%ASSEMBLE_CORRECTED_4PORT_FROM_PAIRS Assemble corrected 4-port DUT from pair submatrices.

freq = pairBlocks(1).freq(:);
nFreq = numel(freq);
S = nan(4, 4, nFreq);
diagAccum = zeros(4, nFreq);
diagWeight = zeros(4, nFreq);

for idxPair = 1:numel(pairBlocks)
    pairBlock = pairBlocks(idxPair);
    ports = pairBlock.ports;
    Spair = pairBlock.S_corrected;
    i = ports(1);
    j = ports(2);

    S(i, j, :) = Spair(1, 2, :);
    S(j, i, :) = Spair(2, 1, :);

    if strcmpi(pairBlock.geometry, 'STRAIGHT')
        weight = config.diagonal_weight_straight;
    else
        weight = config.diagonal_weight_orthogonal;
    end
    diagAccum(i, :) = diagAccum(i, :) + weight * squeeze(Spair(1, 1, :)).';
    diagAccum(j, :) = diagAccum(j, :) + weight * squeeze(Spair(2, 2, :)).';
    diagWeight(i, :) = diagWeight(i, :) + weight;
    diagWeight(j, :) = diagWeight(j, :) + weight;
end

for idxPort = 1:4
    valid = diagWeight(idxPort, :) > 0;
    tmp = nan(1, nFreq);
    tmp(valid) = diagAccum(idxPort, valid) ./ diagWeight(idxPort, valid);
    S(idxPort, idxPort, :) = reshape(tmp, 1, 1, []);
end

assembled = struct();
assembled.freq = freq;
assembled.S = S;
assembled.z0 = config.z0;
end
