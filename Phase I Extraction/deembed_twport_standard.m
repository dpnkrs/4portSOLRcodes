function sActual = deembed_twport_standard(sMeasured, sEmbedding)
%DEEMBED_TWPORT_STANDARD Remove E from both sides of a measured 2-port.

numFreq = size(sMeasured, 3);
sActual = zeros(2, 2, numFreq);

for idx = 1:numFreq
    tMeasured = s_to_abcd_local(sMeasured(:, :, idx));
    tEmbedding = s_to_abcd_local(sEmbedding(:, :, idx));
    tActual = tEmbedding \ tMeasured / tEmbedding;
    sActual(:, :, idx) = abcd_to_s(tActual, 50);
end
end
