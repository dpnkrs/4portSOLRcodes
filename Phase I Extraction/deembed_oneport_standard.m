function gammaActual = deembed_oneport_standard(gammaMeasured, sEmbedding)
%DEEMBED_ONEPORT_STANDARD Remove one embedding network from a one-port standard.

numFreq = numel(gammaMeasured);
gammaActual = zeros(numFreq, 1);

for idx = 1:numFreq
    s11 = sEmbedding(1, 1, idx);
    s12 = sEmbedding(1, 2, idx);
    s21 = sEmbedding(2, 1, idx);
    s22 = sEmbedding(2, 2, idx);
    deltaGamma = gammaMeasured(idx) - s11;
    denominator = (s12 * s21) + (s22 * deltaGamma);

    if abs(denominator) < eps
        gammaActual(idx) = NaN;
    else
        gammaActual(idx) = deltaGamma / denominator;
    end
end
end
