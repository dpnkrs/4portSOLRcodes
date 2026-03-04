function embedding = extract_embedding_from_thru(sThru, sLine)
%EXTRACT_EMBEDDING_FROM_THRU Extract half-network E from THRU = E * L * E.

numFreq = size(sThru, 3);
tEmbedding = zeros(2, 2, numFreq);
sEmbedding = zeros(2, 2, numFreq);

previousRoot = [];

for idx = 1:numFreq
    tThru = s_to_abcd_local(sThru(:, :, idx));
    tLine = s_to_abcd_local(sLine(:, :, idx));
    tHalfSquared = tThru / tLine;
    tRoot = expm(0.5 * logm(tHalfSquared));

    if ~isempty(previousRoot)
        % The matrix half-network remains sign-ambiguous. After the first
        % point, pick the branch that stays closest to the previous sample.
        if norm(previousRoot - tRoot, 'fro') > norm(previousRoot + tRoot, 'fro')
            tRoot = -tRoot;
        end
    end

    tEmbedding(:, :, idx) = tRoot;
    sEmbedding(:, :, idx) = abcd_to_s(tRoot, 50);
    previousRoot = tRoot;
end

embedding = struct();
embedding.T = tEmbedding;
embedding.S = sEmbedding;
end
