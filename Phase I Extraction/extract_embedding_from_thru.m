function embedding = extract_embedding_from_thru(sThru, sLine, options)
%EXTRACT_EMBEDDING_FROM_THRU Extract half-network E from THRU = E * L * E.

if nargin < 3
    options = struct();
end

jumpThreshold = get_opt(options, 'jump_threshold', 0.20);
blendFloor = get_opt(options, 'blend_floor', 0.25);
enableJumpBlend = get_opt(options, 'enable_jump_blend', true);
useSqrtCandidate = get_opt(options, 'use_sqrt_candidate', false);

numFreq = size(sThru, 3);
tEmbedding = zeros(2, 2, numFreq);
sEmbedding = zeros(2, 2, numFreq);

previousRoot = [];

for idx = 1:numFreq
    tThru = s_to_abcd_local(sThru(:, :, idx));
    tLine = s_to_abcd_local(sLine(:, :, idx));
    tHalfSquared = tThru / tLine;
    tRootLog = expm(0.5 * logm(tHalfSquared));
    candidates = cat(3, tRootLog, -tRootLog);
    if useSqrtCandidate
        tRootSqrt = sqrtm(tHalfSquared);
        candidates = cat(3, candidates, tRootSqrt, -tRootSqrt);
    end

    if isempty(previousRoot)
        tRoot = pick_first_root(candidates);
    else
        [tRoot, jumpNorm] = pick_continuous_root(candidates, previousRoot);
        if enableJumpBlend && isfinite(jumpNorm) && jumpNorm > jumpThreshold
            % Soft clamp large isolated jumps to reduce branch artifacts.
            blend = max(blendFloor, jumpThreshold / jumpNorm);
            tRoot = blend * tRoot + (1 - blend) * previousRoot;
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

function val = get_opt(options, fieldName, defaultVal)
if isfield(options, fieldName)
    val = options.(fieldName);
else
    val = defaultVal;
end
end

function tRoot = pick_first_root(candidates)
numCandidates = size(candidates, 3);
scores = inf(numCandidates, 1);
for idx = 1:numCandidates
    cand = candidates(:, :, idx);
    if all(isfinite(cand), 'all')
        scores(idx) = norm(cand - eye(2), 'fro');
    end
end
[~, bestIdx] = min(scores);
tRoot = candidates(:, :, bestIdx);
end

function [tRoot, jumpNorm] = pick_continuous_root(candidates, previousRoot)
numCandidates = size(candidates, 3);
scores = inf(numCandidates, 1);
for idx = 1:numCandidates
    cand = candidates(:, :, idx);
    if all(isfinite(cand), 'all')
        scores(idx) = norm(cand - previousRoot, 'fro');
    end
end
[jumpNorm, bestIdx] = min(scores);
tRoot = candidates(:, :, bestIdx);
end
