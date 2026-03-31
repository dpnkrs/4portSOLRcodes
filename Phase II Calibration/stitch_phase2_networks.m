function stitched = stitch_phase2_networks(networks, interpMethod)
%STITCH_PHASE2_NETWORKS Stitch bandwise networks onto one continuous grid.

if nargin < 2 || isempty(interpMethod)
    interpMethod = 'pchip';
end

allFreq = [];
for idx = 1:numel(networks)
    allFreq = [allFreq; networks{idx}.freq(:)]; %#ok<AGROW>
end
freq = unique(allFreq);

nPorts = size(networks{1}.S, 1);
nFreq = numel(freq);
S = nan(nPorts, nPorts, nFreq);
weights = zeros(1, nFreq);

for idxNet = 1:numel(networks)
    net = networks{idxNet};
    mask = freq >= min(net.freq) & freq <= max(net.freq);
    if ~any(mask)
        continue;
    end
    for row = 1:nPorts
        for col = 1:nPorts
            vals = squeeze(net.S(row, col, :));
            realInterp = interp1(net.freq, real(vals), freq(mask), interpMethod, 'extrap');
            imagInterp = interp1(net.freq, imag(vals), freq(mask), interpMethod, 'extrap');
            interpVals = complex(realInterp, imagInterp);
            current = squeeze(S(row, col, mask));
            if all(isnan(current))
                S(row, col, mask) = reshape(interpVals, 1, 1, []);
            else
                merged = current;
                nanMask = isnan(merged);
                merged(nanMask) = 0;
                S(row, col, mask) = reshape(merged + interpVals, 1, 1, []);
            end
        end
    end
    weights(mask) = weights(mask) + 1;
end

for row = 1:nPorts
    for col = 1:nPorts
        vals = squeeze(S(row, col, :)).';
        valid = weights > 0;
        vals(valid) = vals(valid) ./ weights(valid);
        S(row, col, :) = reshape(vals, 1, 1, []);
    end
end

stitched = struct();
stitched.freq = freq(:);
stitched.S = S;
stitched.z0 = networks{1}.z0;
end
