function lineModel = build_uniform_line_2port(freq, lineLengthM, epsEff, zc, alphaNpPerM, z0)
%BUILD_UNIFORM_LINE_2PORT Build an explicit uniform-line 2-port model.

c0 = 299792458;
beta = 2 * pi * freq .* sqrt(epsEff) ./ c0;
gamma = alphaNpPerM + 1i * beta;

numFreq = numel(freq);
T = zeros(2, 2, numFreq);
S = zeros(2, 2, numFreq);

for idx = 1:numFreq
    gl = gamma(idx) * lineLengthM;
    A = cosh(gl);
    B = zc * sinh(gl);
    C = sinh(gl) / zc;
    D = cosh(gl);
    T(:, :, idx) = [A, B; C, D];
    S(:, :, idx) = abcd_to_s(T(:, :, idx), z0);
end

lineModel = struct();
lineModel.freq = freq;
lineModel.length_m = lineLengthM;
lineModel.eps_eff = epsEff;
lineModel.zc = zc;
lineModel.alpha_np_per_m = alphaNpPerM;
lineModel.gamma = gamma;
lineModel.T = T;
lineModel.S = S;
end
