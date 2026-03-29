function average_s12_s21_s2p(inputFile, outputFile)
%AVERAGE_S12_S21_S2P Replace S12/S21 in an s2p with their average.
%
% Usage:
%   average_s12_s21_s2p('file.s2p')
%   average_s12_s21_s2p('in.s2p', 'out.s2p')
%
% If outputFile is omitted or empty, the input file is overwritten.

if nargin < 1 || isempty(inputFile)
    error('An input .s2p file is required.');
end

if nargin < 2 || isempty(outputFile)
    outputFile = inputFile;
end

data = read_touchstone_file(inputFile);
S = data.S;

s21 = squeeze(S(2, 1, :));
s12 = squeeze(S(1, 2, :));
sAvg = 0.5 * (s21 + s12);

S(2, 1, :) = reshape(sAvg, 1, 1, []);
S(1, 2, :) = reshape(sAvg, 1, 1, []);

write_touchstone_2port(outputFile, data.freq, S, data.z0);

fprintf('Averaged S12/S21 written to: %s\n', outputFile);
end
