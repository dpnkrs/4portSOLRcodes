function [freq, S_param] = read_s_parameter_file(filename)
%READ_S_PARAMETER_FILE Reads S-parameter data from a .s1p or .s2p Touchstone file
%   using MATLAB's RF Toolbox functions (sparameters and rfparam).
%
%   [freq, S_param] = read_s_parameter_file(filename)
%
%   Inputs:
%       filename - Path to the .s1p or .s2p file.
%
%   Outputs:
%       freq     - Column vector of frequencies in Hz.
%       S_param  - Complex matrix of S-parameters.
%                  For .s1p, it's an Nx1 matrix (S11).
%                  For .s2p, it's an Nx4 matrix, with columns ordered as [S11, S21, S12, S22].
%
%   Notes:
%       Requires the RF Toolbox. This method is generally more robust for
%       standard Touchstone files than manual parsing.

try
    % Use the sparameters constructor to read the Touchstone file
    s_params_obj = sparameters(filename);

    % Extract frequencies (sparameters object stores in Hz)
    freq = s_params_obj.Frequencies;

    % Determine number of ports
    num_ports = s_params_obj.NumPorts;

    % Extract S-parameters based on the number of ports
    if num_ports == 1
        % For .s1p, extract S11
        S_param = rfparam(s_params_obj, 1, 1);
    elseif num_ports == 2
        % For .s2p, extract S11, S21, S12, S22
        S11 = rfparam(s_params_obj, 1, 1);
        S21 = rfparam(s_params_obj, 2, 1);
        S12 = rfparam(s_params_obj, 1, 2);
        S22 = rfparam(s_params_obj, 2, 2);
        S_param = [S11, S21, S12, S22];
    else
        % For other port counts, return an error or handle as needed
        error('Unsupported number of ports (%d) in S-parameter file. Only 1-port and 2-port files are supported by this function.', num_ports);
    end

catch ME
    % Catch specific errors related to RF Toolbox or file reading
    if strcmp(ME.identifier, 'MATLAB:UndefinedFunction') && (contains(ME.message, 'sparameters') || contains(ME.message, 'rfparam'))
        error('RF Toolbox not found or required functions (sparameters, rfparam) are undefined. Please ensure you have the RF Toolbox installed and licensed, and that it is added to your MATLAB path.');
    else
        error('Error reading S-parameter file %s: %s', filename, ME.message);
    end
end

end