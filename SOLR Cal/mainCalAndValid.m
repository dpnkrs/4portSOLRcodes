% main_solr_simulation.m
% This script simulates SOLR calibration and applies it.

clear; clc; close all;

%% 1. Define File Paths for Measured Calibration Standards
% You need to provide your actual file paths here.
% Ensure these files exist and are in .s1p or .s2p Touchstone format.

% Measured 1-port standards (S11 only)
% Assume you measured Short, Open, Load on Port 1, then Short, Open, Load on Port 2.
% If your VNA saves S11 and S22 in a single .s2p for a one-port measurement,
% adjust the `read_s_parameter_file` function or the loading part here.

% Example file names (replace with your actual data)
filepath_short_p1 = 'data/short.s1p'; % Measured S11 of short on Port 1
filepath_open_p1 = 'data/open.s1p';   % Measured S11 of open on Port 1
filepath_load_p1 = 'data/load.s1p';   % Measured S11 of load on Port 1

filepath_short_p2 = 'data/short.s1p'; % Measured S22 of short on Port 2
filepath_open_p2 = 'data/open.s1p';   % Measured S22 of open on Port 2
filepath_load_p2 = 'data/load.s1p';   % Measured S22 of load on Port 2

filepath_reciprocal_thru = 'data/thruOarc_dc170.s2p'; % Measured S11, S21, S12, S22 of reciprocal standard

filepath_dut_raw = 'data/thruOdiag_dc170.s2p'; % Raw measurement of a DUT to be corrected

% Corrected standards
filepath_S_char = 'data/tug_mtrl/Short-tugmtrl.s2p'; % Measured S11 of short on Port 1
filepath_O_char = 'data/tug_mtrl/Open-tugmtrl.s2p';   % Measured S11 of open on Port 1
filepath_L_char = 'data/tug_mtrl/Load-tugmtrl.s2p';   % Measured S11 of load on Port 1
% filepath_S_char = 'data/solr_thruOarc/solr_thruOarcShort.s2p'; % Measured S11 of short on Port 1
% filepath_O_char = 'data/solr_thruOarc/solr_thruOarcOpen.s2p';   % Measured S11 of open on Port 1
% filepath_L_char = 'data/solr_thruOarc/solr_thruOarcLoad.s2p';   % Measured S11 of load on Port 1


%% 2. Load Measured S-Parameter Data
disp('Loading measured calibration standard data...');
[freq_s_p1, S11_S_p1] = read_s_parameter_file(filepath_short_p1);
[~, S11_O_p1] = read_s_parameter_file(filepath_open_p1);
[~, S11_L_p1] = read_s_parameter_file(filepath_load_p1);

[~, S22_S_p2] = read_s_parameter_file(filepath_short_p2); % Assuming S22_S_p2 is from a 1-port measurement on P2
[~, S22_O_p2] = read_s_parameter_file(filepath_open_p2);
[~, S22_L_p2] = read_s_parameter_file(filepath_load_p2);

[freq_rec_thru, S_reciprocal] = read_s_parameter_file(filepath_reciprocal_thru);

% Check if frequencies match for all standards
if ~isequal(freq_s_p1, freq_rec_thru)
    error('Frequencies of calibration standards do not match. Ensure all measurements are at the same frequency points.');
end
freq = freq_s_p1; % Use this as the common frequency vector

% Ensure S22 for Port 2 standards are correct (e.g., if .s1p they will be S11, so map to S22)
% If S22_S_p2, S22_O_p2, S22_L_p2 were loaded as S11 from .s1p, they need to be treated as S22.
% The `read_s_parameter_file` assumes S11 for .s1p, so this is fine.

fprintf('Data loaded for %d frequency points.\n', length(freq));

%% 3. Calculate Estimated Thru Delay for SOLR
% recirpocal thru DUT length is approx 329.2 um, epsilon effective is 3.7245
L_thru = 315e-6;%329.2e-6; 
eps_eff_thru = 3.7245;
c0 = physconst('LightSpeed');
thru_delay_est = (L_thru * sqrt(eps_eff_thru)) / c0; 
fprintf('Estimated Thru Delay: %.4f ps\n', thru_delay_est * 1e12);

%% 4. Calculate Error Terms (SOLR Calibration)
disp('Calculating 8-term error coefficients...');
% error_terms = calculate_error_terms_solr_NonIdeal(freq, ...
%                                           S11_S_p1, S11_O_p1, S11_L_p1, ...
%                                           S22_S_p2, S22_O_p2, S22_L_p2, ...
%                                           S_reciprocal, filepath_S_char, ...
%                                           filepath_O_char, filepath_L_char);
% error_terms = calculate_error_terms_solr(freq, ...
%                                           S11_S_p1, S11_O_p1, S11_L_p1, ...
%                                           S22_S_p2, S22_O_p2, S22_L_p2, ...
%                                           S_reciprocal);
% error_terms = calculate_error_terms_solr_PhAmbg(freq, S11_S_p1, S11_O_p1, S11_L_p1, ...
%                                          S22_S_p2, S22_O_p2, S22_L_p2, ...
%                                          S_reciprocal, thru_delay_est);

error_terms = calculate_error_terms_solr_NonIdeal_PhAmbg(freq, ...
                                          S11_S_p1, S11_O_p1, S11_L_p1, ...
                                          S22_S_p2, S22_O_p2, S22_L_p2, ...
                                          S_reciprocal, filepath_S_char, ...
                                          filepath_O_char, filepath_L_char, thru_delay_est);

disp('Error coefficients calculated.');

%% 5. Apply Calibration to a Raw DUT Measurement
disp('Applying calibration to raw DUT data...');
[freq_dut, S_dut_raw] = read_s_parameter_file(filepath_dut_raw);

if ~isequal(freq_dut, freq)
    error('DUT measurement frequencies do not match calibration frequencies.');
end

S_dut_corrected = zeros(size(S_dut_raw)); % Initialize corrected S-parameters
for i = 1:length(freq)
    S_raw_matrix = [S_dut_raw(i,1), S_dut_raw(i,3); S_dut_raw(i,2), S_dut_raw(i,4)]; % [S11, S12; S21, S22]
    S_dut_corrected_matrix = apply_error_correction_8term(S_raw_matrix, error_terms, i);
    S_dut_corrected(i,:) = [S_dut_corrected_matrix(1,1), S_dut_corrected_matrix(2,1), S_dut_corrected_matrix(1,2), S_dut_corrected_matrix(2,2)]; % [S11, S21, S12, S22]
end
disp('DUT data corrected.');

%% 6. Validation: Correct the Calibration Standards Themselves
% This is a crucial step to check the calibration's accuracy.
% Corrected S-parameters of the standards should match their ideal values.

disp('Validating calibration by correcting calibration standards...');

% Ideal Standard Values (adjust for your specific cal kit and definition)
% For perfect standards:
S_ideal_S = -1; % Magnitude 1, Angle -180 deg
S_ideal_O = 1;  % Magnitude 1, Angle 0 deg
S_ideal_L = 0;  % Magnitude 0, Angle 0 deg (perfect 50 Ohm)
% Ideal Reciprocal Through (assuming negligible loss, for a first check):
% This is the "true" S-param of the reciprocal standard (e.g., a short section of 50-ohm line).
% For a typical SOLR calibration, the "true" S-parameters of the reciprocal through are NOT necessarily ideal (e.g., S11, S22 are small but non-zero, S12, S21 are close to 1 but with phase/loss).
% For validation, you'd ideally have the *measured* true S-parameters of your reciprocal standard.
% For simplicity in this example, let's assume it's a very good through.
S_ideal_Reciprocal_Thru = [0, 1; 1, 0]; % [S11, S12; S21, S22] for ideal through

S_corrected_S_p1 = zeros(length(freq), 1);
S_corrected_O_p1 = zeros(length(freq), 1);
S_corrected_L_p1 = zeros(length(freq), 1);
S_corrected_S_p2 = zeros(length(freq), 1);
S_corrected_O_p2 = zeros(length(freq), 1);
S_corrected_L_p2 = zeros(length(freq), 1);
S_corrected_Reciprocal_Thru = zeros(size(S_reciprocal));

for i = 1:length(freq)
    % For 1-port standards, we can effectively treat them as a 2-port with open circuited opposite port
    % Or, more precisely, use only the S11 or S22 correction if available (not in 8-term which is 2-port).
    % To apply 8-term correction to 1-port measurements, you'd typically make them "pseudo 2-port"
    % by assuming the other port is an ideal match. This is a simplification.
    
    % A simpler way for validation of reflection terms:
    % Use the measured S11_S_p1, and manually apply the 1-port error model implied by the 8-term model.
    % The error_terms are designed for a 2-port DUT.
    
    % For validation, we should take the measured standards and "correct" them using the derived error terms.
    % This will show how well the calibration works.

    % Correcting the measured Short standard (Port 1):
    % Assume S22_meas=0, S12_meas=0, S21_meas=0 for a 1-port reflection measurement
    S_meas_short_p1 = [S11_S_p1(i), 0; 0, 0]; % Pseudo 2-port for 1-port short measurement
    S_corrected_short_p1_matrix = apply_error_correction_8term(S_meas_short_p1, error_terms, i);
    S_corrected_S_p1(i) = S_corrected_short_p1_matrix(1,1);
    
    % Correcting the measured Open standard (Port 1):
    S_meas_open_p1 = [S11_O_p1(i), 0; 0, 0];
    S_corrected_open_p1_matrix = apply_error_correction_8term(S_meas_open_p1, error_terms, i);
    S_corrected_O_p1(i) = S_corrected_open_p1_matrix(1,1);
    
    % Correcting the measured Load standard (Port 1):
    S_meas_load_p1 = [S11_L_p1(i), 0; 0, 0];
    S_corrected_load_p1_matrix = apply_error_correction_8term(S_meas_load_p1, error_terms, i);
    S_corrected_L_p1(i) = S_corrected_load_p1_matrix(1,1);

    % Correcting the measured Short standard (Port 2):
    S_meas_short_p2 = [0, 0; 0, S22_S_p2(i)]; % Pseudo 2-port for 1-port short measurement on Port 2
    S_corrected_short_p2_matrix = apply_error_correction_8term(S_meas_short_p2, error_terms, i);
    S_corrected_S_p2(i) = S_corrected_short_p2_matrix(2,2); % Corrected S22

    % Correcting the measured Open standard (Port 2):
    S_meas_open_p2 = [0, 0; 0, S22_O_p2(i)];
    S_corrected_open_p2_matrix = apply_error_correction_8term(S_meas_open_p2, error_terms, i);
    S_corrected_O_p2(i) = S_corrected_open_p2_matrix(2,2);

    % Correcting the measured Load standard (Port 2):
    S_meas_load_p2 = [0, 0; 0, S22_L_p2(i)];
    S_corrected_load_p2_matrix = apply_error_correction_8term(S_meas_load_p2, error_terms, i);
    S_corrected_L_p2(i) = S_corrected_load_p2_matrix(2,2);

    % Correcting the measured Reciprocal Through:
    S_meas_rec_thru = [S_reciprocal(i,1), S_reciprocal(i,3); S_reciprocal(i,2), S_reciprocal(i,4)];
    S_corrected_rec_thru_matrix = apply_error_correction_8term(S_meas_rec_thru, error_terms, i);
    S_corrected_Reciprocal_Thru(i,:) = [S_corrected_rec_thru_matrix(1,1), S_corrected_rec_thru_matrix(2,1), S_corrected_rec_thru_matrix(1,2), S_corrected_rec_thru_matrix(2,2)];
end

disp('Calibration standards corrected for validation.');

%% 7. Plotting Results
disp('Plotting results...');

figure;
subplot(2,2,1);
plot(freq/1e9, 20*log10(abs(S_dut_raw(:,1))), 'r--', 'DisplayName', 'Raw S11'); hold on;
plot(freq/1e9, 20*log10(abs(S_dut_corrected(:,1))), 'b-', 'DisplayName', 'Corrected S11');
title('DUT S11 (Magnitude)'); xlabel('Frequency (GHz)'); ylabel('Magnitude (dB)'); legend; grid on;

subplot(2,2,2);
% plot(freq/1e9, angle(S_dut_raw(:,1))*180/pi, 'r--', 'DisplayName', 'Raw S11'); hold on;
% plot(freq/1e9, angle(S_dut_corrected(:,1))*180/pi, 'b-', 'DisplayName', 'Corrected S11');
% title('DUT S11 (Phase)'); xlabel('Frequency (GHz)'); ylabel('Phase (degrees)'); legend; grid on;
plot(freq/1e9, rad2deg(unwrap(angle(S_dut_raw(:,1)))), 'r--'); hold on;
plot(freq/1e9, rad2deg(unwrap(angle(S_dut_corrected(:,1)))), 'b-');
title('DUT S11 Phase (Unwrapped)'); xlabel('Frequency (GHz)'); ylabel('Degrees'); grid on;

subplot(2,2,3);
plot(freq/1e9, 20*log10(abs(S_dut_raw(:,2))), 'r--', 'DisplayName', 'Raw S21'); hold on;
plot(freq/1e9, 20*log10(abs(S_dut_corrected(:,2))), 'b-', 'DisplayName', 'Corrected S21');
title('DUT S21 (Magnitude)'); xlabel('Frequency (GHz)'); ylabel('Magnitude (dB)'); legend; grid on;

subplot(2,2,4);
% plot(freq/1e9, angle(S_dut_raw(:,2))*180/pi, 'r--', 'DisplayName', 'Raw S21'); hold on;
% plot(freq/1e9, angle(S_dut_corrected(:,2))*180/pi, 'b-', 'DisplayName', 'Corrected S21');
% title('DUT S21 (Phase)'); xlabel('Frequency (GHz)'); ylabel('Phase (degrees)'); legend; grid on;
plot(freq/1e9, rad2deg(angle(S_dut_raw(:,2))), 'r--', 'DisplayName', 'Raw'); hold on;
plot(freq/1e9, rad2deg(angle(S_dut_corrected(:,2))), 'b-', 'DisplayName', 'Corrected');
title('DUT S21 Phase (Unwrapped)'); xlabel('Frequency (GHz)'); ylabel('Degrees'); grid on;

sgtitle('Raw vs. Corrected DUT S-Parameters');

% Plotting Validation Results for Standards
figure;
subplot(2,3,1);
plot(freq/1e9, 20*log10(abs(S_corrected_S_p1)), 'b-', 'DisplayName', 'Corrected Short S11'); hold on;
plot(freq/1e9, 20*log10(abs(S_ideal_S * ones(size(freq)))), 'r--', 'DisplayName', 'Ideal Short S11 (Mag)');
title('Corrected Short S11 (Port 1)'); xlabel('Frequency (GHz)'); ylabel('Magnitude (dB)'); legend; grid on;
ylim([-0.5, 0.5]); % For ideal Short, it should be 0 dB.

subplot(2,3,2);
plot(freq/1e9, 20*log10(abs(S_corrected_O_p1)), 'b-', 'DisplayName', 'Corrected Open S11'); hold on;
plot(freq/1e9, 20*log10(abs(S_ideal_O * ones(size(freq)))), 'r--', 'DisplayName', 'Ideal Open S11 (Mag)');
title('Corrected Open S11 (Port 1)'); xlabel('Frequency (GHz)'); ylabel('Magnitude (dB)'); legend; grid on;
ylim([-0.5, 0.5]);

subplot(2,3,3);
plot(freq/1e9, 20*log10(abs(S_corrected_L_p1)), 'b-', 'DisplayName', 'Corrected Load S11'); hold on;
plot(freq/1e9, 20*log10(abs(S_ideal_L * ones(size(freq)))), 'r--', 'DisplayName', 'Ideal Load S11 (Mag)');
title('Corrected Load S11 (Port 1)'); xlabel('Frequency (GHz)'); ylabel('Magnitude (dB)'); legend; grid on;
ylim([-40, 0]); % For ideal Load, it should be -inf dB (showing residual error)

subplot(2,3,4);
plot(freq/1e9, 20*log10(abs(S_corrected_Reciprocal_Thru(:,2))), 'b-', 'DisplayName', 'Corrected Recip Thru S21'); hold on;
plot(freq/1e9, 20*log10(abs(S_ideal_Reciprocal_Thru(2,1) * ones(size(freq)))), 'r--', 'DisplayName', 'Ideal Thru S21 (Mag)');
title('Corrected Reciprocal Thru S21'); xlabel('Frequency (GHz)'); ylabel('Magnitude (dB)'); legend; grid on;
ylim([-1, 1]); % For ideal Thru, it should be 0 dB.

subplot(2,3,5);
plot(freq/1e9, 20*log10(abs(S_corrected_Reciprocal_Thru(:,1))), 'b-', 'DisplayName', 'Corrected Recip Thru S11'); hold on;
plot(freq/1e9, 20*log10(abs(S_ideal_Reciprocal_Thru(1,1) * ones(size(freq)))), 'r--', 'DisplayName', 'Ideal Thru S11 (Mag)');
title('Corrected Reciprocal Thru S11'); xlabel('Frequency (GHz)'); ylabel('Magnitude (dB)'); legend; grid on;
ylim([-40, 0]); % For ideal Thru, S11 should be -inf dB.

subplot(2,3,6);
plot(freq/1e9, 20*log10(abs(S_corrected_Reciprocal_Thru(:,4))), 'b-', 'DisplayName', 'Corrected Recip Thru S22'); hold on;
plot(freq/1e9, 20*log10(abs(S_ideal_Reciprocal_Thru(2,2) * ones(size(freq)))), 'r--', 'DisplayName', 'Ideal Thru S22 (Mag)');
title('Corrected Reciprocal Thru S22'); xlabel('Frequency (GHz)'); ylabel('Magnitude (dB)'); legend; grid on;
ylim([-40, 0]);

sgtitle('Calibration Validation: Corrected Standards vs. Ideal Values');


%% 8. Determine Valid Frequency Range (Based on Residual Error)

disp('Analyzing residual errors to determine valid frequency range...');

% Define acceptable error thresholds (tune these based on your application)
max_reflection_error_dB = -10; % e.g., anything worse than -30dB return loss is poor
max_transmission_ripple_dB = 0.5; % e.g., +- 0.1 dB variation from ideal 0 dB
max_phase_ripple_deg = 5; % e.g., +- 5 degrees variation from ideal 0 degrees

% Residual errors for reflection standards (Port 1)
residual_S_p1 = 20*log10(abs(S_corrected_S_p1)) - 20*log10(abs(S_ideal_S));
residual_O_p1 = 20*log10(abs(S_corrected_O_p1)) - 20*log10(abs(S_ideal_O));
residual_L_p1 = 20*log10(abs(S_corrected_L_p1)) - 20*log10(abs(S_ideal_L)); % Note: log10(0) is -Inf, handle this.
residual_L_p1(isinf(residual_L_p1)) = -100; % Replace -Inf with a large negative number for plotting

% For Load, it's better to look at the corrected magnitude directly and see if it's below a threshold
corrected_load_magnitude_dB = 20*log10(abs(S_corrected_L_p1));
corrected_load_magnitude_dB(isinf(corrected_load_magnitude_dB)) = -100; % Handle -Inf

% Residual errors for transmission standard
residual_S21_rec_thru = abs(20*log10(abs(S_corrected_Reciprocal_Thru(:,2))) - 20*log10(abs(S_ideal_Reciprocal_Thru(2,1))));
residual_S12_rec_thru = abs(20*log10(abs(S_corrected_Reciprocal_Thru(:,3))) - 20*log10(abs(S_ideal_Reciprocal_Thru(1,2))));
residual_S11_rec_thru_dB = 20*log10(abs(S_corrected_Reciprocal_Thru(:,1))); % Should be very low
residual_S22_rec_thru_dB = 20*log10(abs(S_corrected_Reciprocal_Thru(:,4))); % Should be very low


% Determine valid frequency range
is_reflection_valid = (corrected_load_magnitude_dB < max_reflection_error_dB) & ...
                      (abs(20*log10(abs(S_corrected_S_p1))) < abs(max_reflection_error_dB)) & ... % Short should be near 0 dB mag
                      (abs(20*log10(abs(S_corrected_O_p1))) < abs(max_reflection_error_dB)); % Open should be near 0 dB mag

is_transmission_valid = (residual_S21_rec_thru < max_transmission_ripple_dB) & ...
                        (residual_S12_rec_thru < max_transmission_ripple_dB) & ...
                        (residual_S11_rec_thru_dB < max_reflection_error_dB) & ... % S11 of thru
                        (residual_S22_rec_thru_dB < max_reflection_error_dB); % S22 of thru

% Combine for overall validity
is_valid_freq = is_reflection_valid & is_transmission_valid;

valid_freq_indices = find(is_valid_freq);

if isempty(valid_freq_indices)
    fprintf('No valid frequency range found within the specified error thresholds.\n');
else
    min_valid_freq = freq(valid_freq_indices(1));
    max_valid_freq = freq(valid_freq_indices(end));
    fprintf('Estimated valid frequency range: %.3f GHz to %.3f GHz\n', min_valid_freq/1e9, max_valid_freq/1e9);
end

% Plotting Residual Errors (for more detailed analysis)
figure;
subplot(3,1,1);
plot(freq/1e9, corrected_load_magnitude_dB, 'b', 'DisplayName', 'Corrected Load S11 (Mag)'); hold on;
plot(freq/1e9, max_reflection_error_dB * ones(size(freq)), 'r--', 'DisplayName', 'Max Acceptable Error');
title('Corrected Load Return Loss (Port 1)'); xlabel('Frequency (GHz)'); ylabel('Magnitude (dB)'); legend; grid on;

subplot(3,1,2);
plot(freq/1e9, residual_S21_rec_thru, 'b', 'DisplayName', 'Reciprocal Thru S21 Residual Error'); hold on;
plot(freq/1e9, max_transmission_ripple_dB * ones(size(freq)), 'r--', 'DisplayName', 'Max Acceptable Ripple');
title('Reciprocal Through S21 Magnitude Residual Error'); xlabel('Frequency (GHz)'); ylabel('Error (dB)'); legend; grid on;

subplot(3,1,3);
plot(freq/1e9, residual_S11_rec_thru_dB, 'b', 'DisplayName', 'Reciprocal Thru S11 (Corrected)'); hold on;
plot(freq/1e9, max_reflection_error_dB * ones(size(freq)), 'r--', 'DisplayName', 'Max Acceptable S11');
title('Reciprocal Through S11 (Corrected)'); xlabel('Frequency (GHz)'); ylabel('Magnitude (dB)'); legend; grid on;

sgtitle('Calibration Validity Assessment');

%% Paper results
figure;
subplot(1,2,1);
plot(freq/1e9, 20*log10(abs(S_dut_raw(:,2))), 'r--', 'DisplayName', 'Raw-Sim S21'); hold on;
plot(freq/1e9, 20*log10(abs(S_dut_corrected(:,2))), 'b-', 'DisplayName', 'Corrected S21');
title('DUT S21 (Magnitude)'); xlabel('Frequency (GHz)'); ylabel('Magnitude (dB)'); legend; grid on;

subplot(1,2,2);
plot(freq/1e9, angle(S_dut_raw(:,2))*180/pi, 'r--', 'DisplayName', 'Raw-Sim S21'); hold on;
plot(freq/1e9, angle(S_dut_corrected(:,2))*180/pi, 'b-', 'DisplayName', 'Corrected S21');
title('DUT S21 (Phase)'); xlabel('Frequency (GHz)'); ylabel('Phase (degrees)'); legend; grid on;

%sgtitle('Raw vs. Corrected DUT S-Parameters');