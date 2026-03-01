function error_terms = calculate_error_terms_solr_PhAmbg(freq, S11_S, S11_O, S11_L, S22_S, S22_O, S22_L, S_reciprocal, thru_delay)
%CALCULATE_ERROR_TERMS_SOLR Calculates the 8 error terms for SOLR calibration
%   with phase ambiguity resolution for t21 using an estimated thru delay.
%
%   Inputs:
%       ... (same as before)
%       thru_delay   - (Optional) Estimated delay of the reciprocal thru in seconds.
%                      Default is 0 if not provided.

if nargin < 9 || isempty(thru_delay)
    thru_delay = 0; 
end

num_freq = length(freq);

% Initialize error terms
e1_00 = zeros(num_freq, 1); t11 = zeros(num_freq, 1); e1_11 = zeros(num_freq, 1);
e2_00 = zeros(num_freq, 1); t22 = zeros(num_freq, 1); e2_11 = zeros(num_freq, 1);
t21 = zeros(num_freq, 1); t12 = zeros(num_freq, 1);

Gamma_S_ideal = -1; Gamma_O_ideal = 1; Gamma_L_ideal = 0;

for i = 1:num_freq
    % --- Reflection Calculations (Port 1) ---
    e1_00(i) = S11_L(i);
    e1_11(i) = ( (S11_S(i) - e1_00(i))/Gamma_S_ideal - (S11_O(i) - e1_00(i))/Gamma_O_ideal ) / ...
               ( (S11_S(i) - e1_00(i)) - (S11_O(i) - e1_00(i)) );
    t11(i) = (S11_S(i) - e1_00(i))/Gamma_S_ideal - e1_11(i) * (S11_S(i) - e1_00(i));

    % --- Reflection Calculations (Port 2) ---
    e2_00(i) = S22_L(i);
    e2_11(i) = ( (S22_S(i) - e2_00(i))/Gamma_S_ideal - (S22_O(i) - e2_00(i))/Gamma_O_ideal ) / ...
               ( (S22_S(i) - e2_00(i)) - (S22_O(i) - e2_00(i)) );
    t22(i) = (S22_S(i) - e2_00(i))/Gamma_S_ideal - e2_11(i) * (S22_S(i) - e2_00(i));

    % --- Transmission tracking (t21) with Phase Ambiguity Resolution ---
    t21_squared = (S_reciprocal(i,2) * t11(i) * t22(i)) / S_reciprocal(i,3);
    
    t21_cand1 = sqrt(t21_squared);
    t21_cand2 = -t21_cand1;

    % RESOLUTION LOGIC:
    % We expect the phase of the corrected Through to be approximately -2*pi*f*delay.
    % In terms of the error term t21, we compare the candidate phase to a rough estimate.
    % Rough estimate: phase(S21m) - phase(expected_thru).
    expected_thru_phase = -2 * pi * freq(i) * thru_delay;
    measured_thru_phase = angle(S_reciprocal(i,2));
    target_t21_phase = measured_thru_phase - expected_thru_phase;

    % Pick the candidate whose phase is closest to the target
    diff1 = abs(angle(exp(1i*(angle(t21_cand1) - target_t21_phase))));
    diff2 = abs(angle(exp(1i*(angle(t21_cand2) - target_t21_phase))));

    if diff1 < diff2
        t21(i) = t21_cand1;
    else
        t21(i) = t21_cand2;
    end
    
    t12(i) = (t11(i) * t22(i)) / t21(i);
end

error_terms.e1_00 = e1_00; error_terms.e1_11 = e1_11; error_terms.t11 = t11;
error_terms.e2_00 = e2_00; error_terms.e2_11 = e2_11; error_terms.t22 = t22;
error_terms.t21 = t21; error_terms.t12 = t12;

end