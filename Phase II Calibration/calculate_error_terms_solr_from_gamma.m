function errorTerms = calculate_error_terms_solr_from_gamma(freq, measShortP1, measOpenP1, measLoadP1, measShortP2, measOpenP2, measLoadP2, reciprocalMeas, gammaShort, gammaOpen, gammaLoad, thruDelay)
%CALCULATE_ERROR_TERMS_SOLR_FROM_GAMMA Solve 8-term SOLR from measured/reference gammas.

if nargin < 12 || isempty(thruDelay)
    thruDelay = 0;
end

nFreq = numel(freq);
fields = {'e1_00', 'e1_11', 't11', 'e2_00', 'e2_11', 't22', 't21', 't12'};
for idx = 1:numel(fields)
    errorTerms.(fields{idx}) = nan(nFreq, 1);
end

for idx = 1:nFreq
    A1 = [ ...
        1, measShortP1(idx) * gammaShort(idx), gammaShort(idx); ...
        1, measOpenP1(idx) * gammaOpen(idx), gammaOpen(idx); ...
        1, measLoadP1(idx) * gammaLoad(idx), gammaLoad(idx)];
    b1 = [measShortP1(idx); measOpenP1(idx); measLoadP1(idx)];

    A2 = [ ...
        1, measShortP2(idx) * gammaShort(idx), gammaShort(idx); ...
        1, measOpenP2(idx) * gammaOpen(idx), gammaOpen(idx); ...
        1, measLoadP2(idx) * gammaLoad(idx), gammaLoad(idx)];
    b2 = [measShortP2(idx); measOpenP2(idx); measLoadP2(idx)];

    if rcond(A1) < 1e-12 || rcond(A2) < 1e-12
        continue;
    end

    sol1 = A1 \ b1;
    sol2 = A2 \ b2;

    errorTerms.e1_00(idx) = sol1(1);
    errorTerms.e1_11(idx) = sol1(2);
    errorTerms.t11(idx) = sol1(3) + sol1(1) * sol1(2);

    errorTerms.e2_00(idx) = sol2(1);
    errorTerms.e2_11(idx) = sol2(2);
    errorTerms.t22(idx) = sol2(3) + sol2(1) * sol2(2);

    s21m = reciprocalMeas(2, 1, idx);
    s12m = reciprocalMeas(1, 2, idx);
    if abs(s12m) < 1e-15
        continue;
    end

    t21sq = s21m * errorTerms.t11(idx) * errorTerms.t22(idx) / s12m;
    cand1 = sqrt(t21sq);
    cand2 = -cand1;

    targetPhase = angle(s21m) + 2 * pi * freq(idx) * thruDelay;
    diff1 = abs(angle(exp(1i * (angle(cand1) - targetPhase))));
    diff2 = abs(angle(exp(1i * (angle(cand2) - targetPhase))));
    if diff1 <= diff2
        errorTerms.t21(idx) = cand1;
    else
        errorTerms.t21(idx) = cand2;
    end

    if abs(errorTerms.t21(idx)) >= 1e-15
        errorTerms.t12(idx) = (errorTerms.t11(idx) * errorTerms.t22(idx)) / errorTerms.t21(idx);
    end
end
end
