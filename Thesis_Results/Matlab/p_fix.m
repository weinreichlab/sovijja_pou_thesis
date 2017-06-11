function probability_fixation = p_fix(s)

probability_fixation = max(1-exp(-2*s),0);
end

% Input: s = selection coefficient
% Output: probability_fixation = probability of fixation given a certain
% s(i->j).

% Calculates probabilitiy of fixation with a given s(i->j).