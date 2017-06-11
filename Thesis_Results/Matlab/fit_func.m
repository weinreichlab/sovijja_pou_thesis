function malthfit = fit_func(t, r0, IC50)

drug_concentration = log10((11.193*exp(-0.40051*mod(t,168))+ ...
      0.1723*exp(-0.006777*mod(t,168))- ...
      11.364*exp(-0.4146*mod(t,168)))*(1/248710));

malthfit = r0./(1+exp((IC50-drug_concentration)/(-0.6825)));


end

% Input: t = time point of drug cycle schedule. 
%        r0 = R00 of allele
%        IC50 = IC50 of allele
% Output: Fitness of allele at that specific time in drug cycle schedule

% Calculate Malthusian fitness as a function of drug concentration vs. time. 
% Formula from Ogbunugafor et al. 2016.
% Can customize drug_concentration by changing mod(t,n) where n is the
% number of hours per cycle. And also dividing by 168/n. Then multiplying 
% by a constant c.
% New constant c can be calculated by integrating the original
% drug_concentration function (mod 168) and integrating the new function.
% Then deviding the original integral over the new integral to determine
% the ratio in the area under the curve and multiplying by the new
% constant.


