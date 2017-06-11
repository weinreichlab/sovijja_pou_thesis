function staticAvgFit = arithavgfit(r0, IC50)

staticAvgFit = r0./(1+exp((IC50-(-6.40061))/(-0.6825)));

end

% Input: r0 = allele R00, IC50 = allele IC50
% Output: Malthusian fitness value given r0 and IC50

% This calculates the average fitness of an allele i, given we have its
% IC50, when the drug concentration is static. The value for our static
% drug concentration is -6.40061 = average of the log10(Pyrimethamine
% Concentration). This value is obtained by finding the average value when
% there is no variance (IV Drip), which for all of our alleles, is when the
% growth rate is highest. We use the Avg(log10 (drug concentration) because
% the growth rate formula uses log10(drug concentration). Thus, we don't
% log10(average(drug concentration)). 

% Used to verify our Jensen's Inequality findings. 

