function [fitness_environment,geofit] = fit_landscape(N) 
%N = # environments to simulate/store

muR = [1.398, 1.275, 1.227, 0, 1.37, 1.375, 1.397, 1.219, 1.119, 1.184, ... 
    1.306, 1, 1.273, 1.282, 1.45, 1.250]; %Brown et al 2010 R00 of all alleles
sigmaR = [0.0535, 0.0131, 0.0195, 0, 0.0287, 0.0164, 0.0268, 0.0737, 0.0349, 0.0595, ...
    0.0336, 0.0814, 0.0509, 0.0444, 0.0159, 0.0457]; %Brown et al 2010 R00 standard error of all alleles
muI = [-6.286, -5.812, -4.239, 0, -6.046, -5.774, -3.732, -3.55, -5.724, ...
    -5.491, -4.015, -4.6, -5.773, -5.624, -3.587, -3.3]; %Brown et al 2010 IC50 of all alleles
sigmaI = [0.053, 0.013, 0.014, 0, 0.035, 0.019, 0.025, 0.033, 0.029, ... 
    0.029, 0.017, 0.033, 0.028, 0.034, 0.116, 0.033]; %Brown et al 2010 IC50 standard error of all alleles

num_alleles = length(sigmaR); % Number of alleles
geofit = zeros(N, num_alleles); % Creates N x num_alleles matrix to hold calculated fitness values
fitness_environment = zeros(N,num_alleles,2); %Row = trial number, all alleles in columns, first dimension = R00, second dimension = IC50

for trial_num = 1:N    
    for allele_num=1:num_alleles
        fitness_environment(trial_num,allele_num,1) = normrnd(muR(allele_num),sigmaR(allele_num)); %Fitness Landscape save for R00
        fitness_environment(trial_num,allele_num,2) = normrnd(muI(allele_num),sigmaI(allele_num)); %Fitness Landscape save for IC50
        % Define function to calculate Malthusian Growth rate given R00 
        %and IC50 as a function of x (concentration).
        f = @(x) fit_func(x, fitness_environment(trial_num,allele_num,1), ... 
            fitness_environment(trial_num,allele_num,2));
        
        % Calculates average wrightian fitness and stores in matrix where
        % row = trial, columns = allele.
        geofit(trial_num, allele_num) = exp(integral(f, 0, 2016) / 2016); 
        %Geometric Average Wrightian into matrix where row = trial and columns = allele
    end
end
end

% Creates N fitness environments (samplings of R00 and IC50.. defined as fitness_environment)
% and then calculate each allele's average Wrightian growth rate (geofit). 
% 
% Usually only want geofit -> [~,geofit] = fit_landscape(N)
% fitness_environment = 


% ABSTRACT OUT muR, sigmaR etc and make them arguments.

%save('100000RunsGeoFit.mat','FitResults')

