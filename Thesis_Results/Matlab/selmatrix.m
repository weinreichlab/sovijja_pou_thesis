function Selection = selmatrix(trial_num,FitResults)
n = 16; %# alleles
Selection = zeros(n,n,trial_num);

for num_trial = 1:trial_num
    Selection(1,[2,3,5,9],num_trial) = sel_coeff(FitResults(num_trial,1), FitResults(num_trial,[2,3,5,9])); %Calculates selection coefficient for each mutational neighbor i->j (including reversions)
    Selection(2,[1,4,6,10],num_trial) = sel_coeff(FitResults(num_trial,2), FitResults(num_trial,[1,4,6,10]));
    Selection(3,[1,4,7,11],num_trial) = sel_coeff(FitResults(num_trial,3), FitResults(num_trial,[1,4,7,11]));
    Selection(4,[2,3,8,12],num_trial) = sel_coeff(FitResults(num_trial,4), FitResults(num_trial,[2,3,8,12]));
    Selection(5,[1,6,7,13],num_trial) = sel_coeff(FitResults(num_trial,5), FitResults(num_trial,[1,6,7,13]));
    Selection(6,[2,5,8,14],num_trial) = sel_coeff(FitResults(num_trial,6), FitResults(num_trial,[2,5,8,14]));
    Selection(7,[3,5,8,15],num_trial) = sel_coeff(FitResults(num_trial,7), FitResults(num_trial,[3,5,8,15]));
    Selection(8,[4,6,7,16],num_trial) = sel_coeff(FitResults(num_trial,8), FitResults(num_trial,[4,6,7,16]));
    Selection(9,[1,10,11,13],num_trial) = sel_coeff(FitResults(num_trial,9), FitResults(num_trial,[1,10,11,13]));
    Selection(10,[2,9,12,14],num_trial) = sel_coeff(FitResults(num_trial,10), FitResults(num_trial,[2,9,12,14]));
    Selection(11,[3,9,12,15],num_trial) = sel_coeff(FitResults(num_trial,11), FitResults(num_trial,[3,9,12,15]));
    Selection(12,[4,10,11,16],num_trial) = sel_coeff(FitResults(num_trial,12), FitResults(num_trial,[4,10,11,16]));
    Selection(13,[5,9,14,15],num_trial) = sel_coeff(FitResults(num_trial,13), FitResults(num_trial,[5,9,14,15]));
    Selection(14,[6,10,13,16],num_trial) = sel_coeff(FitResults(num_trial,14), FitResults(num_trial,[6,10,13,16]));
    Selection(15,[7,11,13,16],num_trial) = sel_coeff(FitResults(num_trial,15), FitResults(num_trial,[7,11,13,16])); 
    Selection(16,[8,12,14,15],num_trial) = sel_coeff(FitResults(num_trial,16), FitResults(num_trial,[8,12,14,15])); 
end
end

% Creates a matrix of selection coefficients. Entry (i,j) = selection
% coefficient of allele i to allele j. Uses sel_coeff function. 
% Input = # of landscapes and fitness matrix (3rd dimension (# of trials) must be same
% size).
