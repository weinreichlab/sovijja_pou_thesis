% Our (incomplete) version of trying to calculate probabilities of
% trajectories. Refer to thesis to look at which future directions this
% script should take in order to calculate accurate probabilities of
% trajectories. Must include indirect trajectories, and it's probably best
% to store these probabilities in an alternate data structure. Because
% right now, we can't really accurately average the probabilities due to
% certain trajectories being available in some landscapes while unavailable
% in others. 

muR = [1.398, 1.275, 1.227, 0, 1.37, 1.375, 1.397, 1.219, 1.119, 1.184, 1.306, 1, 1.273, 1.282, 1.45, 1.250]; %Brandon R0 of all alleles
sigmaR = [0.0535, 0.0131, 0.0195, 0, 0.0287, 0.0164, 0.0268, 0.0737, 0.0349, 0.0595, 0.0336, 0.0814, 0.0509, 0.0444, 0.0159, 0.0457]; %Brandon R0 standard error of all alleles
muI = [-6.286, -5.812, -4.239, 0, -6.046, -5.774, -3.732, -3.55, -5.724, -5.491, -4.015, -4.6, -5.773, -5.624, -3.587, -3.3]; %Brandon IC50 of all alleles (in log)
sigmaI = [0.053, 0.013, 0.014, 0, 0.035, 0.019, 0.025, 0.033, 0.029, 0.029, 0.017, 0.033, 0.028, 0.034, 0.116, 0.033]; %Brandon IC50 standard error of all alleles

load('8000RunsGeoFit.mat');
A = ['0000', '0001', '0010', '0011', '0100', '0101', '0110', '0111', '1000', '1001', '1010', '1011', '1100', '1101', '1110', '1111']; % All alleles
N = 10; % Number of iterations

Ptrajectory1Step = zeros(4,1,N); %Creates matrix to hold trajectory probabilities for peaks of one step
Ptrajectory2Step = zeros(6,2,N); %Creates matrix to hold trajectory probabilities for peaks of two steps
Ptrajectory3Step = zeros(4,6,N); %Creates matrix to hold trajectory probabilities for peaks of three steps
Ptrajectory4Step = zeros(1,12,N); %Creates matrix to hold trajectory probabilities for peaks of four steps

Selection = selmatrix(N,AFitResults8000); % Calculates selection matrix
Pfix = arrayfun(@(x) p_fix(x), Selection); %Calculates pfix for each mutational neighbor (including reversions)
ConditionedPfix = zeros(size(Pfix,1),size(Pfix,2),size(Pfix,3)); %Initializes Conditioned probability matrix

for trial_num = 1:N

    for rownum = 1:size(Pfix,1)
    ConditionedPfix(rownum,:,trial_num) = Pfix(rownum,1:16,trial_num)/ ...
        sum(Pfix(rownum,:,trial_num)); 
    %Calculates pfix conditioned on mutation occurring. The sum of all
    %forward mutations from allele i should equal to 1.
    end
    ConditionedPfix(isnan(ConditionedPfix)) = 0;

    %Ptrajectory1Step(1,1,trial_num) = ConditionedPfix(1,9,trial_num);
    %Ptrajectory1Step(2,1,trial_num) = ConditionedPfix(1,5,trial_num);
    %Ptrajectory1Step(3,1,trial_num) = ConditionedPfix(1,3,trial_num);
    %Ptrajectory1Step(4,1,trial_num) = ConditionedPfix(1,2,trial_num);
    
    %Ptrajectory2Step(1,1,trial_num) = ConditionedPfix(1,9,trial_num)*ConditionedPfix(9,13,trial_num);
    %Ptrajectory2Step(1,2,trial_num) = ConditionedPfix(1,5,trial_num)*ConditionedPfix(5,13,trial_num);
    %Ptrajectory2Step(2,1,trial_num) = ConditionedPfix(1,9,trial_num)*ConditionedPfix(9,11,trial_num);
    %Ptrajectory2Step(2,2,trial_num) = ConditionedPfix(1,3,trial_num)*ConditionedPfix(3,11,trial_num);
    Ptrajectory2Step(3,1,trial_num) = ConditionedPfix(1,2,trial_num)*ConditionedPfix(2,10,trial_num);
    Ptrajectory2Step(3,2,trial_num) = ConditionedPfix(1,9,trial_num)*ConditionedPfix(9,10,trial_num);
    %Ptrajectory2Step(4,1,trial_num) = ConditionedPfix(1,3,trial_num)*ConditionedPfix(3,7,trial_num);
    %Ptrajectory2Step(4,2,trial_num) = ConditionedPfix(1,5,trial_num)*ConditionedPfix(5,7,trial_num);
    %Ptrajectory2Step(5,1,trial_num) = ConditionedPfix(1,2,trial_num)*ConditionedPfix(2,6,trial_num);
    %Ptrajectory2Step(5,2,trial_num) = ConditionedPfix(1,5,trial_num)*ConditionedPfix(5,6,trial_num);
    %Ptrajectory2Step(6,1,trial_num) = ConditionedPfix(1,2,trial_num)*ConditionedPfix(2,4,trial_num);
    %Ptrajectory2Step(6,2,trial_num) = ConditionedPfix(1,3,trial_num)*ConditionedPfix(3,4,trial_num);
    
    Ptrajectory3Step(1,1,trial_num) = ConditionedPfix(1,9,trial_num)*ConditionedPfix(9,13,trial_num)*ConditionedPfix(13,15,trial_num);
    Ptrajectory3Step(1,2,trial_num) = ConditionedPfix(1,9,trial_num)*ConditionedPfix(9,11,trial_num)*ConditionedPfix(11,15,trial_num);
    Ptrajectory3Step(1,3,trial_num) = ConditionedPfix(1,5,trial_num)*ConditionedPfix(5,13,trial_num)*ConditionedPfix(13,15,trial_num);
    Ptrajectory3Step(1,4,trial_num) = ConditionedPfix(1,5,trial_num)*ConditionedPfix(5,7,trial_num)*ConditionedPfix(7,15,trial_num);
    Ptrajectory3Step(1,5,trial_num) = ConditionedPfix(1,3,trial_num)*ConditionedPfix(3,11,trial_num)*ConditionedPfix(11,15,trial_num);
    Ptrajectory3Step(1,6,trial_num) = ConditionedPfix(1,3,trial_num)*ConditionedPfix(3,7,trial_num)*ConditionedPfix(7,15,trial_num);
    %Ptrajectory3Step(2,1,trial_num) = ConditionedPfix(1,9,trial_num)*ConditionedPfix(9,13,trial_num)*ConditionedPfix(13,14,trial_num);
    %Ptrajectory3Step(2,2,trial_num) = ConditionedPfix(1,9,trial_num)*ConditionedPfix(9,10,trial_num)*ConditionedPfix(10,14,trial_num);
    %Ptrajectory3Step(2,3,trial_num) = ConditionedPfix(1,5,trial_num)*ConditionedPfix(5,13,trial_num)*ConditionedPfix(13,14,trial_num);
    %Ptrajectory3Step(2,4,trial_num) = ConditionedPfix(1,5,trial_num)*ConditionedPfix(5,7,trial_num)*ConditionedPfix(7,14,trial_num);
    %Ptrajectory3Step(2,5,trial_num) = ConditionedPfix(1,2,trial_num)*ConditionedPfix(2,10,trial_num)*ConditionedPfix(10,14,trial_num);
    %Ptrajectory3Step(2,6,trial_num) = ConditionedPfix(1,2,trial_num)*ConditionedPfix(2,7,trial_num)*ConditionedPfix(7,14,trial_num);
    %Ptrajectory3Step(3,1,trial_num) = ConditionedPfix(1,9,trial_num)*ConditionedPfix(9,11,trial_num)*ConditionedPfix(11,12,trial_num);
    %Ptrajectory3Step(3,2,trial_num) = ConditionedPfix(1,9,trial_num)*ConditionedPfix(9,10,trial_num)*ConditionedPfix(10,12,trial_num);
    %Ptrajectory3Step(3,3,trial_num) = ConditionedPfix(1,3,trial_num)*ConditionedPfix(3,11,trial_num)*ConditionedPfix(11,12,trial_num);
    %Ptrajectory3Step(3,4,trial_num) = ConditionedPfix(1,3,trial_num)*ConditionedPfix(3,4,trial_num)*ConditionedPfix(4,12,trial_num);
    %Ptrajectory3Step(3,5,trial_num) = ConditionedPfix(1,2,trial_num)*ConditionedPfix(2,10,trial_num)*ConditionedPfix(10,12,trial_num);
    %Ptrajectory3Step(3,6,trial_num) = ConditionedPfix(1,2,trial_num)*ConditionedPfix(2,4,trial_num)*ConditionedPfix(4,12,trial_num);
    %Ptrajectory3Step(4,1,trial_num) = ConditionedPfix(1,5,trial_num)*ConditionedPfix(5,7,trial_num)*ConditionedPfix(7,8,trial_num);
    %Ptrajectory3Step(4,2,trial_num) = ConditionedPfix(1,5,trial_num)*ConditionedPfix(5,6,trial_num)*ConditionedPfix(6,8,trial_num);
    %Ptrajectory3Step(4,3,trial_num) = ConditionedPfix(1,3,trial_num)*ConditionedPfix(3,7,trial_num)*ConditionedPfix(7,8,trial_num);
    %Ptrajectory3Step(4,4,trial_num) = ConditionedPfix(1,3,trial_num)*ConditionedPfix(3,4,trial_num)*ConditionedPfix(4,8,trial_num);
    %Ptrajectory3Step(4,5,trial_num) = ConditionedPfix(1,2,trial_num)*ConditionedPfix(2,6,trial_num)*ConditionedPfix(6,8,trial_num);
    %Ptrajectory3Step(4,6,trial_num) = ConditionedPfix(1,2,trial_num)*ConditionedPfix(2,4,trial_num)*ConditionedPfix(4,8,trial_num);
    %Ptrajectory3Step(1,1:6,trial_num) = ConditionedPfix(1,[9,9,5,5,3,3],trial_num)*ConditionedPfix([9,9,5,5,3,3],[13,11,13,7,11,7],trial_num)*ConditionedPfix([13,11,13,7,11,7],15,trial_num);
    %Ptrajectory3Step(2,1:6,trial_num) = ConditionedPfix(1,[9,9,5,5,2,2],trial_num)*ConditionedPfix([9,9,5,5,2,2],[13,10,13,7,10,7],trial_num)*ConditionedPfix([13,10,13,7,10,7],14,trial_num);    
    %Ptrajectory3Step(3,1:6,trial_num) = ConditionedPfix(1,[9,9,3,3,2,2],trial_num)*ConditionedPfix([9,9,3,3,2,2],[11,10,11,4,10,4],trial_num)*ConditionedPfix([11,10,11,4,10,4],12,trial_num);    
    %Ptrajectory3Step(4,1:6,trial_num) = ConditionedPfix(1,[5,5,3,3,2,2],trial_num)*ConditionedPfix([5,5,3,3,2,2],[7,6,7,4,6,4],trial_num)*ConditionedPfix([7,6,7,4,6,4],8,trial_num);    
        
%    Ptrajectory4Step(1,1,trial_num) = ConditionedPfix(1,9,trial_num)*ConditionedPfix(9,13,trial_num)*ConditionedPfix(13,15,trial_num)*ConditionedPfix(15,16,trial_num); %Calculates trajectory value
%     Ptrajectory4Step(1,2,trial_num) = ConditionedPfix(1,9,trial_num)*ConditionedPfix(9,13,trial_num)*ConditionedPfix(13,14,trial_num)*ConditionedPfix(14,16,trial_num);
%    Ptrajectory4Step(1,3,trial_num) = ConditionedPfix(1,9,trial_num)*ConditionedPfix(9,11,trial_num)*ConditionedPfix(11,15,trial_num)*ConditionedPfix(15,16,trial_num);
%     Ptrajectory4Step(1,4,trial_num) = ConditionedPfix(1,9,trial_num)*ConditionedPfix(9,11,trial_num)*ConditionedPfix(11,12,trial_num)*ConditionedPfix(12,16,trial_num);
%     Ptrajectory4Step(1,5,trial_num) = ConditionedPfix(1,9,trial_num)*ConditionedPfix(9,10,trial_num)*ConditionedPfix(10,14,trial_num)*ConditionedPfix(14,16,trial_num);
%     Ptrajectory4Step(1,6,trial_num) = ConditionedPfix(1,9,trial_num)*ConditionedPfix(9,10,trial_num)*ConditionedPfix(10,12,trial_num)*ConditionedPfix(12,16,trial_num);
%    Ptrajectory4Step(1,7,trial_num) = ConditionedPfix(1,5,trial_num)*ConditionedPfix(5,13,trial_num)*ConditionedPfix(13,15,trial_num)*ConditionedPfix(15,16,trial_num);
%     Ptrajectory4Step(1,8,trial_num) = ConditionedPfix(1,5,trial_num)*ConditionedPfix(5,13,trial_num)*ConditionedPfix(13,14,trial_num)*ConditionedPfix(14,16,trial_num);
%   Ptrajectory4Step(1,9,trial_num) = ConditionedPfix(1,5,trial_num)*ConditionedPfix(5,7,trial_num)*ConditionedPfix(7,15,trial_num)*ConditionedPfix(15,16,trial_num);
%     Ptrajectory4Step(1,10,trial_num) = ConditionedPfix(1,5,trial_num)*ConditionedPfix(5,7,trial_num)*ConditionedPfix(7,8,trial_num)*ConditionedPfix(8,16,trial_num);
%     Ptrajectory4Step(1,11,trial_num) = ConditionedPfix(1,5,trial_num)*ConditionedPfix(5,6,trial_num)*ConditionedPfix(6,14,trial_num)*ConditionedPfix(14,16,trial_num);
%     Ptrajectory4Step(1,12,trial_num) = ConditionedPfix(1,5,trial_num)*ConditionedPfix(5,6,trial_num)*ConditionedPfix(6,8,trial_num)*ConditionedPfix(8,16,trial_num);
%    Ptrajectory4Step(1,13,trial_num) = ConditionedPfix(1,3,trial_num)*ConditionedPfix(3,11,trial_num)*ConditionedPfix(11,15,trial_num)*ConditionedPfix(15,16,trial_num);
%     Ptrajectory4Step(1,14,trial_num) = ConditionedPfix(1,3,trial_num)*ConditionedPfix(3,11,trial_num)*ConditionedPfix(11,12,trial_num)*ConditionedPfix(12,16,trial_num);
%    Ptrajectory4Step(1,15,trial_num) = ConditionedPfix(1,3,trial_num)*ConditionedPfix(3,7,trial_num)*ConditionedPfix(7,15,trial_num)*ConditionedPfix(15,16,trial_num);
%     Ptrajectory4Step(1,16,trial_num) = ConditionedPfix(1,3,trial_num)*ConditionedPfix(3,7,trial_num)*ConditionedPfix(7,8,trial_num)*ConditionedPfix(8,16,trial_num);
%     Ptrajectory4Step(1,17,trial_num) = ConditionedPfix(1,3,trial_num)*ConditionedPfix(3,4,trial_num)*ConditionedPfix(4,12,trial_num)*ConditionedPfix(12,16,trial_num);
%     Ptrajectory4Step(1,18,trial_num) = ConditionedPfix(1,3,trial_num)*ConditionedPfix(3,4,trial_num)*ConditionedPfix(4,8,trial_num)*ConditionedPfix(8,16,trial_num);
%     Ptrajectory4Step(1,19,trial_num) = ConditionedPfix(1,2,trial_num)*ConditionedPfix(2,10,trial_num)*ConditionedPfix(10,14,trial_num)*ConditionedPfix(14,16,trial_num);
%     Ptrajectory4Step(1,20,trial_num) = ConditionedPfix(1,2,trial_num)*ConditionedPfix(2,10,trial_num)*ConditionedPfix(10,12,trial_num)*ConditionedPfix(12,16,trial_num);
%     Ptrajectory4Step(1,21,trial_num) = ConditionedPfix(1,2,trial_num)*ConditionedPfix(2,6,trial_num)*ConditionedPfix(6,14,trial_num)*ConditionedPfix(14,16,trial_num);
%     Ptrajectory4Step(1,22,trial_num) = ConditionedPfix(1,2,trial_num)*ConditionedPfix(2,6,trial_num)*ConditionedPfix(6,8,trial_num)*ConditionedPfix(8,16,trial_num);
%     Ptrajectory4Step(1,23,trial_num) = ConditionedPfix(1,2,trial_num)*ConditionedPfix(2,4,trial_num)*ConditionedPfix(4,12,trial_num)*ConditionedPfix(12,16,trial_num);
%     Ptrajectory4Step(1,24,trial_num) = ConditionedPfix(1,2,trial_num)*ConditionedPfix(2,4,trial_num)*ConditionedPfix(4,8,trial_num)*ConditionedPfix(8,16,trial_num);

end






