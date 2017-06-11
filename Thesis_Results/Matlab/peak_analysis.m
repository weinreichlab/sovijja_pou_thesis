function [AllelePeakFrequencyFraction,AllelePeakFraction,PeaksMean,NumPeaks] = peak_analysis(FitResults)

N = size(FitResults,1);
num_alleles = 16;
Neighbors = [2,3,5,9;1,4,6,10;1,4,7,11;2,3,8,12;...
    1,6,7,13;2,5,8,14;3,5,8,15;4,6,7,16;1,10,11,13;2,9,12,14;...
    3,9,12,15;4,10,11,16;5,9,14,15;6,10,13,16;7,11,13,16;8,12,14,15];
    % All mutationally adjacent neighbors for each allele. 
Peaks = zeros(16,8,N); %Rows = Alleles, Col = # of peaks sum, N = trial_num
PeaksIntermediate = zeros(16,N,N); %Holds whether allele is peak or not in each trial
PeakCounter = zeros(1,N); %Counts total number of peaks per trial
AllelePeakRecorder = zeros(16,1); %Records # of times an allele is a peak
NumPeaks = zeros(1,8); %Records # of times there are i peaks

for trial_num = 1:N
    for allele_num=1:num_alleles
        if FitResults(trial_num,allele_num) > max(FitResults(trial_num,Neighbors(allele_num,:))); %If allele_num fitness > all mutant neighbors, classify as peak
            PeaksIntermediate(allele_num,trial_num,trial_num) = PeaksIntermediate(allele_num,trial_num,trial_num) + 1; %For every allele, tallies everytime it is a peak, 3 Dimensions to correlate to 3 Dimensions of Peaks
            PeakCounter(trial_num) = PeakCounter(trial_num) + 1; %Tallies total number of peaks per trial
            AllelePeakRecorder(allele_num) = AllelePeakRecorder(allele_num) + 1; %Tallies total number of times an allele is a peak
        end
    end
    NumPeaks(PeakCounter(trial_num)) = NumPeaks(PeakCounter(trial_num)) + 1;
    Peaks(:,PeakCounter(trial_num),trial_num) = PeaksIntermediate(:,trial_num,trial_num);    
    AllelePeakFrequency = sum(Peaks,3); %When there are n peaks, how often is allele i a peak
end

NumPeaks = NumPeaks ./ N;
AllelePeakFraction = AllelePeakRecorder(:)./sum(AllelePeakRecorder(:)); %Fraction of times in total that allele i is a peak

AllelePeakFrequencyFraction = AllelePeakFrequency; %Gives fraction of time that given n peaks, allele i is a peak
PeakColumnSums = sum(AllelePeakFrequencyFraction); %Sums columns
for i = 1:numel(PeakColumnSums)
    AllelePeakFrequencyFraction(:,i) = AllelePeakFrequencyFraction(:,i)./PeakColumnSums(i); %Divides each column by sum of column
end
PeaksMean = mean(Peaks,3);

end


% Input: % FitResults = Input Geofit Results matrix
% Output:% AllelePeakFrequencyFraction = Fraction of time that given an
         %   environment has column j peaks, that allele i is a peak
         % AllelePeakFraction = Fraction of time in total that allele i = peak
         % PeaksMean = Mean of peak occurences across all trials.
         % NumPeaks = # of times there are i peaks
         
% This function outputs various statistics of the peaks of the fitness
% landscape environment. This is specific for the 16 allele environment
% that is combinatorially complete. 