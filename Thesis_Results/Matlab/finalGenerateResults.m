%Create a new script that creates __ Fitness landscapes and then from those
% fitness landscapes, use those to create ____ geofit runs. 

%


%% Pre-Setup: Initializing/Storing the Fitness Landscapes

% Could only hold ~8000 trials per run due to memory issues. It took around
% 10+ minutes per 8000 trials to create and store the runs. Thus, to run
% more than 100,000 times needed to store multiple variables consolidated
% into one file. The below code was the code used to run and store the
% landscapes used. Un-comment the code below to generate new runs.

% i = 1;
% [~,AFitResults8000] = fit_landscape(8000);
% fprintf('Trial = %d',i) %Prints Trial #
% save('8000RunsGeoFit.mat','AFitResults8000')
% 
% i = i + 1;
% [~,BFitResults8000] = fit_landscape(8000);
% fprintf('Trial = %d',i) %Prints Trial #
% save('8000RunsGeoFit.mat','AFitResults8000','BFitResults8000')
% 
% i = i + 1;
% [~,CFitResults8000] = fit_landscape(8000);
% fprintf('Trial = %d',i) %Prints Trial #
% save('8000RunsGeoFit.mat','AFitResults8000','BFitResults8000','CFitResults8000')
% 
% i = i + 1;
% [~,DFitResults8000] = fit_landscape(8000);
% fprintf('Trial = %d',i) %Prints Trial #
% save('8000RunsGeoFit.mat','AFitResults8000','BFitResults8000','CFitResults8000','DFitResults8000')
% 
% i = i + 1;
% [~,EFitResults8000] = fit_landscape(8000);
% fprintf('Trial = %d',i) %Prints Trial #
% save('8000RunsGeoFit.mat','AFitResults8000','BFitResults8000','CFitResults8000','DFitResults8000','EFitResults8000')
% 
% i = i + 1;
% [~,FFitResults8000] = fit_landscape(8000);
% fprintf('Trial = %d',i) %Prints Trial #
% save('8000RunsGeoFit.mat','AFitResults8000','BFitResults8000','CFitResults8000','DFitResults8000','EFitResults8000','FFitResults8000')
% 
% i = i + 1;
% [~,GFitResults8000] = fit_landscape(8000);
% fprintf('Trial = %d',i) %Prints Trial #
% save('8000RunsGeoFit.mat','AFitResults8000','BFitResults8000','CFitResults8000','DFitResults8000','EFitResults8000','FFitResults8000','GFitResults8000')
% 
% i = i + 1;
% [~,HFitResults8000] = fit_landscape(8000);
% fprintf('Trial = %d',i) %Prints Trial #
% save('8000RunsGeoFit.mat','AFitResults8000','BFitResults8000','CFitResults8000','DFitResults8000','EFitResults8000','FFitResults8000','GFitResults8000','HFitResults8000')
% 
% i = i + 1;
% [~,IFitResults8000] = fit_landscape(8000);
% fprintf('Trial = %d',i) %Prints Trial #
% save('8000RunsGeoFit.mat','AFitResults8000','BFitResults8000','CFitResults8000','DFitResults8000','EFitResults8000','FFitResults8000','GFitResults8000','HFitResults8000','IFitResults8000')
% 
% i = i + 1;
% [~,JFitResults8000] = fit_landscape(8000);
% fprintf('Trial = %d',i) %Prints Trial #
% save('8000RunsGeoFit.mat','AFitResults8000','BFitResults8000','CFitResults8000','DFitResults8000','EFitResults8000','FFitResults8000','GFitResults8000','HFitResults8000','IFitResults8000','JFitResults8000')
% 
% i = i + 1;
% [~,KFitResults8000] = fit_landscape(8000);
% fprintf('Trial = %d',i) %Prints Trial #
% save('8000RunsGeoFit.mat','AFitResults8000','BFitResults8000','CFitResults8000','DFitResults8000','EFitResults8000','FFitResults8000','GFitResults8000','HFitResults8000','IFitResults8000','JFitResults8000','KFitResults8000')
% 
% i = i + 1;
% [~,LFitResults8000] = fit_landscape(8000);
% fprintf('Trial = %d',i) %Prints Trial #
% save('8000RunsGeoFit.mat','AFitResults8000','BFitResults8000','CFitResults8000','DFitResults8000','EFitResults8000','FFitResults8000','GFitResults8000','HFitResults8000','IFitResults8000','JFitResults8000','KFitResults8000','LFitResults8000')
% 
% i = i + 1;
% [~,MFitResults8000] = fit_landscape(8000);
% fprintf('Trial = %d',i) %Prints Trial #
% save('8000RunsGeoFit.mat','AFitResults8000','BFitResults8000','CFitResults8000','DFitResults8000','EFitResults8000','FFitResults8000','GFitResults8000','HFitResults8000','IFitResults8000','JFitResults8000','KFitResults8000','LFitResults8000','MFitResults8000')

%% Determine Average/STDev Fitness for Results

load('8000RunsGeoFit.mat');
%Find average fitness for all alleles over these 100,000 trials

AMean = mean(AFitResults8000,1);
BMean = mean(BFitResults8000,1);
CMean = mean(CFitResults8000,1);
DMean = mean(DFitResults8000,1);
EMean = mean(EFitResults8000,1);
FMean = mean(FFitResults8000,1);
GMean = mean(GFitResults8000,1);
HMean = mean(HFitResults8000,1);
IMean = mean(IFitResults8000,1);
JMean = mean(JFitResults8000,1);
KMean = mean(KFitResults8000,1);
LMean = mean(LFitResults8000,1);
MMean = mean(MFitResults8000,1);

overallMean = (AMean + BMean + CMean + DMean + EMean + FMean + GMean + ...
    HMean + IMean + JMean + KMean + LMean + MMean) ./ 13;

% Consolidates standard deviation into one matrix to calculate overall
% standard deviation. 

fitResultConcat = vertcat(AFitResults8000,BFitResults8000,CFitResults8000, ...
    DFitResults8000,EFitResults8000,FFitResults8000,GFitResults8000, ...
    HFitResults8000,IFitResults8000,JFitResults8000,KFitResults8000, ...
    LFitResults8000,MFitResults8000);

stdReal = std(fitResultConcat)';

%% Determine Peak Information from the Fitness Landscapes
% Because the peak_analysis function also takes a little while to run, we
% calculated the peak information then saved them as variables. We then
% calculate the average of each peak information. peakFreq gives us that 
% when there are n peaks, the fraction of time of how often allele i is a
% peak. 

% [Aa,Ab,~,Ad] = peak_analysis(AFitResults8000);
% [Ba,Bb,~,Bd] = peak_analysis(BFitResults8000);
% [Ca,Cb,~,Cd] = peak_analysis(CFitResults8000);
% [Da,Db,~,Dd] = peak_analysis(DFitResults8000);
% [Ea,Eb,~,Ed] = peak_analysis(EFitResults8000);
% [Fa,Fb,~,Fd] = peak_analysis(FFitResults8000);
% [Ga,Gb,~,Gd] = peak_analysis(GFitResults8000);
% [Ha,Hb,~,Hd] = peak_analysis(HFitResults8000);
% [Ia,Ib,~,Id] = peak_analysis(IFitResults8000);
% [Ja,Jb,~,Jd] = peak_analysis(JFitResults8000);
% [Ka,Kb,~,Kd] = peak_analysis(KFitResults8000);
% [La,Lb,~,Ld] = peak_analysis(LFitResults8000);
% [Ma,Mb,~,Md] = peak_analysis(MFitResults8000);
% 
% 
% peakFreq = (Aa + Ba + Ca + Da + Ea + Fa + Ga + Ha + ...
%      Ia + Ja + Ka + La + Ma) ./ 13;
% peakMean = (Ab + Bb + Cb + Db + Eb + Fb + Gb + Hb + ...
%      Ib + Jb + Kb + Lb + Mb) ./ 13;
%  peakNum = (Ad + Bd + Cd + Dd + Ed + Fd + Gd + Hd + ...
%      Id + Jd + Kd + Ld + Md) ./ 13;

 % save('PeakInfo.mat','peakFreq','peakMean','peakNum')




