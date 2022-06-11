% script to process data and create figure panels for 
% Calaim*, Dehmelt*, Goncalves* & Machens 2021 (e.g., those summarised in 
% the figure "Network response to natural perturbations", panels E and F)


% Process all data from scratch every time the code is run? (Default: true)
liveupdate = true;

% %   target = ''\\172.25.250.112\fdehmelt\Delay\noisedata36\';
% % target = [mfilefolder,filesep,'data',filesep,'full',filesep];
% target = 'C:\Users\dehmelt\Desktop\FD code collection\noisedata36\';
target = [mfilefolder,filesep,'data',filesep];
if liveupdate
  CDGM2021_ProcessFile(target);
end
[result,~,data] = CDGM2021_PlotNoiseData(target,'plot');