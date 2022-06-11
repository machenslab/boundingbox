% script to process data and create figure panels for % Calaim*, Dehmelt*, 
% Goncalves* & % Machens 2021 ("Robustness to noise for different signal 
% dimensionalities (...) and different redundancies (...)")


% Process all data from scratch every time the code is run? (Default: true)
liveupdate = true;

% target = '\\172.25.250.112\fdehmelt\Delay\noisedata35\';
target = [mfilefolder,filesep,'data',filesep];
if liveupdate
  CDGM2021_ProcessFile(target);
end
[result,~,data] = CDGM2021_PlotSuppNoiseData(target,'plot');