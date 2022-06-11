% code to reproduce a figure of Calaim*, Dehmelt*, Goncalves* & Machens 2021
% ("Synaptic transmission delays cause uninformed spikes, 
%  but networks with high-dimensional inputs are less affected.")

clear all
close all

trialnumber = 10;  % number of trials for each combination of dimensionality, redundancy)
% Note that because of the various sources of variability (input, decoder tuning, etc), 
% meaningful results generally require a large number of trials.

simulate = true;  % Simulate new trials from scratch...
analyse = true;   % ...and/or process existing data files?

% Each simulated trial is saved to a separate *.mat file in the appropriate subfolder
% of the "data" folder ("full" for default box, full connectivity networks; 
% "wide" for wide box networks; "reduced" for the reduced connectivity network).
% Once processed, these results are combined into a single file per folder,
% datasummary.mat. Each time data is processed, this file is overwritten.

mfilefolder = fileparts(mfilename('fullpath'));
if simulate
  
  list.dimensionality = [20];
  list.redundancy = [2 5 10 20];
  delay = 1e-3;  % in seconds, so "1e-3" is 1msec
  
  run([mfilefolder,filesep,'subfunctions',filesep,'CDGM2021_SimulateDelay.m'])
  
end

if analyse
  
  run([mfilefolder,filesep,'subfunctions',filesep,'CDGM2021_ProcessDelayData.m'])
  run([mfilefolder,filesep,'subfunctions',filesep,'CDGM2021_MakeFigure.m'])
  
end