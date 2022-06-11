% code to reproduce a figure of Calaim*, Dehmelt*, Goncalves* & Machens 2021
% ("Networks response to natural perturbations (...) noise level (...)")

clear all
close all

trialnumber = 100;  % number of trials for each combination of voltage noise, redundancy.
% Note that because of the various sources of variability (input, decoder tuning, etc), 
% meaningful results generally require a large number of trials.

simulate = true;  % Simulate new trials from scratch...
analyse = true;   % ...and/or process existing data files?

% Each simulated trial is saved to a separate *.mat file in the "data" folder.
% Once processed, these results are combined into a single file per folder,
% datasummary.mat. Each time data is processed, this file is overwritten.

mfilefolder = fileparts(mfilename('fullpath'));
if simulate
  
  list.dimensionality = 2;
  list.redundancynarrow = [3 10 50];  % redundancy of narrow-box networks (T=0.5)
  list.redundancywide = 50;           % redundancy of wide-box networks (T=0.7)
  list.noisestd = [0 .02 .04 .08 .16 .32 .64 1.25 2 3];
  
  run([mfilefolder,filesep,'subfunctions',filesep,'CDGM2021_SimulateNoise.m'])
  
end

if analyse
  
  run([mfilefolder,filesep,'subfunctions',filesep,'CDGM2021_ProcessNoiseData.m'])
  run([mfilefolder,filesep,'subfunctions',filesep,'CDGM2021_MakeNoiseFigure.m'])
  
end