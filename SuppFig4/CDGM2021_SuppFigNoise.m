% code to reproduce a supplementary figure of Calaim*, Dehmelt*, Goncalves* & 
% Machens 2021 ("Robustness to noise for different signal dimensionalities (...)
% and different redundancies (...)")

clear all
close all

trialnumber = 4;  % number of trials for each combination of voltage noise, redundancy.
% Note that because of the various sources of variability (input, decoder tuning, etc), 
% meaningful results generally require a large number of trials.

simulate = true;  % Simulate new trials from scratch...
analyse = true;   % ...and/or process existing data files?

% Each simulated trial is saved to a separate *.mat file in the "data" folder.
% Once processed, these results are combined into a single file per folder,
% datasummary.mat. Each time data is processed, this file is overwritten.


mfilefolder = fileparts(mfilename('fullpath'));
if simulate
  
  list.dimensionality = [5 10 20 50];
  list.redundancynarrow = [3 10 50];  % redundancy of narrow-box networks (T=0.5)
  list.redundancywide = 50;           % redundancy of wide-box networks (T=0.7)
  list.noisestd = [0 .02 .04 .08 .16 .32 .64 1.25 2 3];
  
%   % alternative set for incomplete, but faster simulations:
%   list.dimensionality = [50];
%   list.noisestd = [0 .32 3];

  
  run([mfilefolder,filesep,'subfunctions',filesep,'CDGM2021_SimulateNoise.m'])
  
end

if analyse
  
  run([mfilefolder,filesep,'subfunctions',filesep,'CDGM2021_ProcessSuppNoiseData.m'])
  run([mfilefolder,filesep,'subfunctions',filesep,'CDGM2021_MakeSuppNoiseFigure.m'])
  
end