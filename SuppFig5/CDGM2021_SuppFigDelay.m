% code to reproduce a Supplementary Figure of Calaim*, Dehmelt*, Goncalves & Machens 2021
% ("Single trials of delayed and undelayed networks for intermediate dimensionalities")
% original file name: exampleFigure20200808_delay.m

clear all
close all

dimension = 20;
redundancy = 5;
delay = 1e-3;

seedstring = num2str(cputime);  % same seed for all cases
% seedstring = [];                % different seed for each case

% figures 1/2: fully connected default box, undelayed/delayed
CDGM2021_SimulateSCN(dimension,redundancy,0,'full',seedstring)
% figures 3/4: fully connected default box, undelayed/delayed
CDGM2021_SimulateSCN(dimension,redundancy,delay,'full',seedstring)
% figures 5/6: fully connected adaptive box, undelayed/delayed
CDGM2021_SimulateSCN(dimension,redundancy,delay,'wide',seedstring)
% figures 7/8: reduced excitation default box, undelayed/delayed
CDGM2021_SimulateSCN(dimension,redundancy,delay,'reduced',seedstring)
% % figures 9/10: inhibition-only default box, undelayed/delayed
% CalaimDehmeltMachens_scnSimulation(dimension,redundancy,delay,'none',seedstring)

% figures 11, 12, 13, 14, 15: excerpt
figureCircularDelay20200210(1,2)
figureCircularDelay20200210(3,4)
figureCircularDelay20200210(5,6)
figureCircularDelay20200210(7,8)
% figureCircularDelay20200210(9,10)

% figure 11: merger
figureMergeExample20200411(9,10,11,12)
updateMergeExample20200411(13)
% figureMergeExample20200411(11,12,13,14,15)
% updateMergeExample20200411(16)

close(1:12)

% remove higher-dimensional readouts and targets to improve visibility:
fg = figure(13);
ax = fg.Children([end-12,end-8,end-4,end]);
for k = 1:numel(ax)
  ob = ax(k).Children;
  [ob(6:end-4).Visible] = deal('off');
  if dimension > 3
    [ob(4:5).Color] = deal(ob(ceil(end*1/3)).Color);
    [ob(end-3:end-2).Color] = deal(ob(ceil(end*2/3)).Color);
  end
end
