% script to process data and create figure panels for Calaim*, Dehmelt*, Goncalves* & Machens 2021
% ("Synaptic transmission delays cause uninformed spikes, 
% but networks with high-dimensional inputs are less affected.")
% original file name: figureDelay20200510.m (parts thereof)


clear all
close all

mfilefolder = [fileparts(mfilename('fullpath')),filesep,'..'];  % top-level folder

% Hide partial figures until data processing is complete.
set(0,'DefaultFigureVisible','off');

papercolour = true;
fontsize = 11;

% Process all data from scratch every time the code is run? (Default: true)
liveupdate = true;

%   target = '\\172.25.250.112\fdehmelt\Delay\boxdata12\delay 1.0msec\full\';
target = [mfilefolder,filesep,'data',filesep,'full',filesep];
if liveupdate, CDGM2021_ProcessFile(target), end
[~,dimhandlefull,~] = CDGM2021_PlotData(target,'plot');
%   target = '\\172.25.250.112\fdehmelt\Delay\boxdata12\delay 1.0msec\wide\';
target = [mfilefolder,filesep,'data',filesep,'wide',filesep];
if liveupdate, CDGM2021_ProcessFile(target), end
[~,dimhandlewide,~] = CDGM2021_PlotData(target,'plot');
%   target = '\\172.25.250.112\fdehmelt\Delay\boxdata12\delay 1.0msec\reduced\';
target = [mfilefolder,filesep,'data',filesep,'reduced',filesep];
if liveupdate, CDGM2021_ProcessFile(target), end
[~,dimhandlered,~] = CDGM2021_PlotData(target,'plot');

  ....................................
  
% % %   target = '\\172.25.250.112\fdehmelt\Delay\boxdata12\merged\full\random\';
% %   target = [mfilefolder,filesep,'data',filesep,'vardel',filesep,'full',filesep];
% %   if liveupdate, CDGM2021_ProcessFile(target), end
% %   [~,delhandlefull,~] = CDGM2021_PlotData(target,'plot');
% % %   target = '\\172.25.250.112\fdehmelt\Delay\boxdata12\merged\wide\random\';
% %   target = [mfilefolder,filesep,'data',filesep,'vardel',filesep,'wide',filesep];
% %   if liveupdate, CDGM2021_ProcessFile(target), end
% %   [~,delhandlewide,~] = CDGM2021_PlotData(target,'plot');
% % %   target = '\\172.25.250.112\fdehmelt\Delay\boxdata12\merged\reduced\random\';
% %   target = [mfilefolder,filesep,'data',filesep,'vardel',filesep,'reduced',filesep];
% %   if liveupdate, CDGM2021_ProcessFile(target), end
% %   [~,delhandlered,~] = CDGM2021_PlotData(target,'plot');


% Turn partial figures visible again in order to combine them into one:
% set(0,'DefaultFigureVisible','on');
partialfigurelist = findobj(allchild(groot), 'flat', 'type', 'figure');
for k = 1:numel(partialfigurelist)
  figure(k)
end