function [par,his,hon,performance,networkerror,pro,option] = scnDelay(dimensionchoice,coveragechoice,perturbationchoice,networkchoice,varargin)

  % PREAMBLE 1: Initialise RNG and disable unnecessary warnings.

  rng('default')
  if numel(varargin) == 1 && ~isempty(varargin{1})
    seedstring = varargin{1};
  else
    seedstring = num2str(cputime);
  end
  seedstring(seedstring=='.') = [];
  rng(str2double(seedstring))

  warning('off','MATLAB:interp1:UsePCHIP')
  % The line above deactivates the following cbrewer warning:
  % Warning: INTERP1(...,'CUBIC') will change in a future release. Use INTERP1(...,'PCHIP') instead. 
  % > In interp1>sanitycheckmethod (line 265)
  %   In interp1>parseinputs (line 401)
  %   In interp1 (line 80)
  %   In interpolate_cbrewer (line 31)
  %   In cbrewer (line 101)
  %   In heatmap20180823 (line 21)
  
  if gpuDeviceCount > 0
    warning('off','parallel:gpu:device:DeviceLibsNeedsRecompiling')
  end
  % The lines above deactivate the following gpu-related warning:
  % Warning: The CUDA driver must recompile the GPU libraries because your device is   
  % more recent than the libraries. Recompiling can take several minutes.
  % > In Delay20190630 (line 19)
  %   In run (line 91)
  
  
  
  % PREAMBLE 2: Choose which features to include or exclude.
  
  
  boxchoicelist = readOptimalBoxSize;                   % look up table of numerical estimates
  excitationchoicelist = readOptimalExcitationRemoved;  % look up table of numerical estimates
% %   boxchoicelist = [10, 10, 3e-4, .055]
%   manualsize = .55
%   boxchoicelist = [dimensionchoice, coveragechoice, perturbationchoice, manualsize];
%   boxchoicelist = [20, 5, 1e-3, 1.00];
% %   excitationchoicelist = [50, 2, 1e-3, .7]
  
  % FD20200423: increase box size by another 10%, just to be safe
  boxchoicelist(:,4) = boxchoicelist(:,4) * 1.1;
    
  option.network = networkchoice;
  
  switch option.network
    
    case 'full'
      boxchoice = 0.55;
      option.excitation = true;
      option.lowerbound = false;
      option.removestrongexcitation = false; % only applied to test network, not reference network,
      % which is always full connectivity now - and always has box size 0.55 here!      

    case 'wide'
      
      if numel(dimensionchoice) == 1 && numel(coveragechoice) == 1
        row = find(boxchoicelist(:,1) == dimensionchoice & ...
                   boxchoicelist(:,2) == coveragechoice & ...
                   boxchoicelist(:,3) == perturbationchoice);
       if isempty(row)
%          error('No box size was predefined for the given combination of parameters.')
         warning(['No box size was predefined for the given combination of parameters. ', ...
                  'Using default box size instead.'])
         boxchoice = 0.55;
       else
         boxchoice = boxchoicelist(row,4);
       end
      else
        error('This specific version can only test one type of network at a time.')
      end
  
      option.excitation = true;
      option.lowerbound = false;
      option.removestrongexcitation = false; % only applied to test network, not reference network,
      % which is always full connectivity now - and always has box size 0.55 here!      
      
    case 'narrownoise'
      
%       option.network = 'full';  % continue from here on as in 'full' case (e.g., normal recurrence)
      boxchoice = 0.5;
%       boxchoice = 0.55
      
      option.excitation = true;
      option.lowerbound = false;
      option.removestrongexcitation = false; % only applied to test network, not reference network,
      % which is always full connectivity now - and always has box size 0.55 here!
      
    case 'widenoise'
      
%       option.network = 'full';  % continue from here on as in 'full' case (e.g., normal recurrence)
      boxchoice = 0.7;
      
      option.excitation = true;
      option.lowerbound = false;
      option.removestrongexcitation = false; % only applied to test network, not reference network,
      % which is always full connectivity now - and always has box size 0.55 here!      
      
    case 'reduced'
      boxchoice = 0.55;
      option.excitation = true;
      option.lowerbound = false;
      option.removestrongexcitation = true; % only applied to test network, not reference network,
      % which is always full connectivity now - and always has box size 0.55 here!
      
      if numel(dimensionchoice) == 1 && numel(coveragechoice) == 1
        row = find(excitationchoicelist(:,1) == dimensionchoice & ...
                   excitationchoicelist(:,2) == coveragechoice & ...
                   excitationchoicelist(:,3) == perturbationchoice);
       if isempty(row)
         warning(['No fraction of excitation to be cut was predefined for the given ', ...
                  'combination of parameters. Using default removal fraction instead.'])
         option.fractionremoved = 0.8;
       else
         option.fractionremoved = excitationchoicelist(row,4);
       end
      else
        error('This specific version can only test one type of network at a time.')
      end
      
    case 'none'
      boxchoice = 0.55;
      option.excitation = false;
      option.lowerbound = true;
%       option.lowerbound = false;
      option.removestrongexcitation = false; % only applied to test network, not reference network,
      % which is always full connectivity now - and always has box size 0.55 here!      
      
  end
    
  option.savetofile = true;
  option.showplot = true;
  
  
  % choose decoder type:
  
%   option.weight = 'random';
%   option.weight = 'random plus cardinal';
  option.weight = 'equidistant';
%   option.weight = 'manually equidistant in 2D';


  % choose input type:
  
%   option.input = 'frozen';
  option.input = 'sine 2D';
%   option.input = 'sine ND';
%   option.input = 'random constant';
%   option.input = 'smooth ramp normalised';
%   option.input = 'smooth ramp interpolated';


  % choose refractory period:
  
%   option.refraction = .002;
  option.refraction = 2e-4;
%   option.refraction = []

  

  % PREAMBLE 3: Read desired parameter ranges.
  
  dimensionlist = dimensionchoice;
  coveragelist = coveragechoice;
  perturbationlist = perturbationchoice;
  polytopestretch = 2*boxchoice;
  numrepeat = 1;

  perturbationlist = perturbationlist(randperm(numel(perturbationlist))); % more balanced sampling in case of crashes
  coveragelist = coveragelist(randperm(numel(coveragelist))); % more balanced sampling in case of crashes
  dimensionlist = dimensionlist(randperm(numel(dimensionlist))); % more balanced sampling in case of crashes
  
  
  
  % PREAMBLE 4: Assemble all possible combinations of parameters.

  parameterset = NaN(3,numel(dimensionlist),numel(coveragelist),numel(perturbationlist));
  for dimensionindex  = 1:numel(dimensionlist)
    for coverageindex = 1:numel(coveragelist)
      for perturbationindex  = 1:numel(perturbationlist)
        parameterset(:,dimensionindex,coverageindex,perturbationindex) = ...
          [dimensionlist(dimensionindex); coveragelist(coverageindex); perturbationlist(perturbationindex)];
      end
    end
  end

  setlist = reshape(parameterset,[3 numel(parameterset)/3]);
  
  
  
  % THE ACTUAL SIMULATIONS:

  format shortg
  
%   [par,his,hon,honReference,ran] = deal([]); % To avoid errors in parfor
%   % because workers cannot access "the file containing ran", because to 
%   % them, it looks like a function.
  
  for setindex = 1:size(setlist,2)

    tic

    numdimension = setlist(1,setindex);
    coverage = setlist(2,setindex);
    numneuron = coverage*numdimension;    

    perturbation = setlist(3,setindex);
    
    for repeatindex = 1:numrepeat

      if repeatindex > 1   % pick a new seed after each trial
        rng('shuffle');
      end
      randomseed = randi(1e5);

      
      % Pick new decoding weights before each trial (note: in return, frozen input):
      switch option.weight
        case 'random'
          rng(randomseed)
          feedforward = randn(numneuron,numdimension);
          feedforward = feedforward ./ repmat(sqrt(sum(feedforward.^2,2)), [1 numdimension]);
          
        case 'random plus cardinal'
          feedforward(1:numdimension,:) = diag(ones(numdimension,1));
          if numneuron >= 2*numdimension
            feedforward(numdimension+1:2*numdimension,:) = -diag(ones(numdimension,1));
          end
          
        case 'equidistant'
          strength = .01;
          steps = 1e3;
          [weightcell,~,~] = findWeight(numdimension,numneuron,steps,strength);
          feedforward = weightcell{1} / 1; % <--- manual scaling
          
        case 'manually equidistant in 2D'
          % manual equidistance plus noise in 2D
          increment = 360/numneuron;
          offset = -72;
          range = mod(offset+(increment:increment:360), 360);
%           randn('seed',1234567)
          range = range + [6*randn([1 numel(range)-1]) 0];
          feedforward = [sind(range)', +cosd(range)'];

        otherwise
          error(['''option.weight'' must be one of the following: ''random'', ''random ', ...
                 'plus cardinal'', ''equidistant'', ''manually equidistant in 2D'''])
      
      end
      
      
      [~,~,honReference(setindex,repeatindex),~,~] = ...
        runSCN(option, feedforward, 1.1, randomseed, 0);  % reference always has box size 0.55

      [par(setindex,repeatindex), his(setindex,repeatindex), hon(setindex,repeatindex),~,pro] = ...
        runSCN(option, feedforward, polytopestretch, randomseed, perturbation);

      [performance(setindex,repeatindex),networkerror(setindex,repeatindex)] = ...
        NFPCperformance(hon(setindex,repeatindex).target, ...
                        hon(setindex,repeatindex).readout, ...
                        honReference(setindex,repeatindex).readout);

      iterationscompleted = setindex;
      iterationsintotal = size(setlist,2);
      disp([num2str(iterationscompleted),' of ',num2str(iterationsintotal),' completed.'])
      
      if option.savetofile
        saveSCN(par(setindex,repeatindex), his(setindex,repeatindex), ...
                hon(setindex,repeatindex), honReference(setindex,repeatindex), ...
                performance(setindex,repeatindex), networkerror(setindex,repeatindex), ...
                option)
      end
      
      
%       % HOUSEKEEPING: optionally, reduce active memory usage (after saving data)
%
%       parfield = fieldnames(par);
%       hisfield = fieldnames(his);
%       honfield = fieldnames(hon);
%       hrffield = fieldnames(honReference);
%       for fieldindex = 1:numel(parfield)
%         par(setindex,repeatindex).(parfield{fieldindex}) = [];
%       end
%       for fieldindex = 1:numel(hisfield)
%         his(setindex,repeatindex).(hisfield{fieldindex}) = [];
%       end
%       for fieldindex = 1:numel(honfield)
%         hon(setindex,repeatindex).(honfield{fieldindex}) = [];
%       end
%       for fieldindex = 1:numel(hrffield)
%         honReference(setindex,repeatindex).(honfield{fieldindex}) = [];
%       end
      
    end
  end
end



function [inputweight,lasttrial,laststrength] = findWeight(NDrange,NCrange,NS,strength)
     
  inputweight = cell(numel(NDrange),numel(NCrange));   % input weights chosen
  lasttrial    = NaN(numel(NDrange),numel(NCrange));   % no. of trials needed
  laststrength = NaN(numel(NDrange),numel(NCrange));   % repulsion strength needed
   
  display(' ')
   
  for j = 1:numel(NDrange)
    for k = 1:numel(NCrange)
       
      fprintf('\b\b')
      display(['Now computing set ',num2str(k+(j-1)*numel(NCrange)), ...
         ' of ',num2str(numel(NCrange)*numel(NDrange)),':'])
 
      ND = NDrange(j);
      NC = NCrange(k);
       
      if ND == 1
         
        positive = +ones(1,ceil(NC/2));
        negative = -ones(1,floor(NC/2));
        weightsample = NaN(NC,1);
        weightsample(1:2:end) = positive;
        weightsample(2:2:end) = negative;
         
      else
         
  %         weightsample = weightShell(NC*2^ND,ND,strength);
        [weightsample,pos,rep] = weightShell(NC,ND,NS,strength);
%         weightShellDraw(pos,rep)
         
      end
 
      % output all weights
      inputweight{j,k} = weightsample;
      display(['--- Done for ',num2str(NC),' cells, ',num2str(ND),' dimensions. ---'])
%       display(' ')
       
%       fprintf('\n')
%       display(repmat(' ',[1 14]))
                  
    end
  end
   
  display('--- All done. ---')
  fprintf('\n')
  display(repmat(' ',[1 14+20]))
   
end
 
 
  
function [finalpos,pos,rep] = weightShell(NC,ND,NS,initialstrength)
 
%   NS = 1e4;
   
  pos = NaN(NC,ND,NS);
%   pos(:,:,1) = ones(NC,ND,1) + 1e-1*rand(NC,ND,1);
  pos(:,:,1) = rand(NC,ND,1) - 0*.5;  % important change: centre around zero
  pos(:,:,1) = pos(:,:,1)./repmat(sqrt(sum(pos(:,:,1).^2,2)),[1 ND 1]);
  rep = NaN(NC,ND,NS);
   
  strength = initialstrength;   % initialise repulsion strength
  flipcounter = 1;              % initialise counter for changes to strength
  convergence = 0;              % initialise convergence flag as false
  checkpoint = 35;       % when to check whether variance has begun to decrease noticeably
%   checkpointConv = 20;          % when to check whether variance has begun to approach zero
  checkrange = 5;               % how many steps each are used to compute early/later variance
   
  % pos = position of weights over time, in NC x ND x NS
  % distance = distance vectors between all pairs of points, in NC x ND x NC
  % euclidean = abs. value of distances betw. pairs, in NC x ND x NC (for further computation)
  % geodesic = inner product between all pairs of points, NC x ND x NC
  % direction = unit distance vectors, i.e. normalised "distance", NC x ND x NC
  % repulsion = sum of distance vectors, each weighted by euclidean, in NC x ND
  % definition of euclidean calls "min" to eliminate "Inf" contribution from self-pairing
 
%   display(strength)
 
  attempt = 0;
  while ~convergence
     
    for t = 1:NS
 
  %       pos(:,:,t) = max(pos(:,:,t),1e-6*ones(NC,ND));
 
      % (2) randomize slightly, then normalize length to force weights onto unit hypersphere
      pos(:,:,t) = pos(:,:,t) + 0e-6*abs(rand(NC,ND,1));
      pos(:,:,t) = pos(:,:,t) ./ repmat(sqrt(sum(pos(:,:,t).^2,2)),[1 ND]);
 
      % (3b) compute distance measures
      distance  = repmat(permute(pos(:,:,t),[3 2 1]),[NC 1 1]) - repmat(pos(:,:,t),[1 1 NC]);
      euclidean = repmat(sqrt(sum(distance.^2,2)),[1 ND 1]);
      geodesic  = repmat(permute(real(acos(pos(:,:,t)*pos(:,:,t)')),[1 3 2]),[1 ND 1]);
      direction = distance ./ repmat(sqrt(sum(distance.^2,2)),[1 ND 1]);
 
      % (3b) sanitize distance measures: minimally positive distance between a point and itself
      [euclidean(euclidean==0), geodesic(geodesic==0)] = deal(1e10*rand(1));
      direction(isnan(direction)) = deal(0);
 
      % (4) compute repulsion, impose upper bound for stability, update positions
      power = -1;
      rep(:,:,t) = strength * squeeze(sum((geodesic.^power).*direction,1))';
      pos(:,:,t+1) = pos(:,:,t) + rep(:,:,t);
 
      % check if the paramaters are way off, i.e., repulsion much too weak
      if t == checkpoint 
        earlymaxvar = max(sqrt(sum(var(pos(:,:,1:checkrange),[],3).^2,2)));
        latermaxvar = max(sqrt(sum(var(pos(:,:,checkpoint+1-checkrange:checkpoint),[],3).^2,2)));
        tooslow = latermaxvar > .9*earlymaxvar;
        if tooslow
          display('Too slow. Adapting parameters...')
          strength = strength * 10^flipcounter;
          flipcounter = -sign(flipcounter) * (abs(flipcounter)+1);    % => *1e1,*1e-2,*1e3,*1e-4...
%           display(strength)
          break
        else
          display('Looking good so far, keeping parameters.')
        end
 
      end
 
      % check for convergence, and quit happily if encountered; else continue up to max. NS steps
      if ~mod(t,20) && t >= 40
  %       penulmaxvar = max(sqrt(sum(var(pos(:,:,t-checkpointConv + (1-checkrange:0)),[],3).^2,2)));
        finalmaxvar = max(sqrt(sum(var(pos(:,:,t+(1-checkrange:0)),[],3).^2,2)));
        slowenough = finalmaxvar < 1e-1*latermaxvar;
        closeenough = sqrt(sum(mean(pos(:,:,t),1).^2,2)) < .01;
        convergence = slowenough && closeenough;
        if convergence
          display(['Good enough convergence after ',num2str(t),' steps.'])
          display(['Repulsion strength was ',num2str(strength),'.'])
          break
        end
      end
 
      finalstep = t;
 
    end
 
    if ~convergence && ~tooslow
      display(['No convergence after ',num2str(NS),' steps. Adapting parameters...'])
      strength = strength * 10^flipcounter;
      increment = .5;
      flipcounter = -sign(flipcounter) * (abs(flipcounter)+increment);
%       flipcounter = -sign(flipcounter) * (abs(flipcounter)+1);    % => *1e1,*1e-2,*1e3,*1e-4...
  %     finalstep = NS;
    end
    
    attempt = attempt+1;
    if attempt == 10
      break
    end
         
  end
   
  pos = pos(:,:,1:finalstep);
  finalpos = squeeze(pos(:,:,end));
 
% %     if somecriteriontomeasureconvergence
% %     if (pos(:,:,end-1) - pos(:,:,end)) < 1e-3 * (pos(:,:,20) - pos(:,:,21))
%       convergence = 1;
% %     else
% %       strength = strength * 10^flipcounter;
% %       flipcounter = -sign(flipcounter) * (abs(flipcounter)+1);    % => *1e1, *1e-2, *1e3, *1e-4
% %     end
 
%       normrep = repmat(sum(repulsion.^2,2),[1 ND]);
%       cutrep = repulsion./normrep;
%       replimit = 10;
%       repulsion(normrep>replimit) = cutrep(normrep>replimit);      
 
end
 
 
 
function weightShellDraw(pos,rep)
   
  f = figure(1);
%   set(gcf,'Color',[1 1 1])
%   set(gcf,'Color',[1 1 1],'Position',[360,278,443,420])
  set(gcf,'Color',[1 1 1],'Position',[50,240,440,420])
  clf
   
  NC = size(pos,1);
  ND = size(pos,2);
  NS = size(pos,3);
 
  % plot unit circle
  x = -1:.01:1;
  ypos = sqrt(1-x.^2);
  yneg = -ypos;
 
  plotstep = 1;
  tracecolour = winter(1+ceil(1.2*NS/plotstep));
  colourcounter = 0;
  greycolour = .75*[1 1 1];
     
  if ND == 2
     
    folder = ['weightevolution',date2str('second')];
    mkdir(folder);
       
    for t = 1:plotstep:NS
 
      plot(x,ypos,'Color',.6*[1 1 1])
      xlabel('dimension 1','FontSize',12)
      ylabel('dimension 2','FontSize',12)
      hold on
      plot(x,yneg,'Color',greycolour)
      set(gca,'XLim',[-1.5 1.5],'YLim',[-1.5 1.5],'XTick',[-1 0 1],'YTick',[-1 0 1])
      p = scatter(pos(:,1,t),pos(:,2,t),'o');
      for j = 1:NC
        r(j) = plot([pos(j,1,t),pos(j,1,t)+rep(j,1,t)], ...
                    [pos(j,2,t),pos(j,2,t)+rep(j,2,t)]);
        n(j) = plot([pos(j,1,t)+rep(j,1,t),0], ...
                    [pos(j,2,t)+rep(j,2,t),0],'--','Color',greycolour);
      end
      hold off
 
      colourcounter = colourcounter + 1;
      colour = tracecolour(1+mod(colourcounter,size(tracecolour,1)+1),:);
      set(p,'MarkerEdgeColor',colour)
      set(r,'Color',colour)
      drawnow
      pause(.01)
      saveas(f,[folder,'/weightevolution2d_step',num2str(t),'.fig'])
 
    end
   
  elseif ND == 3
     
    for t = 1:plotstep:NS
 
      % plot unit sphere
      [x,y,z] = sphere(200);
      s = surfl(x,y,z);
      % s = surfl(max(x,zeros(size(x))),max(y,zeros(size(y))),max(z,zeros(size(z))));
%       colormap(gray)
%       set(s,'FaceAlpha',.75)
      colormap([1 1 1])
      set(s,'FaceAlpha',.7)
      shading flat
       
      xlabel('dimension 1','FontSize',12)
      ylabel('dimension 2','FontSize',12)
      zlabel('dimension 3','FontSize',12)
      hold on
      axis equal
      set(gca,'XLim',[-1.2 1.2],'YLim',[-1.2 1.2],'ZLim',[-1.2 1.2], ...
              'XTick',[-1 0 1],'YTick',[-1 0 1],'ZTick',[-1 0 1])
      p = scatter3(pos(:,1,t),pos(:,2,t),pos(:,3,t),'o');
      for j = 1:NC
        r(j) = plot3([pos(j,1,t),pos(j,1,t)+rep(j,1,t)], ...
                     [pos(j,2,t),pos(j,2,t)+rep(j,2,t)], ...
                     [pos(j,3,t),pos(j,3,t)+rep(j,3,t)]);
%         n(j) = plot3([pos(j,1,t)+effrep(j,1,t),0], ...
%                      [pos(j,2,t)+effrep(j,2,t),0], ...
%                      [pos(j,3,t)+effrep(j,3,t),0],'--','Color',greycolour);
      end
      hold off
 
      colourcounter = colourcounter + 1;
      colour = tracecolour(1+mod(colourcounter,size(tracecolour,1)+1),:);
      set(p,'MarkerEdgeColor','none','MarkerFaceColor',colour,'SizeData',70)
      set(r,'Color',colour)
      drawnow
      pause(.01)
 
    end
     
  else
     
    axis off
    tt = text(.5,.5,'Weight evolution will only be displayed in 2D and 3D.');
    set(tt,'HorizontalAlignment','center')
     
  end
 
%   figure(1)
%   % set(gcf,'Color',[1 1 1],'Position',[1361,-175,784,769])
%   clf
% 
%   % plot unit sphere
%   [x,y,z] = sphere(200);
%   s = surfl(x,y,z);
%   % s = surfl(max(x,zeros(size(x))),max(y,zeros(size(y))),max(z,zeros(size(z))));
%   % colormap(gray)
%   % set(s,'FaceAlpha',.3)
%   colormap([1 1 1])
%   set(s,'FaceAlpha',.7)
%   set(gca,'XLim',[-1 1],'YLim',[-1 1],'ZLim',[-1 1])
%   shading flat
% 
%   xlabel('dimension 1','FontSize',12)
%   ylabel('dimension 2','FontSize',12)
%   zlabel('dimension 3','FontSize',12)
% 
%   hold on
% 
%   plotstep = 10;
%   tracecolour = winter(1+NS/plotstep);
%   colourcounter = 0;
% 
%   % for t = [1:5,6:plotstep:NS]
%   for t = 1:plotstep:NS
% 
%     p = scatter3(pos(:,1,t),pos(:,2,t),pos(:,3,t),'o');
%     colourcounter = colourcounter + 1;
%     colour = tracecolour(1+mod(colourcounter,size(tracecolour,1)+1),:);
%     set(p,'MarkerEdgeColor',colour)
%     pause(.001)
% 
%   end
% 
%   hold off
 
end



function [par,his,hon,col,pro] = runSCN(option,feedforward,polytopestretch,randomseed,perturbation)
  
  if isempty(option) || ~isfield(option,'excitation')
    option.excitation = true;  % by default, include lateral excitation
  end
      
  if isempty(option) || ~isfield(option,'lowerbound')
    option.lowerbound = false;  % by default, do not impose a lower bound on voltage
  end
  
  if ~isempty(option) && isfield(option,'refraction')
    if ~isempty(option.refraction)
      par.refraction = option.refraction;
    else
      par.refraction = 0;
    end
  else
    par.refraction = 0;
  end      
      
  if isempty(feedforward)
    par.coverage        = 20;
    par.dimension       = 2;
    par.networksize = par.coverage*par.dimension;
    par.feedforward = randn(par.networksize,par.dimension);
    par.feedforward = par.feedforward ./ sqrt(sum(par.feedforward.^2,2)); % normalisation
  else
    par.feedforward = feedforward;
    par.coverage = size(par.feedforward,1)/size(par.feedforward,2);
    par.dimension = size(par.feedforward,2);
  end
  
  if isempty(polytopestretch)
    par.polytopestretch = 1;
  else
    par.polytopestretch = polytopestretch;
  end
  
  switch option.network
    case {'narrownoise','widenoise'}
      par.delay = 0;
      par.voltage.noisestd = perturbation;
    otherwise
      par.delay = perturbation;
      par.voltage.noisestd = .32;
  end
  
  par.timestep         = 1e-4;
%   if par.delay > 0
%     par.timestep = min([par.timestep,par.delay/2]);
%   end
  par.duration         = 10.5;
  par.decoder.tau      = .01;  % 10 msec, so lambda = 100/sec
  par.decoder.norm     = 1;
  par.l1norm           = 0;
  par.l2norm           = 0*1e-6;
  par.input.magnitude  = 1;
  par.input.interval   = par.duration;
  
  if isempty(randomseed)
    randomseed = randi(1e5);
  end
  
  [par,his]       = setupSCN(option,par,randomseed,option.excitation,option.input);
  [par,his,pro]   = simulateSCN(par,his,option.lowerbound);
  [par,his,hon]   = analyseSCN(par,his);

  if option.showplot
    col           = colourSCN(par);
    [par,his,hon] = plotSCN(par,his,hon,col,option);
  else
    col = [];
  end

end



function [par,his] = setupSCN(option,par,randomseed,excitationpresent,inputtype)
  
  rng(randomseed);

  par.stepdelay = ceil(par.delay/par.timestep);
  par.steprefraction = ceil(par.refraction/par.timestep);
  
% %   par.honeymoon = 2*par.decoder.tau;
  par.honeymoon = .5;
% %   par.honeymoon = .9*par.duration;
%   par.honeymoon = par.timestep;
  
  par.honeymoon = min([par.honeymoon, par.duration*3/4]);
  his.realtime = par.timestep:par.timestep:par.duration;
  par.numtimestep = ceil(par.duration/par.timestep);
  par.input.interval = min(par.input.interval, par.duration);
  par.numinterval = ceil(par.duration/par.input.interval);
  par.numintervaltimestep = ceil(par.input.interval/par.timestep);
%   par.input.value = par.input.magnitude * (1+rand(par.dimension,par.numinterval))/2;
%   par.input.value = par.input.magnitude * (1+randn(par.dimension,par.numinterval))/2;
  par.input.value = par.input.magnitude * (1+ones(par.dimension,par.numinterval))/2;
  par.honeymoonduration = par.honeymoon * par.numinterval;
  
%   inputtype = 'constant';
%   inputtype = 'sine ND';
  switch inputtype
    case 'frozen'
      ladder = .5 : (1/(par.dimension-1)) : 1.5;
      his.input = repmat(ladder',[1 par.numintervaltimestep]);
    case 'random constant'
      par.input.value = par.input.magnitude * (1+randn(par.dimension,par.numinterval))/2;
      his.input = repmat(par.input.value,[1 1 par.numintervaltimestep]);
      his.input = permute(his.input,[1 3 2]);
      his.input = reshape(his.input,[par.dimension,par.numinterval*par.numintervaltimestep]);
      his.input = his.input(:,1:par.numtimestep);
    case 'sine 2D'
      period = 0.1;
      period = 0.5;
      for dimension = 1:par.dimension
        for interval = 1:par.numinterval
          switch dimension
            case 1
              phaseshift = 0;
            case 2
              phaseshift = pi/2;
            otherwise
              phaseshift = pi/2;
          end
          amplitude = 2*par.input.magnitude;
          currentinterval = (interval-1)*par.numintervaltimestep + (1:par.numintervaltimestep);
          switch dimension
            case {1,2}
              his.input(dimension,currentinterval) = ...
                amplitude * sin((1:par.numintervaltimestep)*par.timestep/period * 2*pi + phaseshift);
            otherwise
              his.input(dimension,currentinterval) = ...
                .7*amplitude * 2*(rand(1)-.5) * ones(1,par.numintervaltimestep);
          end
        end
      end

      % normalise, only then apply amplitude
      range = 1:2;
      his.input(range,:) = his.input(range,:) * amplitude / norm(his.input(range,1));
%       if par.dimension > 2
%         for range = 3:par.dimension
%           his.input(range,:) = his.input(range,:) * amplitude / norm(his.input(range,1));
%         end
%       end

    case 'sine ND'
      period = 0.1;
      for dimension = 1:par.dimension
        for interval = 1:par.numinterval
          switch dimension
            case 1
              phaseshift = 0;
            case 2
              phaseshift = pi/2;
            otherwise
              phaseshift = pi/2;
          end
          amplitude = par.input.magnitude;
          currentinterval = (interval-1)*par.numintervaltimestep + (1:par.numintervaltimestep);
%           switch dimension
%             case {1,2}
              his.input(dimension,currentinterval) = ...
                amplitude * sin((1:par.numintervaltimestep)*par.timestep/period * 2*pi + phaseshift);
%             otherwise
%               his.input(dimension,currentinterval) = ...
%                 .7*amplitude * 2*(rand(1)-.5) * ones(1,par.numintervaltimestep);
%           end
        end
      end

      % normalise, only then apply amplitude
      range = 1:par.dimension;
      his.input(range,:) = his.input(range,:) * amplitude / norm(his.input(range,1));
      
    case 'smooth ramp normalised' % as previously used for the paper
      
      if par.numinterval == 1
        
        par.input.signalnoisestd = .5;
        par.input.ramptime = .4;
        par.input.windowtime = 1;

        numrampstep = ceil(par.input.ramptime/par.timestep);

        % ramp up to random target, no signal noise; then stay there
        ramptarget = par.input.magnitude * randn(par.dimension, 1);
        his.input = ramptarget * ones(1,par.numintervaltimestep);
        his.input(:,1:numrampstep) = ramptarget * (0:1/(numrampstep-1):1);
        
        % prepare smoothed signal noise
        slowvariability = par.input.signalnoisestd * randn(par.dimension, par.numintervaltimestep);
        windownumstep = ceil(par.input.windowtime/par.timestep);
        for smooth = 1:2
          slowvariability = cumsum(slowvariability,2);
          slowvariability = movmean(slowvariability,windownumstep,2);
        end
        normeachstep = sqrt(sum(slowvariability.^2,1));
        slowvariability = slowvariability / max(normeachstep);
        
        % add signal noise to signal (including during ramp-up); 1e-4 is NC's usual dt
%         his.input = his.input + slowvariability*sqrt(par.timestep/1e-4);
        his.input = his.input + slowvariability;
%         his.input = his.input * par.timestep;

      end
        
    case 'smooth ramp interpolated' % as used for the paper
      
      if par.numinterval == 1
        
        
        par.input.signalnoisestd = .5;
        par.input.ramptime = .4;
        par.input.windowtime = 1;

        numrampstep = ceil(par.input.ramptime/par.timestep);

        % ramp up to random target, no signal noise; then stay there
        ramptarget = par.input.magnitude * randn(par.dimension, 1);
        his.input = ramptarget * ones(1,par.numintervaltimestep);
        his.input(:,1:numrampstep) = ramptarget * (0:1/(numrampstep-1):1);
        
        % prepare smoothed signal noise
        timestepnuno = 1e-4;
        numtimestepnuno = par.numintervaltimestep * par.timestep/timestepnuno;
        slowvariability = par.input.signalnoisestd * randn(par.dimension, numtimestepnuno);
        windownumstep = ceil(par.input.windowtime/timestepnuno);
        for smooth = 1:2
          slowvariability = cumsum(slowvariability,2);
          slowvariability = movmean(slowvariability,windownumstep,2);
        end
        % do not normalise per time step, but limit max value per dimension
        maxeachdim = max(abs(slowvariability),[],2);
        slowvariability = 0.5 * slowvariability./repmat(maxeachdim,[1 size(slowvariability,2)]);
        % interpolate to account for Nuno's fixed time step
        nunostep = 0:timestepnuno:par.duration-timestepnuno;
        realstep = 0:par.timestep:par.duration-par.timestep;
        for dimindex = 1:size(slowvariability,1)
          scaledvariability(dimindex,:) = interpn(nunostep, slowvariability(dimindex,:), realstep);
        end
        % hack to avoid NaN at the end:
        scaledvariability(:,end) = scaledvariability(:,end-1);
        % add signal noise to signal (including during ramp-up)
        his.input = his.input + scaledvariability;
        
      end
      
    otherwise
      
  end

  
  par.networksize = ceil(par.dimension*par.coverage);
  par.decoder.weight = ones(par.dimension,par.networksize);
  par.voltage.tau = par.decoder.tau;

  par.recurrent = - par.feedforward * par.feedforward' - par.l2norm * diag(par.networksize);
  
  if ~excitationpresent
    %%% MODIFICATION: REMOVE RECURRENT EXCITATION
    par.recurrent(par.recurrent>0) = 0;
  end
  
%   removestrongexcitation = true;
  if excitationpresent && option.removestrongexcitation
%     removepercentage = .2; % percentage of excitatory weights to be removed
    removepercentage = option.fractionremoved; % percentage of excitatory weights to be removed
    [~,weightindex] = sortrows(par.recurrent(:),'descend');
    % both negative and positive weights are included here, so percentage
    % defined above is twice what it should be for this vector:
    par.recurrent(weightindex(1:round(removepercentage/2*numel(weightindex)))) = 0;
  end

  par.threshold = 1/2 * sum(par.feedforward.^2,2) + par.l1norm + par.l2norm;
  par.threshold = par.threshold * par.polytopestretch; % choose a wider box if desired
  
  his.voltage(:,1) = zeros(par.networksize,1);
  his.spike = zeros(par.networksize,par.numtimestep);
%   his.arrival = zeros(par.networksize,par.numtimestep,par.networksize);
  his.delayedspike = zeros(par.networksize,par.numtimestep);
  his.refraction = zeros(par.networksize,par.numtimestep);

end



function [par,his,pro] = simulateSCN(par,his,lowerboundpresent)

%   rng(123); % to ensure same frozen noise each time

  if lowerboundpresent
    lowerbound = -par.threshold;
  end

  % FD20200224: delay probe
  pro.continuousSpiking = false([par.numtimestep 1]);
  entity = 0;
  
  for t = 1:par.numtimestep-1

    displaySimulationProgressSCN(t,par.numtimestep,.2);
    
    % First, check for above-threshold neurons...
    abovethreshold = find(his.voltage(:,t) >= par.threshold);
    
    while ~isempty(abovethreshold)
      
      % ...and keep checking until none are left.
      abovethreshold = find(his.voltage(:,t) >= par.threshold);
      
      % Second, ignore above-threshold neurons that are currently in refraction.
      if sum((his.refraction(:,t))) > 0
        inactive = find(his.refraction(:,t));
        abovethreshold = setdiff(abovethreshold,inactive);
      end
      
      if ~isempty(abovethreshold) % Skip this during final iteration of while loop.
        
        % Third, allow the one neuron furthest above threshold to spike.
        [~,furthestabove] = max(his.voltage(abovethreshold,t));
        spiker = abovethreshold(furthestabove);
        
        % FD20200224: delay probe
%         if t > par.honeymoon/par.timestep % don't probe during honeymoon
          pro.continuousSpiking(t) = true;
          if t == 1 || (pro.continuousSpiking(t-1) == 0 && (sum(his.spike(:,t)) == 0))
            entity = entity + 1;
            pro.numCrossTotal(entity) = 0;
            pro.numCrossFirst(entity) = 0;
            pro.crossListTotal{entity} = [];
            pro.crossListFirst{entity} = [];
            pro.seizureStart(entity) = t; % only valid if "pro" aligned with "his"
%             pro.seizureStart(entity) = t - par.honeymoon/par.timestep; % "pro" aligned with "hon"
          end
          pro.numCrossTotal(entity) = pro.numCrossTotal(entity) + 1;
          pro.crossListTotal{entity} = [pro.crossListTotal{entity}, spiker];

          continuousNumStep = t - find(diff(pro.continuousSpiking)==1,1,'last');
          if continuousNumStep*par.timestep <= par.delay
  %           voltagenow = sortrows(his.voltage(:,t))
            pro.numCrossFirst(entity) = pro.numCrossFirst(entity) + 1;
            pro.crossListFirst{entity} = [pro.crossListFirst{entity}, spiker];
          else  % sum up the results and compare to undelayed case
  %           if pro.numCrossTotal(entity) > 1
  %             voltageBeforeJump = his.voltage(:,t-1);
  %             readout = voltageBeforeJump' / par.feedforward'
  %             norm(readout)
  %             pro.crossListFirst{entity}
  %           end
          end
%         end
        
        his.spike(spiker,t) = 1;  % record spike for posterity
        his.delayedspike(spiker,t+par.stepdelay) = 1;  % create subsequent delayed recurrence
        his.refraction(spiker,t+(0:par.steprefraction)) = 1;  % create subsequent refraction
        % Beginning refraction at the current time step t+0 prevents multiple
        % spikes from the same neuron during a single time step and, by extension,
        % infinite ping-pong loops during simulation.

        if par.delay == 0
          
          % Fourth (in the absence of delays), full instantaneous lateral recurrence.
          % If applicable, enforce lower limit on voltage.
          instantEffect = par.recurrent(:,spiker);
          his.voltage(:,t) = his.voltage(:,t) + instantEffect;
          if lowerboundpresent
            his.voltage(:,t) = max(his.voltage(:,t), lowerbound);
          end
          
        else
          
          % Fourth (with finite delays), effect instantaneous self-reset.
          % If applicable, enforce lower limit on voltage.
          instantEffect = par.recurrent(spiker,spiker);
          his.voltage(spiker,t) = his.voltage(spiker,t) + instantEffect;
          if lowerboundpresent
            his.voltage(:,t) = max(his.voltage(:,t), lowerbound);
          end

        end
                
      end
            
      % Repeat checking for spikes until there are no more.
      
%       % DEBUGGING ONLY: Allow at most one spike per time step.
%       break

    end
    
    
    % FD20200224: delay probe
    % if there was no spike this step, and there is no remaining delayed
    % recurrence, end recording the current phase of "continuous spiking"
    % Otherwise, continue tracking the current phase or "epoch":
    if pro.continuousSpiking(t) && sum(sum(his.delayedspike(:,t:end))) > 0
      pro.continuousSpiking(t+1) = true; % Yes, "true"!
    end


    % Fifth (and only with finite delays), delayed lateral recurrence.
    % For simultaneous spikes, all such recurrence is applied at once, not in any order.
    % If applicable, enforce lower limit on voltage.
    if par.delay > 0
      originlist = find(his.delayedspike(:,t));
      if ~isempty(originlist)
        for sender = originlist'
          recipient = setdiff(1:par.networksize,sender);
          delayedEffect = par.recurrent(recipient,sender);
          his.voltage(recipient,t) = his.voltage(recipient,t) + delayedEffect;
        end
      end
      % Apply lower bound only after excitation and inhibition had a chance to balance out.
      if lowerboundpresent
        his.voltage(:,t) = max(his.voltage(:,t), lowerbound);
      end
    end
    
    % Finally, apply voltage leak, input currents and voltage. Move on to next time step.
    % If applicable, enforce lower limit on voltage.
    dvdt = 1/par.decoder.tau * (- his.voltage(:,t) + par.feedforward * his.input(:,t));
%     dxdt = diff(his.input(:,[t t+1]),[],2) / par.timestep;
%     dvdt = - 1/par.decoder.tau * his.voltage(:,t) ...
%            + par.feedforward * (1/par.decoder.tau*his.input(:,t) + dxdt);     
    voltagenoise = par.voltage.noisestd * randn(par.networksize,1);
    his.voltage(:,t+1) = his.voltage(:,t) + dvdt*par.timestep + voltagenoise*sqrt(par.timestep);
    if lowerboundpresent
      his.voltage(:,t+1) = max(his.voltage(:,t+1), lowerbound);
    end   
    
  end
  
end



function displaySimulationProgressSCN(step,numstep,progressinterval)

  if mod(step,numstep*progressinterval) == 0
    disp(['Progress: ',num2str(round(100*step/numstep)),'%'])
  elseif step == numstep-1
    disp('Progress: 100%')
  end
    
end



function [par,his,hon] = analyseSCN(par,his)
  
  parallelcomputing = false;

  his.rate = sum(his.spike,2)/par.duration;
%   kernel = exp(-(par.timestep:par.timestep:10*par.decoder.tau)/par.decoder.tau);
%   kernel = exp(-his.realtime(1:ceil(end/10))/par.decoder.tau);
  kernel = exp(-his.realtime/par.decoder.tau);
%   % kernel = par.timestep * exp(-his.realtime/par.decoder.tau)/par.decoder.tau;

  if parallelcomputing

    parin = his.spike; % "parfor" recommends this assignment outside to reduce communication overhead
    if gpuDeviceCount > 0
      parin = gpuArray(parin);
      kernel = gpuArray(kernel);
      warning('off','parallel:gpu:device:DeviceLibsNeedsRecompiling')
    end
    parfor k = 1:par.networksize
        parout(k,:) = conv(parin(k,:),kernel,'full');
    end
  %   parout = gather(parout);
    his.psp = parout; % "parfor" requires this assignment outside(!) to classify variable

    parin2 = his.input;
    if gpuDeviceCount > 0
      parin2 = gpuArray(parin2);
      warning('off','parallel:gpu:device:DeviceLibsNeedsRecompiling')
    end
    parfor k = 1:par.dimension
      parout2(k,:) = conv(parin2(k,:),kernel,'full') / sum(kernel);
  %     parout2(k,:) = conv(parin2(k,:),kernel,'full');
    end
    his.target = parout2;
    
  else
    
    for k = 1:par.networksize
      his.psp(k,:) = conv(his.spike(k,:),kernel,'full');
  %     his.psp(k,:) = fastconvolution(his.spike(k,:),kernel);
    end
    
    % for k = 1:par.dimension
    %   his.target(k,:) = conv(his.input(k,:),kernel,'same') * par.decoder.tau;
    % end
    for k = 1:par.dimension
      his.target(k,:) = conv(his.input(k,:),kernel,'full') / sum(kernel);
    end    
    
  end
  
  his.target = his.target(:,1:par.numtimestep);
  his.psp = his.psp(:,1:par.numtimestep);
  his.psp = gather(his.psp); % added to avoid bug on gpu21
  his.readout = par.feedforward'*his.psp;
  
  % readout scaling:
  timeaveragednorm = mean(sqrt(sum(his.readout.^2,1)));
  scalingfactor = (timeaveragednorm + (par.polytopestretch-1)/2) / timeaveragednorm;
  his.readout = his.readout * scalingfactor;
  % bugfix: In MATLAB, 0 (readout) * Inf (scalingfactor) = NaN. We need it to return 0 instead:
  his.readout(isnan(his.readout)) = 0;
  
  his.target = gather(his.target); % added to avoid bug on gpu21
  
  his.error = his.target-his.readout;
  
  [par,his,hon] = honeymoonSCN(par,his);

end



function [par,his,hon] = honeymoonSCN(par,his)

  % identify honeymoon period
  his.honeymoon.x = [];
  his.honeymoon.y = [];
  his.honeymoon.timeindex = [];
  for k = 1:par.numinterval
    newx = (k-1)*par.input.interval + [0 0 par.honeymoon*[1 1]];
    his.honeymoon.x = [his.honeymoon.x NaN newx];
    his.honeymoon.y = [his.honeymoon.y NaN -10 1 1 -10];
    newindex = (k-1)*par.numintervaltimestep + (1:ceil(par.honeymoon/par.timestep));
    his.honeymoon.timeindex = [his.honeymoon.timeindex; newindex'];
  end

  % disregard honeymoon period for subsequent analysis - not before
  hon.spike = his.spike;
  hon.spike(:,his.honeymoon.timeindex) = 0;
  hon.rate = sum(hon.spike,2)/par.honeymoonduration;
  hon.error = his.error;
  hon.error(:,his.honeymoon.timeindex) = NaN;
  
  % not necessary, but convenient:
  hon.spike    = his.spike(:,his.honeymoon.timeindex(end):end);
  hon.readout  = his.readout(:,his.honeymoon.timeindex(end):end);
  hon.target   = his.target(:,his.honeymoon.timeindex(end):end);
  hon.error    = his.error(:,his.honeymoon.timeindex(end):end);
  hon.realtime = his.realtime(his.honeymoon.timeindex(end):end);

  for int = 1:par.numinterval

    segment(int).step = (int-1)*par.numintervaltimestep + ...
                        ((ceil(par.honeymoon/par.timestep)+1) : par.numintervaltimestep);
    segment(int).spike = his.spike(:,segment(int).step);

    % pool across neurons within interval
    segment(int).cvpool = [];
    for neuron = 1:par.networksize
      singlespike = segment(int).spike(neuron,:);
      singlesisi = diff(find(singlespike)); % one neuron, one interval
      if numel(singlesisi(~isnan(singlesisi))) > 0
        segment(int).std(neuron) = std(singlesisi);
        segment(int).mean(neuron) = mean(singlesisi);
        segment(int).median(neuron) = median(singlesisi);
        newcv = segment(int).std(neuron) / segment(int).mean(neuron);
        segment(int).cvpool = [segment(int).cvpool newcv];
      end
    end
  
  end
  % pool across intervals
  hon.cv = [segment.cvpool];

end



function col = colourSCN(par)

  col.honeymoon = .9*[1 1 1];
%   col.voltage = flipud(cbrewer('seq','Greens',ceil(1.5*par.networksize)));
%   col.readout = flipud(cbrewer('seq','Purples',ceil(1.5*par.dimension)));
  col.voltage = flipud(cbrewer('seq','Purples',ceil(1.5*par.networksize)));
  if par.delay > 0
    col.readout = flipud(cbrewer('seq','Blues',ceil(1.5*par.dimension)));
  else
    col.readout = flipud(cbrewer('seq','Greens',ceil(1.5*par.dimension)));
  end
  col.target = col.readout * .9;
%   col.target = ones(size(col.readout));
  col.spike = .1*[1 1 1];
  col.input = col.readout;
  col.error = col.readout;
  col.weight = col.voltage(ceil(par.networksize/2),:);
  col.histogram = col.voltage(ceil(par.networksize*1/4),:);
  
end



function saveSCN(par,his,hon,honReference,performance,networkerror,option)

%   hon_partial.cv = hon.cv;
%   hon_partial.numspike = sum(hon.spike,2);
%   hon_partial.excess = sum(hon.spike,2)-sum(honReference.spike,2);
%   hon_partial.rate = hon.rate;
% %   hon_partial.readout = hon.readout;
% %   hon_partial.target = hon.target;
  hon_partial = hon;
  clearvars his hon honReference
  
  hon = hon_partial;
  for k = 1:numel(par), par(k).recurrent = []; end
  elapsedtime = toc;

%   savefolder = 'C:\Users\fdehmelt\Desktop\SingleNetworkFigure\';
  runtime = clock;
  timestring = ['_', ...
                num2str(runtime(1),'%.4u'), ...
                num2str(runtime(2),'%.2u'), ...
                num2str(runtime(3),'%.2u'), ...
                '_', ...
                num2str(runtime(4),'%.2u'), ...
                num2str(runtime(5),'%.2u'), ...
                num2str(runtime(6),'%02.0f')];
  switch option.network
    case 'full'
      leadstring = 'scnDel';
    case 'wide'
      leadstring = 'wide';
    case 'reduced'
      leadstring = 'redExc';
    case {'narrowNoise','wideNoise'}
      leadstring = 'scnNoi';
    case 'none'
      leadstring = 'iolb';
    otherwise
      leadstring = 'scnDel';
  end
  savefile = [leadstring, ...
    '_',num2str(par.dimension),'d', ...
    '_',num2str(par.coverage),'c', ...
    '_',num2str(par.delay),'th', ...
    '_',num2str(par.polytopestretch),'b', ...
    '_',num2str(par.duration),'s',timestring,'.mat'];
  while exist(savefile,'file') == 2
    savefile = [savefile(1:end-4),'a','.mat'];
  end
  
  save(savefile)
  
end



function z = fastconvolution(x,y)

  lengthafterconvolution = length(x) + length(y) - 1;
  smallestpower = pow2(nextpow2(lengthafterconvolution));    % Find smallest power of 2 that is > Ly
  X = fft(x,smallestpower);
  Y = fft(y,smallestpower);
  Z = X.*Y;        	           % 
  z = real(ifft(Z,smallestpower));
  z = z(1:lengthafterconvolution);

end



function [performance,networkerror] = ...
  NFPCperformance(target,readout,referencereadout)
  
  errorDead      = NFPCerror(target, zeros(size(target)));
  errorTrial     = NFPCerror(target, readout);
  errorReference = NFPCerror(target, referencereadout);
  
  performance = (errorDead - errorTrial) / (errorDead - errorReference);
  
  networkerror.errorDead      = errorDead;
  networkerror.errorTrial     = errorTrial;
  networkerror.errorReference = errorReference;
  
end



function meanerrornorm = NFPCerror(target,readout)

  codingerror = readout - target;
  l2normoferror = sqrt(sum(codingerror.^2,1));
  meanerrornorm = mean(l2normoferror);
  
end



function b = cutdowntosize(a,b)

  if size(b,2) > numel(a)
    b = b(:,1:numel(a));
  end

end



function boxchoicelist = readOptimalBoxSize
  
  % dimension, coverage, delay, box size
%   % updated on 20200414
%   boxchoicelist = [ 2,   2, 1e-3, 0.864; ...
%                     2,   5, 1e-3, 1.18; ...
%                     2,  10, 1e-3, 1.68; ...
%                     2,  20, 1e-3, 2.51; ...
%                     2,  50, 1e-3, 2.65; ...
%                     5,   2, 1e-3, 1.06; ...
%                     5,   5, 1e-3, 1.18; ...
%                     5,  10, 1e-3, 1.44; ...
%                     5,  20, 1e-3, 1.86; ...
%                     5,  50, 1e-3, 2.27; ...
%                    10,   2, 1e-3, 1.24; ...
%                    10,   5, 1e-3, 1.24; ...
%                    10,  10, 1e-3, 1.44; ...
%                    10,  20, 1e-3, 1.51; ...
%                    10,  50, 1e-3, 2.05; ...
%                    20,   2, 1e-3, 0.864; ...
%                    20,   5, 1e-3, 1.18; ...
%                    20,  10, 1e-3, 1.44; ...
%                    20,  20, 1e-3, 1.68; ...
%                    20,  50, 1e-3, 1.86; ...
%                    50,   2, 1e-3, 0.825; ...
%                    50,   5, 1e-3, 1.12; ...
%                    50,  10, 1e-3, 1.24; ...
%                    50,  20, 1e-3, 1.44; ...
%                    50,  50, 1e-3, 1.68];

  % updated on 20200414
  boxchoicelist = [ 2,   2, 1e-3, 1.06; ...
                    2,   5, 1e-3, 1.24; ...
                    2,  10, 1e-3, 1.86; ...
                    2,  20, 1e-3, 2.65; ...
                    2,  50, 1e-3, 3.07; ...
                    5,   2, 1e-3, 1.24; ...
                    5,   5, 1e-3, 1.24; ...
                    5,  10, 1e-3, 1.51; ...
                    5,  20, 1e-3, 1.94; ...
                    5,  50, 1e-3, 2.27; ...
                   10,   2, 1e-3, 1.24; ...
                   10,   5, 1e-3, 1.24; ...
                   10,  10, 1e-3, 1.44; ...
                   10,  20, 1e-3, 1.68; ...
                   10,  50, 1e-3, 2.15; ...
                   20,   2, 1e-3, 1.24; ...
                   20,   5, 1e-3, 1.24; ...
                   20,  10, 1e-3, 1.44; ...
                   20,  20, 1e-3, 1.68; ...
                   20,  50, 1e-3, 1.86; ...
                   50,   2, 1e-3, 1.01; ...
                   50,   5, 1e-3, 1.18; ...
                   50,  10, 1e-3, 1.3; ...
                   50,  20, 1e-3, 1.44; ...
                   50,  50, 1e-3, 1.68];

  % updated on 20200419 (wrong 1msec delay for 5M20R fixed on 20200424)
  boxchoicelist = [boxchoicelist; ...
                    2,   2, 2e-3, 1.36; ...
                    2,   5, 2e-3, 1.86; ...
                    2,  10, 2e-3, 2.15; ...
                    2,  20, 2e-3, 2.78; ...
                    2,  50, 2e-3, 3.77; ...
                    5,   2, 2e-3, 1.76; ...
                    5,   5, 2e-3, 1.76; ...
                    5,  10, 2e-3, 2.27; ...
                    5,  20, 2e-3, 2.4; ...
                    5,  50, 2e-3, 2.78; ...
                   10,   2, 2e-3, 2.05; ...
                   10,   5, 2e-3, 2.1; ...
                   10,  10, 2e-3, 2.39; ...
                   10,  20, 2e-3, 2.78; ...
                   10,  50, 2e-3, 2.78; ...
                   20,   2, 2e-3, 2.05; ...
                   20,   5, 2e-3, 2.51; ...
                   20,  10, 2e-3, 2.78; ...
                   20,  20, 2e-3, 2.78; ...
                   20,  50, 2e-3, 2.78; ...
                   50,   2, 2e-3, 2.78; ...
                   50,   5, 2e-3, 2.78; ...
                   50,  10, 2e-3, 2.78; ...
                   50,  20, 2e-3, 2.78; ...
                   50,  50, 2e-3, 2.78];

	% Remember: the reference network now always has box size 0.55!
  
end



function excitationchoicelist = readOptimalExcitationRemoved

  % dimension, coverage, delay, fraction of excitation to be removed
  % updated on 20200414
  excitationchoicelist = [ 2,   2, 1e-3, 0.39; ...
                           2,   5, 1e-3, 0.68; ...
                           2,  10, 1e-3, 0.68; ...
                           2,  20, 1e-3, 0.68; ...
                           2,  50, 1e-3, 0.70; ...
                           5,   2, 1e-3, 0.4; ...
                           5,   5, 1e-3, 0.4; ...
                           5,  10, 1e-3, 0.4; ...
                           5,  20, 1e-3, 0.5; ...
                           5,  50, 1e-3, 0.5; ...
                          10,   2, 1e-3, 0.25; ...
                          10,   5, 1e-3, 0.38; ...
                          10,  10, 1e-3, 0.38; ...
                          10,  20, 1e-3, 0.38; ...
                          10,  50, 1e-3, 0.38; ...
                          20,   2, 1e-3, 0.2; ...
                          20,   5, 1e-3, 0.27; ...
                          20,  10, 1e-3, 0.27; ...
                          20,  20, 1e-3, 0.27; ...
                          20,  50, 1e-3, 0.27; ...
                          50,   2, 1e-3, 0.1; ...
                          50,   5, 1e-3, 0.1; ...
                          50,  10, 1e-3, 0.18; ...
                          50,  20, 1e-3, 0.18; ...
                          50,  50, 1e-3, 0.18];
  
  % updated on 20200419
  excitationchoicelist = [excitationchoicelist; ...
                           2,   2, 2e-3, 1; ...
                           2,   5, 2e-3, .9; ...
                           2,  10, 2e-3, .71; ...
                           2,  20, 2e-3, .7; ...
                           2,  50, 2e-3, .7; ...
                           5,   2, 2e-3, 1; ...
                           5,   5, 2e-3, 1; ...
                           5,  10, 2e-3, .59; ...
                           5,  20, 2e-3, .52; ...
                           5,  50, 2e-3, .52; ...
                          10,   2, 2e-3, 1; ...
                          10,   5, 2e-3, 0.82; ...
                          10,  10, 2e-3, 0.61; ...
                          10,  20, 2e-3, 0.4; ...
                          10,  50, 2e-3, 0.4; ...
                          20,   2, 2e-3, 1; ...
                          20,   5, 2e-3, .58; ...
                          20,  10, 2e-3, .39; ...
                          20,  20, 2e-3, .26; ...
                          20,  50, 2e-3, .1; ...
                          50,   2, 2e-3, .6; ...
                          50,   5, 2e-3, .31; ...
                          50,  10, 2e-3, .29; ...
                          50,  20, 2e-3, .17; ...
                          50,  50, 2e-3, .1];
  
end



function [par,his,hon,col] = plotSCN(par,his,hon,col,option)
  
  % quick fix for ceil/floor rounding errors when setting array size
  his.voltage = cutdowntosize(his.realtime,his.voltage);
  his.readout = cutdowntosize(his.realtime,his.readout);
  his.target  = cutdowntosize(his.realtime,his.target);
  his.error   = cutdowntosize(his.realtime,his.error);
  his.input   = cutdowntosize(his.realtime,his.input);
  
  % sort neurons by first-dimension weight (for prettiness only)
%   [~,sortindex] = sortrows(par.feedforward(:,1));
%   [~,sortindex] = sortrows(atan(par.feedforward(:,2)./par.feedforward(:,1)));
  angle = -1i*log(par.feedforward(:,1) + 1i*par.feedforward(:,2));
  [~,sortindex] = sortrows(real(angle));
  his.spike = his.spike(sortindex,:);  
  
  
  fig = figure();
  set(gcf,'Color',[1 1 1],'Position',[5 45 1500 950])

  a3 = axes('Position',[.25 .74 .25 .22]);
  y3 = [min(min(his.voltage)) - .1*abs(min(min(his.voltage))), ...
        max(max(his.voltage)) + .1*abs(max(max(his.voltage)))];
  hold on
  p32 = plot(his.realtime,his.voltage');
  p31 = area(his.honeymoon.x,y3(2)*his.honeymoon.y,y3(1));
  % p31 = area([0 par.honeymoon], y3(2)*[1 1], y3(1));
  % set(p31,'FaceColor',col.honeymoon,'EdgeColor','none')
  set(p31,'FaceColor',[0 0 0],'FaceAlpha',.1,'EdgeColor','none')
  hold off
  for line = 1:numel(p32), set(p32(line),'Color',col.voltage(line,:)), end
  title('voltage')
  xlabel('time (sec)')
  set(a3,'Layer','top','Box','on')
  set(a3,'YLim',y3);

  a4 = axes('Position',[.55 .74 .25 .22]);
  train.markerlength = .8;
  y4 = [1-train.markerlength*3/4, par.networksize+train.markerlength*3/4];
  [train.cell,train.timeindex] = find(his.spike);
  train.time = his.realtime(train.timeindex)';
  train.x = [train.time, train.time, NaN(size(train.time))];
  train.x = reshape(train.x',[3*numel(train.timeindex),1]);
  train.y = [train.cell - train.markerlength/2, ...
             train.cell + train.markerlength/2, NaN(size(train.cell))];
  train.y = reshape(train.y',[3*numel(train.cell),1]);
  hold on
  % p41 = area([0 par.honeymoon], y4(2)*[1 1], y4(1));
  p41 = area(his.honeymoon.x,y4(2)*his.honeymoon.y,y4(1));
  p42 = plot(train.x,train.y);
  hold off
  set(p42,'Color',col.spike)
  title('spikes')
  xlabel('time (sec)')
  colormap('gray')
  set(p41,'FaceColor',col.honeymoon,'EdgeColor','none')
  % set(a4,'XLim',[his.realtime(1), his.realtime(end)])
  set(a4,'XLim',[0, his.realtime(end)])
  set(a4,'YLim',y4)
  set(a4,'Layer','top','Box','on')

  a5 = axes('Position',[.85 .74 .10 .22]);
  b5 = barh(his.rate);
  title('firing rate')
  xlabel('rate over trial (Hz)')
  set(a5,'YLim',y4)
  set(b5,'FaceColor',.6*[1 1 1])
  set(b5,'EdgeColor',get(b5,'FaceColor'))
  set(a5,'TickLength',[0 0])
  set(a5,'XGrid','on')

  a6 = axes('Position',[.85 .42 .10 .22]);
  b6 = barh(sortrows(his.rate));
  title('firing rate (sorted)')
  xlabel('rate over trial (Hz)')
  set(a6,'YLim',y4)
  set(b6,'FaceColor',.6*[1 1 1])
  set(b6,'EdgeColor',get(b6,'FaceColor'))
  set(a6,'TickLength',[0 0])
  set(a6,'XGrid','on')

  a7 = axes('Position',[.55 .42 .25 .22]);
  minimum = min([min(min(his.readout)) min(min(his.target))]);
  maximum = max([max(max(his.readout)) max(max(his.target))]);
  y7 = [minimum - .1*minimum, maximum + .1*maximum];
  hold on
  % p71 = area([0 par.honeymoon], y7(2)*[1 1], y7(1));
  p72 = plot(his.realtime,his.readout');
  p73 = plot(his.realtime,his.target');
  p71 = area(his.honeymoon.x,y7(2)*his.honeymoon.y,y7(1));
  % set(p71,'FaceColor',col.honeymoon,'EdgeColor','none')
  set(p71,'FaceColor',[0 0 0],'FaceAlpha',.1,'EdgeColor','none')
  for line = 1:numel(p72), set(p72(line),'Color',col.readout(line,:)), end
  for line = 1:numel(p73), set(p73(line),'Color',col.target(line,:)), end
  set(p73,'LineWidth',1)
  p7372 = [p73,p72]';
  uistack(p7372(:),'bottom')
  hold off
  title('readout')
  xlabel('time (sec)')
  ylabel(option.network)
  set(get(a7,'YLabel'),'Visible','off')
  set(a7,'Layer','top','Box','on')
  set(a7,'YLim',y7)

  a8 = axes('Position',[.25 .42 .25 .22]);
  y8 = [min(min(his.error)) - .1*abs(min(min(his.error))), ...
        max(max(his.error)) + .1*abs(max(max(his.error)))];
  hold on
  p82 = plot(his.realtime,his.error');
  p81 = area(his.honeymoon.x,y8(2)*his.honeymoon.y,y8(1));
  % p81 = area([0 par.honeymoon], y8(2)*[1 1], y8(1));
  % set(p81,'FaceColor',col.honeymoon,'EdgeColor','none')
  set(p81,'FaceColor',[0 0 0],'FaceAlpha',.1,'EdgeColor','none')
  hold off
  for line = 1:numel(p82), set(p82(line),'Color',col.error(line,:)), end
  title('decoding error')
  xlabel('time (sec)')
  set(a8,'Layer','top','Box','on')
  set(a8,'YLim',y8);

  a9 = axes('Position',[.25 .10 .186 .22]);
  title('max. error histogram')
  xlabel({'maximum absolute difference over dimensions','(pooled across individual steps)'})
  % xlabel('$\sqrt{\frac{\sum_{k=1}^N(x_k-\hat{x}_k)^2}{N}}$  [incorrect label]','Interpreter','latex','FontSize',14)
  numbin = 180;
  data = max(abs(hon.error.^2),[],1); % "hon" (not "his") excludes honeymoon periods
  xmin = min(data(data>0));
  xmax = max(data(data>0));
  if isempty(xmin),  xmin = 1e-9;  end
  if isempty(xmax),  xmax = 1;     end
%     if xmin == xmax,  xmax = xmin + 1e-9*abs(xmin);  end
  edge = logspace(log10(xmin),log10(xmax),numbin);
  % edge = logspace(-.1,.1,numbin);
  [count,edge] = histcounts(data,edge);
  normalisedcount = count/sum(count);
  hold on
  p91 = histogram(data,edge,'Normalization','probability');
  % bugfix: for values == 0, "stairs" does not draw vertical lines
  bugfix = 1e-30;
%     p92 = stairs(edge(1:end-1),normalisedcount+bugfix);
  p92 = stairs(edge(1:end-1),p91.Values+bugfix);
  hold off
  box off
  set(a9,'XScale','log','YScale','log')
  set(a9,'YLim',[1e-5 1e-1])
  set(p91,'FaceColor',col.histogram,'EdgeColor','none')
  set(p92,'Color','k')

  a10 = axes('Position',[.517 .10 .186 .22]);
  title('rate histogram')
  xlabel({'firing rate (Hz)','(pooled across neurons and intervals)'})
  numbin = 30;
  data = hon.rate; % "hon" (not "his") excludes honeymoon periods
  xmin = min(data(data>0));
  xmax = max(data(data>0));
  if isempty(xmin),  xmin = 1e-9;  end
  if isempty(xmax),  xmax = 1;     end
%     if xmin == xmax,  xmax = xmin + 1e-9*abs(xmin);  end
  edge = logspace(log10(xmin),log10(xmax),numbin);
  [count,edge] = histcounts(data,edge);
  normalisedcount = count/sum(count);
  hold on
  p101 = histogram(data,edge,'Normalization','probability');
  % bugfix: for values == 0, "stairs" does not draw vertical lines
  bugfix = 1e-30;
%     p102 = stairs(edge(1:end-1),normalisedcount+bugfix);
  p102 = stairs(edge(1:end-1),p101.Values+bugfix);
  hold off
  box off
  set(a10,'XScale','log','YScale','log')
  set(a10,'XLim',[1/par.input.interval 1/par.timestep],'YLim',[1e-2 3e-1])
  set(p101,'FaceColor',col.histogram,'EdgeColor','none')
  set(p102,'Color','k')


  a11 = axes('Position',[.784 .10 .186 .22]);
  title('CV histogram')
  xlabel({'coefficient of variation','(pooled across neurons and intervals)'})
  numbin = 30;
  data = hon.cv; % "hon" (not "his") excludes honeymoon periods
  xmin = min(data(data>0));
  xmax = max(data(data>0));
  if isempty(xmin),  xmin = 1e-9;  end
  if isempty(xmax),  xmax = 1;     end
%     if xmin == xmax,  xmax = xmin + 1e-9*abs(xmin);  end
  edge = logspace(log10(xmin),log10(xmax),numbin);
  [count,edge] = histcounts(data,edge);
  normalisedcount = count/sum(count);
  hold on
  p111 = histogram(data,edge,'Normalization','probability');
  % bugfix: for values == 0, "stairs" does not draw vertical lines
  bugfix = 1e-30;
%     p112 = stairs(edge(1:end-1),normalisedcount+bugfix);
  p112 = stairs(edge(1:end-1),p111.Values+bugfix);
  hold off
  box off
  set(a11,'XScale','log','YScale','lin')
  set(p111,'FaceColor',col.histogram,'EdgeColor','none')
  set(p112,'Color','k')


  a12 = axes('Position',[.05 .10 .15 .40]);
  title('simulation parameters')
  t12(1) = text(0,-1,['dimension = ',num2str(par.dimension)]);
  t12(2) = text(0,-2,['coverage = ',num2str(par.coverage)]);
  t12(3) = text(0,-3,['decoder \tau (sec) = ',num2str(par.decoder.tau)]);
  t12(4) = text(0,-4,['time step (sec) = ',num2str(par.timestep)]);
  t12(5) = text(0,-5,['delay (sec) = ',num2str(par.delay)]);
  t12(6) = text(0,-6,['L1 = ',num2str(par.l1norm)]);
  t12(7) = text(0,-7,['L2 = ',num2str(par.l2norm)]);
  t12(8) = text(0,-8,['noise st.dev. = ',num2str(par.voltage.noisestd)]);
  t12(9) = text(0,-9,['duration (sec) = ',num2str(par.duration)]);
  t12(10) = text(0,-10,['input segment (sec) = ',num2str(par.input.interval)]);
  t12(11) = text(0,-11,['honeymoon period (sec) = ',num2str(par.honeymoon)]);
  t12(12) = text(0,-12,['polytope stretch = ',num2str(par.polytopestretch)]);
  t12(13) = text(0,-13,['input strength = ',num2str(par.input.magnitude)]);
  % t12(14) = text(0,-14,[' = ',num2str(par.)]);
  % t12(15) = text(0,-15,[' = ',num2str(par.)]);
  set(a12,'XLim',[-.1 1],'YLim',[-14.2 0])
  set(a12,'XTick',[],'YTick',[])
  set(a12,'Box','on')


  a2 = axes('Position',[.05 .605 .15 .12]);
  p21 = plot(his.realtime,his.input');
  for line = 1:numel(p21), set(p21(line),'Color',col.input(line,:)), end
  set(p21,'LineWidth',2)
  title('input current')
  xlabel('time (sec)')
  set(a2,'XLim',[0 par.duration])

  if par.dimension <= 3
    a1 = axes('Position',[.05 .8 .15 .15]);
    if par.dimension == 15
      p11 = scatter(par.feedforward(:,1));
    elseif par.dimension == 2
      hold on
      circle.x = -1:.01:1;
      circle.y1 = +sqrt(1-circle.x.^2);
      circle.y2 = -sqrt(1-circle.x.^2);
      plot(circle.x,circle.y1,'-k')
      plot(circle.x,circle.y2,'-k')
      p11 = scatter(par.feedforward(:,1),par.feedforward(:,2));
      hold off
    elseif par.dimension == 3
      p11 = scatter3(par.feedforward(:,1),par.feedforward(:,2),par.feedforward(:,3));
    end
    set(p11,'MarkerEdgeColor',col.weight)
    title('feedforward weights')
    set(gca,'XLim',[-1.1 1.1],'YLim',[-1.1 1.1],'ZLim',[-1.1 1.1])
    axis square
  end
   
  % hgexport(fig,[savefolder,savefile])
  % print([savefolder,savefile],'-dtiffn')
  % saveas(fig,[savefolder,savefile],'tif')
  % save([savefolder,savefile,'.mat'])
  % savefig(fig,[savefolder,savefile])

end