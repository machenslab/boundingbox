function [result,hd,data] = CDGM2021_PlotData(datafolder,plotstring)
  
  tic
  load([datafolder,'summary\datasummary.mat']); % loads structure "data"
  toc
    
  [dimlist, covlist, dellist, boxlist] = deal([]);
  for k = 1:numel(data)
    dimlist = [dimlist, data(k).par.dimension];       %#ok<AGROW>
    covlist = [covlist, data(k).par.coverage];        %#ok<AGROW>
    dellist = [dellist, data(k).par.delay];           %#ok<AGROW>
    boxlist = [boxlist, data(k).par.polytopestretch]; %#ok<AGROW>
  end
  uni.dim = unique(dimlist);
  uni.cov = unique(covlist);
  uni.del = unique(dellist);
  uni.box = unique(boxlist);

  
  [perfplot,dimplot] = deal(NaN(numel(data),numel(uni.box),numel(uni.cov)));
  for k = 1:numel(data)

    dimindex = find(data(k).par.dimension       == uni.dim);
    covindex = find(data(k).par.coverage        == uni.cov);
    delindex = find(data(k).par.delay           == uni.del);
    boxindex = find(data(k).par.polytopestretch == uni.box);

    perfplot(k,dimindex,covindex,delindex,boxindex) = data(k).performance;
    dimplot(k,dimindex,covindex,delindex,boxindex) = data(k).par.dimension;

  end


  quantilelist = [.25; .50; .75];
  
% %   for k = 1:numel(data)
% %     if isfield(data(k),'hon')
% %       cvarflaglist(k) = isfield(data(k).hon,'cv');
% %       rateflaglist(k) = isfield(data(k).hon,'numspike');
% %     else
% %       cvarflaglist(k) = 0;
% %       rateflaglist(k) = 0;
% %     end
% %   end
% %   cvarflag = min(cvarflaglist);
% %   rateflag = min(rateflaglist);
  
  cvarflag = false;
  rateflag = false;

  
  % [linehandle.perf, patchhandle.perf] = deal([]);
  clearvars linehandle.perf patchhandle.perf linehandle.cvar patchhandle.cvar
   
%   variantlist = {'multicoverage','multidelay','multibox'};
%   variantlist = {'multicoverage','multidelay'};
%   variantlist = {'multicoveragelow','multicoveragehigh'};
%   variantlist = {'multicoveragelow','multicoveragehigh','multidim'};
%   variantlist = {'multicoveragelow','multidim'};

  if numel(uni.del) > 1
    variantlist = {'multidim'};  
  else
    variantlist = {'multicoveragehigh'};
  end
  
  for variant = 1:numel(variantlist)
    
    switch variantlist{variant}
      case 'multicoverage'
        displaylist.cov = 1:numel(uni.cov);
        displaylist.del = 1; 
        displaylist.box = NaN;
      case 'multicoveragelow'
        displaylist.cov = 1:numel(uni.cov);
        displaylist.del = 1; 
        displaylist.box = NaN;
      case 'multicoveragehigh'
        displaylist.cov = 1:numel(uni.cov);
        displaylist.del = numel(uni.del); 
        displaylist.box = NaN;
      case 'multidelay'
        displaylist.cov = 1;
        displaylist.del = 1:numel(uni.del);
        displaylist.box = NaN;
      case 'multibox'
        displaylist.cov = 1;
        displaylist.del = 1; 
        displaylist.box = 1:numel(uni.box);
      case 'multidim'
%         displaylist.cov = 1;
        displaylist.cov = ceil(numel(uni.cov)/2);
%         displaylist.del = 1:numel(uni.del); 
        displaylist.dim = 1:numel(uni.dim);
        displaylist.box = NaN;
    end
    
    hd(variant).perfdim.fg = figure();
    if cvarflag, hd(variant).cvardim.fg = figure(); end
    if rateflag, hd(variant).ratedim.fg = figure(); end
    
    switch variantlist{variant}
      case 'multidim'
        [result(variant),hd(variant)] = ...
          makeFigureDel(hd(variant),uni,data,displaylist,quantilelist,cvarflag,rateflag,plotstring);
      otherwise
        [result(variant),hd(variant)] = ...
          makeFigureDim(hd(variant),uni,data,displaylist,quantilelist,cvarflag,rateflag,plotstring);
    end
    
  end
  
end

    
    
function [result,hd] = makeFigureDim(hd,uni,data,displaylist,quantilelist,cvarflag,rateflag,plotstring)

  figure(hd.perfdim.fg)
  hd.perfdim.fg.Color = [1 1 1];
  hd.perfdim.ax(1) = axes();
  hd.perfdim.ax(1).Position = [.07 .13 .90 .80];
  hd.perfdim.ax(1).Title.String = 'performance';
  hold on
  
  if cvarflag
    figure(hd.cvardim.fg)
    hd.cvardim.fg.Color = [1 1 1];
    hd.cvardim.ax(1) = axes();
    hd.cvardim.ax(1).Position = [.07 .13 .90 .80];
    hd.cvardim.ax(1).Title.String = 'coefficient of variation';
    hold on
  end
  
  if rateflag
    figure(hd.ratedim.fg)
    hd.ratedim.fg.Color = [1 1 1];
    hd.ratedim.ax(1) = axes();
    hd.ratedim.ax(1).Position = [.07 .13 .90 .80];
    hd.ratedim.ax(1).Title.String = 'firing rate (Hz)';
    hold on
  end
  
  for covindex = displaylist.cov
    for delindex = displaylist.del
      for boxindex = displaylist.box

% % %         collectgroupdata(covindex,delindex,boxindex);

        groupindex = NaN(numel(uni.dim),numel(data));
        perfgroup = cell(numel(uni.dim),1);
        perfquantile = NaN(numel(quantilelist),numel(uni.dim));
        if cvarflag
          cvargroup = cell(numel(uni.dim),1);
          cvarquantile = NaN(numel(quantilelist),numel(uni.dim));
        end
        if rateflag
          rategroup = cell(numel(uni.dim),1);
          ratequantile = NaN(numel(quantilelist),numel(uni.dim));
        end
        
        for dimindex = 1:numel(uni.dim)
          for k = 1:numel(data)
%             display([covindex,delindex,boxindex])
            if isnan(boxindex)
              groupindex(dimindex,k) = data(k).par.dimension       == uni.dim(dimindex) && ...
                                       data(k).par.coverage        == uni.cov(covindex) && ...
                                       data(k).par.delay           == uni.del(delindex);
            else
              groupindex(dimindex,k) = data(k).par.dimension       == uni.dim(dimindex) && ...
                                       data(k).par.coverage        == uni.cov(covindex) && ...
                                       data(k).par.delay           == uni.del(delindex) && ...
                                       data(k).par.polytopestretch == uni.box(boxindex);
            end
          end
          
          perfgroup{dimindex} = [data(logical(groupindex(dimindex,:))).performance]';
          perfquantile(:,dimindex) = quantile(perfgroup{dimindex}, quantilelist);
%           if uni.dim(dimindex) == 50 && uni.cov(covindex) == 80
%             pause(1)
%           end
          
          if cvarflag
            grouplist = find(groupindex(dimindex,:));
            cvargroup{dimindex} = [];
            for groupmember = 1:numel(grouplist)
              individualcvar = data(grouplist(groupmember)).hon.cv;
              cvargroup{dimindex} = [cvargroup{dimindex}; nanmean(individualcvar)];
%               cvargroup{dimindex} = [cvargroup{dimindex}; individualcvar(:)];
%               cvargroup{dimindex} = [data(logical(groupindex(dimindex,:))).hon.cv]';
            end
            cvarquantile(:,dimindex) = quantile(cvargroup{dimindex}, quantilelist); 
          end
          
          if rateflag
            grouplist = find(groupindex(dimindex,:));
            rategroup{dimindex} = [];
            for groupmember = 1:numel(grouplist)
%               individualrate = data(grouplist(groupmember)).hon.numspike / ...
%                                data(grouplist(groupmember)).par.networksize;
%               individualrate = data(grouplist(groupmember)).hon.numspike / ...
%                                sum(data(grouplist(groupmember)).hon.numspike > 0);
%               individualrate = sum(data(grouplist(groupmember)).hon.excess) / ...
%                                data(grouplist(groupmember)).par.networksize;
                netduration = data(grouplist(groupmember)).par.duration - ...
                              data(grouplist(groupmember)).par.honeymoonduration;
                % population rate:
                individualrate = sum(data(grouplist(groupmember)).hon.numspike) / ...
                                     netduration;
%                 % mean individual neuron rate, all neurons:
%                 individualrate = mean(data(grouplist(groupmember)).hon.numspike) / ...
%                                      netduration;
%                 % population rate, relative to network size (should be same as mean individual):
% %                 individualrate = sum(data(grouplist(groupmember)).hon.numspike) / ...
% %                                      netduration * data(grouplist(groupmember)).par.networksize;
              individualrate(~isfinite(individualrate)) = NaN;
              rategroup{dimindex} = [rategroup{dimindex}; nanmean(individualrate)];
%               rategroup{dimindex} = [rategroup{dimindex}; individualrate(:)];
%               rategroup{dimindex} = [data(logical(groupindex(dimindex,:))).hon.cv]';
            end
            % COMMENT OUT following line if you wish to use RELATIVE rate:
            ratequantile(:,dimindex) = quantile(rategroup{dimindex}, quantilelist);
          end
          
        end
        if isnan(boxindex)
          boxindex = 1;
        end
%         % UNCOMMENT the following loop(s) if you wish to use RELATIVE rate:
%         if rateflag
%           for dimindex = numel(uni.dim):-1:1
%             % OPTIONAL: normalise firing rate to that of first group (e.g. lowest dimensionality):
%             rategroup{dimindex} = rategroup{dimindex} / mean(rategroup{1});
%             ratequantile(:,dimindex) = quantile(rategroup{dimindex}, quantilelist);
%           end
%         end

        
        % for export
        perfquantilecell{covindex,boxindex} = perfquantile;
        if cvarflag, cvarquantilecell{covindex,boxindex} = cvarquantile; end
        if rateflag, ratequantilecell{covindex,boxindex} = ratequantile; end


        placement = 'right';
        switch placement
          case 'left'
            if numel(displaylist.cov)>1
              textypos(covindex,1) = perfquantile(2,1);
              if cvarflag, textypos(covindex,2) = cvarquantile(2,1); end
              if rateflag, textypos(covindex,3) = ratequantile(2,1); end
            elseif numel(displaylist.del)>=1
              textypos(dimindex,1) = perfquantile(2,1);
              if cvarflag, textypos(delindex,2) = cvarquantile(2,1); end
              if rateflag, textypos(delindex,3) = ratequantile(2,1); end
            else
              textypos = NaN;
            end
          case 'right'
            if numel(displaylist.cov)>1
              textypos(covindex,1) = perfquantile(2,end);
              if cvarflag, textypos(covindex,2) = cvarquantile(2,end); end
              if rateflag, textypos(covindex,3) = ratequantile(2,end); end
            elseif numel(displaylist.del)>=1
              textypos(dimindex,1) = perfquantile(2,end);
              if cvarflag, textypos(delindex,2) = cvarquantile(2,end); end
              if rateflag, textypos(delindex,3) = ratequantile(2,end); end
            else
              textypos = NaN;
            end
        end


        count = 0;
        for k = 1:numel(perfgroup)
          count = max(count,numel(perfgroup{k}));
        end
        perfmatrix = NaN(count,numel(perfgroup));
        if cvarflag, cvarmatrix = NaN(count,numel(perfgroup)); end
        for k = 1:numel(perfgroup)
          perfmatrix(1:numel(perfgroup{k}),k) = perfgroup{k};
          if cvarflag, cvarmatrix(1:numel(cvargroup{k}),k) = cvargroup{k}; end
        end

        switch plotstring

          case 'plot'

            figure(hd.perfdim.fg)
            [linehandle.perf(covindex,boxindex), patchhandle.perf(covindex,boxindex)] = ...
              boundedline(uni.dim,perfquantile(2,:),diff(perfquantile,1)');

            if cvarflag
              figure(hd.cvardim.fg)
              [linehandle.cvar(covindex,boxindex), patchhandle.cvar(covindex,boxindex)] = ...
                boundedline(uni.dim,cvarquantile(2,:),diff(cvarquantile,1)');
            end
            
            if rateflag
              figure(hd.ratedim.fg)
              [linehandle.rate(covindex,boxindex), patchhandle.rate(covindex,boxindex)] = ...
                boundedline(uni.dim,ratequantile(2,:),diff(ratequantile,1)');
            end

        end

      end
    end
  end

  switch plotstring

    case 'plot'

      linehandle.perf(~isgraphics(linehandle.perf)) = [];
      patchhandle.perf(~isgraphics(patchhandle.perf)) = [];
      if cvarflag
        linehandle.cvar(~isgraphics(linehandle.cvar)) = [];
        patchhandle.cvar(~isgraphics(patchhandle.cvar)) = [];
      end
      if rateflag
        linehandle.rate(~isgraphics(linehandle.rate)) = [];
        patchhandle.rate(~isgraphics(patchhandle.rate)) = [];
      end
      colourlist = [.0314 .2275 .4863; ...
                    .4 .4 .4; ...
                    .0980 .4118 .6902; ...
                    .4 .4 .4; ...
                    .5216 .7412 .8627];

      colourlist2 = bone(ceil(1.5*size(linehandle.perf,1)));

      darkening = 1.0;  % 1 <= darkening < Inf
      for k = 1:size(colourlist,1)
        [maxvalue,maxposition] = max(colourlist(k,:));
        colourlist(k,:) = colourlist(k,:) / darkening;
        colourlist(k,maxposition) = mean([1 maxvalue]);
      end

    [patchhandle.perf(:).FaceAlpha] = deal(.6);
    hd.perfdim.fg.Children.XScale = 'log';
    for k = 1:size(linehandle.perf,1)
      [linehandle.perf(k,1).Color, patchhandle.perf(k,1).FaceColor] = deal(colourlist(k,:));
      [linehandle.perf(k,2:end).Color, patchhandle.perf(k,2:end).FaceColor] = deal(colourlist2(k,:));
    end
    if cvarflag
      [patchhandle.cvar(:).FaceAlpha] = deal(.6);
        for k = 1:size(linehandle.cvar,1)
          [linehandle.cvar(k,1).Color, patchhandle.cvar(k,1).FaceColor] = deal(colourlist(k,:));
          [linehandle.cvar(k,2:end).Color, patchhandle.cvar(k,2:end).FaceColor] = deal(colourlist2(k,:));
        end
      hd.cvardim.fg.Children.XScale = 'log';
    end
    if rateflag
      [patchhandle.rate(:).FaceAlpha] = deal(.6);
      for k = 1:size(linehandle.rate,1)
        [linehandle.rate(k,1).Color, patchhandle.rate(k,1).FaceColor] = deal(colourlist(k,:));
        [linehandle.rate(k,2:end).Color, patchhandle.rate(k,2:end).FaceColor] = deal(colourlist2(k,:));
      end
      hd.ratedim.fg.Children.XScale = 'log';
    end

      figure(hd.perfdim.fg)
      axes(hd.perfdim.ax(1))

      switch placement
        case 'left'
          textxpos = 1.9;
        case 'right'
          textxpos = 56;
      end

      figure(hd.perfdim.fg)
      if numel(displaylist.cov) > 1
        for covindex = 1:numel(displaylist.cov)
          textstring = '\rho';
          expostyle = returnStringExponent(uni.cov(displaylist.cov(covindex)));
          switch placement
            case 'left'
              text(textxpos,textypos(displaylist.cov(covindex),1), ...
                [textstring,char(8201),'=',char(8201), ...
                num2str(uni.cov(displaylist.cov(covindex)))], ...
                'VerticalAlignment','middle','HorizontalAlignment','right', ...
                'Color',colourlist(covindex,:))
            case 'right'
              text(textxpos,textypos(displaylist.cov(covindex),1), ...
                [textstring,char(8201),'=',char(8201), ...
                num2str(uni.cov(displaylist.cov(covindex)))], ...
                'VerticalAlignment','middle','HorizontalAlignment','left', ...
                'Color',colourlist(covindex,:))
          end
        end
      elseif numel(displaylist.del) >= 1
        for delindex = 1:numel(displaylist.del)
          textstring = '\theta';
          expostyle = returnStringExponent(uni.del(displaylist.del(delindex))*1000);
          text(textxpos,textypos(displaylist.del(delindex),1),[textstring,char(8201),'=',char(8201), ...
	          num2str(uni.del(displaylist.del(delindex))*1000)], ...
            'VerticalAlignment','middle','HorizontalAlignment','right', ...
            'Color',colourlist(delindex,:))
        end
      end
      
      if cvarflag
        figure(hd.cvardim.fg)
        if numel(displaylist.cov) > 1
          for covindex = 1:numel(displaylist.cov)
            textstring = '\rho';
            expostyle = returnStringExponent(uni.cov(displaylist.cov(covindex)));
            text(textxpos,textypos(displaylist.cov(covindex),2),[textstring,char(8201),'=',char(8201), ...
              num2str(uni.cov(displaylist.cov(covindex)))], ...
              'VerticalAlignment','middle','HorizontalAlignment','right', ...
              'Color',colourlist(covindex,:))
          end
        elseif numel(displaylist.del) >= 1
          for delindex = 1:numel(displaylist.del)
            textstring = '\theta';
            expostyle = returnStringExponent(uni.del(displaylist.del(delindex))*1000);
            text(textxpos,textypos(displaylist.del(delindex),2),[textstring,char(8201),'=',char(8201), ...
              num2str(uni.del(displaylist.del(delindex))*1000)], ...
              'VerticalAlignment','middle','HorizontalAlignment','right', ...
              'Color',colourlist(delindex,:))
          end
        end              
      end   
      
      if rateflag
        figure(hd.ratedim.fg)
        if numel(displaylist.cov) > 1
          for covindex = 1:numel(displaylist.cov)
            textstring = '\rho';
            expostyle = returnStringExponent(uni.cov(displaylist.cov(covindex)));
            switch placement
              case 'left'
                text(textxpos,textypos(displaylist.cov(covindex),3), ...
                  [textstring,char(8201),'=',char(8201), ...
                  num2str(uni.cov(displaylist.cov(covindex)))], ...
                  'VerticalAlignment','middle','HorizontalAlignment','right', ...
                  'Color',colourlist(covindex,:))
              case 'right'
                text(textxpos,textypos(displaylist.cov(covindex),3), ...
                  [textstring,char(8201),'=',char(8201), ...
                  num2str(uni.cov(displaylist.cov(covindex)))], ...
                  'VerticalAlignment','middle','HorizontalAlignment','left', ...
                  'Color',colourlist(covindex,:))
            end
          end
        elseif numel(displaylist.del) >= 1
          for delindex = 1:numel(displaylist.del)
            textstring = '\theta';
            expostyle = returnStringExponent(uni.del(displaylist.del(delindex))*1000);
            text(textxpos,textypos(displaylist.del(delindex),3),[textstring,char(8201),'=',char(8201), ...
              num2str(uni.del(displaylist.del(delindex))*1000)], ...
              'VerticalAlignment','middle','HorizontalAlignment','right', ...
              'Color',colourlist(delindex,:))
          end
        end              
      end

      
      xlabel('dimension')
      hd.perfdim.ax(1).XLim = [1.23, max(uni.dim)];
      if cvarflag, hd.cvardim.ax(1).XLim = [1.23, max(uni.dim)]; end
      if rateflag, hd.ratedim.ax(1).XLim = [1.23, max(uni.dim)]; end
      hd.perfdim.ax(1).XTick = [2 5 10 20 50];
      if cvarflag, hd.cvardim.ax(1).XTick = [2 5 10 20 50]; end
      if rateflag, hd.ratedim.ax(1).XTick = [2 5 10 20 50]; end

      % bring to front and arrange
      figure(hd.perfdim.fg);
      if cvarflag, figure(hd.cvardim.fg), end
      if rateflag, figure(hd.ratedim.fg), end
      hd.perfdim.fg.Position = [100 450 560 420];
      if cvarflag, hd.cvardim.fg.Position = [670 450 560 420]; end
      if rateflag, hd.ratedim.fg.Position = [1240 450 560 420]; end

  end

  result.perfquantile = perfquantilecell;
  if cvarflag, result.cvarquantile = cvarquantilecell; end
  if rateflag, result.ratequantile = ratequantilecell; end
  result.uni = uni;

end



function [result,hd] = makeFigureDel(hd,uni,data,displaylist,quantilelist,cvarflag,rateflag,plotstring)

  figure(hd.perfdim.fg)
  hd.perfdim.fg.Color = [1 1 1];
  hd.perfdim.ax(1) = axes();
  hd.perfdim.ax(1).Position = [.07 .13 .90 .80];
  hd.perfdim.ax(1).Title.String = 'performance';
  hold on
  
  if cvarflag
    figure(hd.cvardim.fg)
    hd.cvardim.fg.Color = [1 1 1];
    hd.cvardim.ax(1) = axes();
    hd.cvardim.ax(1).Position = [.07 .13 .90 .80];
    hd.cvardim.ax(1).Title.String = 'coefficient of variation';
    hold on
  end
  
  if rateflag
    figure(hd.ratedim.fg)
    hd.ratedim.fg.Color = [1 1 1];
    hd.ratedim.ax(1) = axes();
    hd.ratedim.ax(1).Position = [.07 .13 .90 .80];
    hd.ratedim.ax(1).Title.String = 'relative firing rate';
    hold on
  end
    
  for covindex = displaylist.cov
    for dimindex = displaylist.dim
      for boxindex = displaylist.box

% % %         collectgroupdata(covindex,dimindex,boxindex);

        groupindex = NaN(numel(uni.del),numel(data));
        perfgroup = cell(numel(uni.del),1);
        perfquantile = NaN(numel(quantilelist),numel(uni.del));
        if cvarflag
          cvargroup = cell(numel(uni.del),1);
          cvarquantile = NaN(numel(quantilelist),numel(uni.del));
        end
        if rateflag
          rategroup = cell(numel(uni.del),1);
          ratequantile = NaN(numel(quantilelist),numel(uni.del));
        end
        
        for delindex = 1:numel(uni.del)
          for k = 1:numel(data) 
            if isnan(boxindex)
              groupindex(delindex,k) = data(k).par.dimension       == uni.dim(dimindex) && ...
                                       data(k).par.coverage        == uni.cov(covindex) && ...
                                       data(k).par.delay           == uni.del(delindex);
            else
              groupindex(delindex,k) = data(k).par.dimension       == uni.dim(dimindex) && ...
                                       data(k).par.coverage        == uni.cov(covindex) && ...
                                       data(k).par.delay           == uni.del(delindex) && ...
                                       data(k).par.polytopestretch == uni.box(boxindex);
            end
          end
          
          perfgroup{delindex} = [data(logical(groupindex(delindex,:))).performance]';
          perfquantile(:,delindex) = quantile(perfgroup{delindex}, quantilelist);
          
          if cvarflag
            grouplist = find(groupindex(delindex,:));
            cvargroup{delindex} = [];
            for groupmember = 1:numel(grouplist)
              individualcvar = data(grouplist(groupmember)).hon.cv;
              cvargroup{delindex} = [cvargroup{delindex}; nanmean(individualcvar)];
%               cvargroup{delindex} = [cvargroup{delindex}; individualcvar(:)];
%               cvargroup{delindex} = [data(logical(groupindex(delindex,:))).hon.cv]';
            end
            cvarquantile(:,delindex) = quantile(cvargroup{delindex}, quantilelist); 
          end
          
          if rateflag
            grouplist = find(groupindex(delindex,:));
            rategroup{delindex} = [];
            for groupmember = 1:numel(grouplist)
              individualrate = data(grouplist(groupmember)).hon.numspike / ...
                               data(grouplist(groupmember)).par.networksize;
%               individualrate = data(grouplist(groupmember)).hon.numspike / ...
%                                sum(data(grouplist(groupmember)).hon.numspike > 0);
%               individualrate = sum(data(grouplist(groupmember)).hon.excess) / ...
%                                data(grouplist(groupmember)).par.networksize;
              netduration = data(grouplist(groupmember)).par.duration - ...
                            data(grouplist(groupmember)).par.honeymoonduration;
              individualrate = sum(data(grouplist(groupmember)).hon.numspike) / ...
                                   netduration;
              individualrate = mean(data(grouplist(groupmember)).hon.numspike) / ...
                                   netduration;
%               individualrate = sum(data(grouplist(groupmember)).hon.numspike) / ...
%                                      netduration * data(grouplist(groupmember)).par.networksize;
              individualrate(~isfinite(individualrate)) = NaN;
              rategroup{delindex} = [rategroup{delindex}; nanmean(individualrate)];
%               rategroup{delindex} = [rategroup{delindex}; individualrate(:)];
%               rategroup{delindex} = [data(logical(groupindex(delindex,:))).hon.cv]';
            end
            % COMMENT OUT following line if you wish to use RELATIVE rate:
            ratequantile(:,delindex) = quantile(rategroup{delindex}, quantilelist);
          end
          
        end
        if isnan(boxindex)
          boxindex = 1;
        end
%         % UNCOMMENT the following loop(s) if you wish to use RELATIVE rate:
%         if rateflag
%           for delindex = numel(uni.del):-1:1
%             % OPTIONAL: normalise firing rate to that of first group:
%             rategroup{delindex} = rategroup{delindex} / mean(rategroup{1});
%             ratequantile(:,delindex) = quantile(rategroup{delindex}, quantilelist);
%           end
%         end
        
        % for export
        perfquantilecell{dimindex,boxindex} = perfquantile;
        if cvarflag, cvarquantilecell{dimindex,boxindex} = cvarquantile; end
        if rateflag, ratequantilecell{dimindex,boxindex} = ratequantile; end


       if numel(displaylist.cov)>1
          textypos(covindex,1) = perfquantile(2,end);
          if cvarflag, textypos(covindex,2) = cvarquantile(2,end); end
          if rateflag, textypos(covindex,3) = ratequantile(2,end); end
        elseif numel(displaylist.dim)>=1
          textypos(dimindex,1) = perfquantile(2,end);
          if cvarflag, textypos(dimindex,2) = cvarquantile(2,end); end
          if rateflag, textypos(dimindex,3) = ratequantile(2,end); end
        else
          textypos = NaN;
        end


        count = 0;
        for k = 1:numel(perfgroup)
          count = max(count,numel(perfgroup{k}));
        end
        perfmatrix = NaN(count,numel(perfgroup));
        if cvarflag, cvarmatrix = NaN(count,numel(perfgroup)); end
        for k = 1:numel(perfgroup)
          perfmatrix(1:numel(perfgroup{k}),k) = perfgroup{k};
          if cvarflag, cvarmatrix(1:numel(cvargroup{k}),k) = cvargroup{k}; end
        end

        switch plotstring

          case 'plot'
%             figure(), boxplot(perfmatrix);
      %       perfmatrix
      %       figure(), boxplot([perfgroup{:}]);

            figure(hd.perfdim.fg)
            [linehandle.perf(dimindex,boxindex), patchhandle.perf(dimindex,boxindex)] = ...
              boundedline(1000*uni.del,perfquantile(2,:),diff(perfquantile,1)');

            if cvarflag
              figure(hd.cvardim.fg)
              [linehandle.cvar(dimindex,boxindex), patchhandle.cvar(dimindex,boxindex)] = ...
                boundedline(1000*uni.del,cvarquantile(2,:),diff(cvarquantile,1)');
            end
            
            if rateflag
              figure(hd.ratedim.fg)
              [linehandle.rate(dimindex,boxindex), patchhandle.rate(dimindex,boxindex)] = ...
                boundedline(1000*uni.del,ratequantile(2,:),diff(ratequantile,1)');
            end

        end

      end
    end
  end

  switch plotstring

    case 'plot'

      linehandle.perf(~isgraphics(linehandle.perf)) = [];
      patchhandle.perf(~isgraphics(patchhandle.perf)) = [];
      if cvarflag
        linehandle.cvar(~isgraphics(linehandle.cvar)) = [];
        patchhandle.cvar(~isgraphics(patchhandle.cvar)) = [];
      end
      if rateflag
        linehandle.rate(~isgraphics(linehandle.rate)) = [];
        patchhandle.rate(~isgraphics(patchhandle.rate)) = [];
      end
      colourlist = winter(2*numel(linehandle.perf));
      
      colourlist(1:3,:) = [.0314 .2275 .4863; ...
                    .0980 .4118 .6902; ...
                    .5216 .7412 .8627];

      colourlist2 = bone(ceil(1.5*size(linehandle.perf,1)));

      darkening = 1.0;  % 1 <= darkening < Inf
      for k = 1:size(colourlist,1)
        [maxvalue,maxposition] = max(colourlist(k,:));
        colourlist(k,:) = colourlist(k,:) / darkening;
        colourlist(k,maxposition) = mean([1 maxvalue]);
      end

    [patchhandle.perf(:).FaceAlpha] = deal(.6);
    hd.perfdim.fg.Children.XScale = 'lin';
    hd.perfdim.fg.Children.YScale = 'lin';
    for k = 1:size(linehandle.perf,1)
      [linehandle.perf(k,1).Color, patchhandle.perf(k,1).FaceColor] = deal(colourlist(k,:));
      [linehandle.perf(k,2:end).Color, patchhandle.perf(k,2:end).FaceColor] = deal(colourlist2(k,:));
    end
    if cvarflag
      [patchhandle.cvar(:).FaceAlpha] = deal(.6);
        for k = 1:size(linehandle.cvar,1)
          [linehandle.cvar(k,1).Color, patchhandle.cvar(k,1).FaceColor] = deal(colourlist(k,:));
          [linehandle.cvar(k,2:end).Color, patchhandle.cvar(k,2:end).FaceColor] = deal(colourlist2(k,:));
        end
      hd.cvardim.fg.Children.XScale = 'lin';
      hd.cvardim.fg.Children.YScale = 'lin';
    end
    if rateflag
      [patchhandle.rate(:).FaceAlpha] = deal(.6);
      for k = 1:size(linehandle.rate,1)
        [linehandle.rate(k,1).Color, patchhandle.rate(k,1).FaceColor] = deal(colourlist(k,:));
        [linehandle.rate(k,2:end).Color, patchhandle.rate(k,2:end).FaceColor] = deal(colourlist2(k,:));
      end
      hd.ratedim.fg.Children.XScale = 'lin';
      hd.ratedim.fg.Children.YScale = 'log';
    end

      figure(hd.perfdim.fg)
      axes(hd.perfdim.ax(1))

      
      textxpos = 5.07;

      figure(hd.perfdim.fg)
      if numel(displaylist.cov) > 1
        for covindex = 1:numel(displaylist.cov)
          textstring = '\rho';
          expostyle = returnStringExponent(uni.cov(displaylist.cov(covindex)));
          text(textxpos,textypos(displaylist.cov(covindex),1),[textstring,char(8201),'=',char(8201), ...
	          num2str(uni.cov(displaylist.cov(covindex)))], ...
            'VerticalAlignment','middle','HorizontalAlignment','right', ...
            'Color',colourlist(covindex,:))
        end
      elseif numel(displaylist.dim) >= 1
        for dimindex = 1:numel(displaylist.dim)
          textstring = 'M';
          expostyle = returnStringExponent(uni.dim(displaylist.dim(dimindex)));
          text(textxpos,textypos(displaylist.dim(dimindex),1),[textstring,char(8201),'=',char(8201), ...
	          num2str(uni.dim(displaylist.dim(dimindex)))], ...
            'VerticalAlignment','middle','HorizontalAlignment','left', ...
            'Color',colourlist(dimindex,:))
        end
      end
      
      if cvarflag
        figure(hd.cvardim.fg)
        if numel(displaylist.cov) > 1
          for covindex = 1:numel(displaylist.cov)
            textstring = '\rho';
            expostyle = returnStringExponent(uni.cov(displaylist.cov(covindex)));
            text(textxpos,textypos(displaylist.cov(covindex),2),[textstring,char(8201),'=',char(8201), ...
              num2str(uni.cov(displaylist.cov(covindex)))], ...
              'VerticalAlignment','middle','HorizontalAlignment','right', ...
              'Color',colourlist(covindex,:))
          end
        elseif numel(displaylist.dim) >= 1
          for dimindex = 1:numel(displaylist.dim)
            textstring = 'M';
            expostyle = returnStringExponent(uni.dim(displaylist.dim(dimindex)));
            text(textxpos,textypos(displaylist.dim(dimindex),2),[textstring,char(8201),'=',char(8201), ...
              num2str(uni.dim(displaylist.dim(dimindex)))], ...
              'VerticalAlignment','middle','HorizontalAlignment','left', ...
              'Color',colourlist(dimindex,:))
          end
        end
      end
      
      if rateflag
        figure(hd.ratedim.fg)
        if numel(displaylist.cov) > 1
          for covindex = 1:numel(displaylist.cov)
            textstring = '\rho';
            expostyle = returnStringExponent(uni.cov(displaylist.cov(covindex)));
            text(textxpos,textypos(displaylist.cov(covindex),3),[textstring,char(8201),'=',char(8201), ...
              num2str(uni.cov(displaylist.cov(covindex)))], ...
              'VerticalAlignment','middle','HorizontalAlignment','right', ...
              'Color',colourlist(covindex,:))
          end
        elseif numel(displaylist.dim) >= 1
          for dimindex = 1:numel(displaylist.dim)
            textstring = 'M';
            expostyle = returnStringExponent(uni.dim(displaylist.dim(dimindex)));
            text(textxpos,textypos(displaylist.dim(dimindex),3),[textstring,char(8201),'=',char(8201), ...
              num2str(uni.dim(displaylist.dim(dimindex)))], ...
              'VerticalAlignment','middle','HorizontalAlignment','left', ...
              'Color',colourlist(dimindex,:))
          end
        end
      end


      xlimit = [0, 1.12*max(uni.del)] * 1000;
      hd.perfdim.ax(1).XLim = xlimit;
      if cvarflag, hd.cvardim.ax(1).XLim = xlimit; end
      if rateflag, hd.ratedim.ax(1).XLim = xlimit; end

      % bring to front and arrange
      figure(hd.perfdim.fg);
      if cvarflag, figure(hd.cvardim.fg), end
      if rateflag, figure(hd.ratedim.fg), end
      hd.perfdim.fg.Position = [100 450 560 420];
      if cvarflag, hd.cvardim.fg.Position = [670 450 560 420]; end
      if rateflag, hd.ratedim.fg.Position = [1240 450 560 420]; end
      
      figure(hd.perfdim.fg);
      axes(hd.perfdim.ax(1));
      hold on
        plot([min(uni.del),max(uni.del)]*1000,[1 1],'--k')
        plot([min(uni.del),max(uni.del)]*1000,[0 0],'--k')
      hold off

      xlabelstring = 'delay length (msec)';
      hd.perfdim.ax(1).XLabel.String = xlabelstring;
      if cvarflag, hd.cvardim.ax(1).XLabel.String = xlabelstring; end
      if rateflag, hd.ratedim.ax(1).XLabel.String = xlabelstring; end

  end

  result.perfquantile = perfquantilecell;
  if cvarflag, result.cvarquantile = cvarquantilecell; end
  if rateflag, result.ratequantile = ratequantilecell; end
  result.uni = uni;

end



function compoundstring = returnStringExponent(numberlist)

  exponent = floor(log10(numberlist));

  factor = numberlist./(10.^exponent);

  % To avoid nasty rounding errors beyond the 5th digit...
  factor = round(factor*1e5)/1e5;

  for k = 1:numel(numberlist)

    if ~isfinite(exponent(k))
      compoundstring{k} = num2str(numberlist(k)); %#ok<AGROW>
    elseif factor(k) == 1
      compoundstring{k} = ['10^{',num2str(exponent(k)),'}']; %#ok<AGROW>
    else
      compoundstring{k} = [num2str(factor(k),'%g\n'),'\cdot','10^{',num2str(exponent(k)),'}']; %#ok<AGROW>
    end

  end
  
end