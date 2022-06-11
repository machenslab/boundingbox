function [result,hd,data] = CDGM2021_PlotSuppNoiseData(datafolder,plotstring)
  
  tic
  load([datafolder,'summary\datasummary.mat']); % loads structure "data"
  toc
    
  [dimlist, covlist, nsdlist, boxlist] = deal([]);
  for k = 1:numel(data)
    dimlist = [dimlist, data(k).par.dimension];       %#ok<AGROW>
    covlist = [covlist, data(k).par.coverage];        %#ok<AGROW>
    nsdlist = [nsdlist, data(k).par.voltage.noisestd]; %#ok<AGROW>
    boxlist = [boxlist, data(k).par.polytopestretch]; %#ok<AGROW>
  end
  uni.dim = unique(dimlist);
  uni.cov = unique(covlist);
  uni.nsd = unique(nsdlist);
  uni.box = unique(boxlist);

  [perfplot,dimplot] = deal(NaN(numel(data),numel(uni.box),numel(uni.cov)));
  for k = 1:numel(data)

    dimindex = find(data(k).par.dimension        == uni.dim);
    covindex = find(data(k).par.coverage         == uni.cov);
    nsdindex = find(data(k).par.voltage.noisestd == uni.nsd);
    boxindex = find(data(k).par.polytopestretch  == uni.box);

    perfplot(k,dimindex,covindex,nsdindex,boxindex) = data(k).performance;
    dimplot(k,dimindex,covindex,nsdindex,boxindex)  = data(k).par.dimension;

  end


  quantilelist = [.25; .50; .75];
  
  for k = 1:numel(data)
    if isfield(data(k),'hon')
      cvarflaglist(k) = isfield(data(k).hon,'cv');
      rateflaglist(k) = isfield(data(k).hon,'numspike');
    else
      cvarflaglist(k) = 0;
      rateflaglist(k) = 0;
    end
  end
  cvarflag = min(cvarflaglist);
  rateflag = min(rateflaglist);

  hd.perfnsd.fg = figure();
  hd.perfnsd.fg.Color = [1 1 1];
  hd.perfnsd.ax(1) = axes();
  hd.perfnsd.ax(1).Position = [.07 .13 .90 .80];
  title('performance')
  hold on
  
  if cvarflag
    hd.cvarnsd.fg = figure();
    hd.cvarnsd.fg.Color = [1 1 1];
    hd.cvarnsd.ax(1) = axes();
    hd.cvarnsd.ax(1).Position = [.07 .13 .90 .80];
    title('coefficient of variation')
    hold on
  end
  
  if rateflag
    hd.ratensd.fg = figure();
    hd.ratensd.fg.Color = [1 1 1];
    hd.ratensd.ax(1) = axes();
    hd.ratensd.ax(1).Position = [.07 .13 .90 .80];
%     title('# excess spikes / neuron')
    title('relative firing rate')
    hold on
  end
  
  disp(uni)
  
  for covindex = 1:numel(uni.cov)
    for dimindex = 1:numel(uni.dim)
      for boxindex = fliplr(1:numel(uni.box))

        groupindex = NaN(numel(uni.nsd),numel(data));
        perfgroup = cell(numel(uni.nsd),1);
        perfquantile = NaN(numel(quantilelist),numel(uni.nsd));
        if cvarflag
          cvargroup = cell(numel(uni.nsd),1);
          cvarquantile = NaN(numel(quantilelist),numel(uni.nsd));
        end
        if rateflag
          rategroup = cell(numel(uni.nsd),1);
          ratequantile = NaN(numel(quantilelist),numel(uni.nsd));
        end
        
        for nsdindex = 1:numel(uni.nsd)
          for k = 1:numel(data)
            groupindex(nsdindex,k) = data(k).par.dimension       == uni.dim(dimindex) && ...
                                     data(k).par.coverage        == uni.cov(covindex) && ...
                                     data(k).par.voltage.noisestd == uni.nsd(nsdindex) && ...
                                     data(k).par.polytopestretch == uni.box(boxindex);
          end

          perfgroup{nsdindex} = [data(logical(groupindex(nsdindex,:))).performance]';
          perfquantile(:,nsdindex) = quantile(perfgroup{nsdindex}, quantilelist);

          if cvarflag
            grouplist = find(groupindex(nsdindex,:));
            cvargroup{nsdindex} = [];
            for groupmember = 1:numel(grouplist)
              individualcvar = data(grouplist(groupmember)).hon.cv;
              cvargroup{nsdindex} = [cvargroup{nsdindex}; mean(individualcvar)];
%               cvargroup{nsdindex} = [cvargroup{nsdindex}; individualcvar(:)];
%               cvargroup{nsdindex} = [data(logical(groupindex(nsdindex,:))).hon.cv]';
            end
            cvarquantile(:,nsdindex) = quantile(cvargroup{nsdindex}, quantilelist);
          end

          if rateflag
            grouplist = find(groupindex(nsdindex,:));
            rategroup{nsdindex} = [];
            for groupmember = 1:numel(grouplist)
%               individualrate = data(grouplist(groupmember)).hon.numspike / ...
%                                data(grouplist(groupmember)).par.networksize;
%               individualrate = data(grouplist(groupmember)).hon.numspike / ...
%                                sum(data(grouplist(groupmember)).hon.numspike > 0);
%               individualrate = sum(data(grouplist(groupmember)).hon.excess) / ...
%                                data(grouplist(groupmember)).par.networksize;
                netduration = data(grouplist(groupmember)).par.duration - ...
                              data(grouplist(groupmember)).par.honeymoonduration;
                individualrate = sum(data(grouplist(groupmember)).hon.numspike) / ...
                                     netduration;
%                 individualrate = sum(data(grouplist(groupmember)).hon.numspike) / ...
%                                      netduration * data(grouplist(groupmember)).par.networksize;
              individualrate(~isfinite(individualrate)) = NaN;
              rategroup{nsdindex} = [rategroup{nsdindex}; mean(individualrate)];
%               rategroup{nsdindex} = [rategroup{nsdindex}; individualrate(:)];
%               rategroup{nsdindex} = [data(logical(groupindex(nsdindex,:))).hon.cv]';
            end
            % COMMENT OUT following line if you wish to use RELATIVE rate:
            ratequantile(:,nsdindex) = quantile(rategroup{nsdindex}, quantilelist);
          end
          
        end
%         % UNCOMMENT the following loop(s) if you wish to use RELATIVE rate:
%         if rateflag
%           for nsdindex = numel(uni.nsd):-1:1
%             % OPTIONAL: normalise firing rate to that of first group:
%             rategroup{nsdindex} = rategroup{nsdindex} / mean(rategroup{1});
%             ratequantile(:,nsdindex) = quantile(rategroup{nsdindex}, quantilelist);
%           end
%         end
        
        % for export
        perfquantilecell{covindex,boxindex,dimindex} = perfquantile;
        if cvarflag, cvarquantilecell{covindex,boxindex,dimindex} = cvarquantile; end
        if rateflag, ratequantilecell{covindex,boxindex,dimindex} = ratequantile; end

           
        textypos(covindex) = perfquantile(2,end);

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
          
            figure(), boxplot(perfmatrix);

            figure(hd.perfnsd.fg)
            [linehandle.perf(covindex,boxindex), patchhandle.perf(covindex,boxindex)] = ...
              boundedline(uni.nsd,perfquantile(2,:),diff(perfquantile,1)'); %#ok<SAGROW>

            if cvarflag
              figure(hd.cvarnsd.fg)
              [linehandle.cvar(covindex,boxindex), patchhandle.cvar(covindex,boxindex)] = ...
                boundedline(uni.nsd,cvarquantile(2,:),diff(cvarquantile,1)'); %#ok<SAGROW>
            end

            if rateflag
              figure(hd.ratensd.fg)
              [linehandle.rate(covindex,boxindex), patchhandle.rate(covindex,boxindex)] = ...
                boundedline(uni.nsd,ratequantile(2,:),diff(ratequantile,1)'); %#ok<SAGROW>
            end
          
        end

      end
    end
  end

  switch plotstring
    case 'plot'
          
    colourlist = winter(ceil(1.5*size(linehandle.perf,1)));
    colourlist2 = bone(ceil(1.5*size(linehandle.perf,1)));

    darkening = 1.0;  % 1 <= darkening < Inf
    for k = 1:size(colourlist,1)
      [maxvalue,maxposition] = max(colourlist(k,:));
      colourlist(k,:) = colourlist(k,:) / darkening;
      colourlist(k,maxposition) = mean([1 maxvalue]);
    end

    [patchhandle.perf(:).FaceAlpha] = deal(.5);
    for k = 1:size(linehandle.perf,1)
      [linehandle.perf(k,1).Color, patchhandle.perf(k,1).FaceColor] = deal(colourlist(k,:));
      [linehandle.perf(k,2:end).Color, patchhandle.perf(k,2:end).FaceColor] = deal(colourlist2(k,:));
    end
    if cvarflag
      [patchhandle.cvar(:).FaceAlpha] = deal(.5);
        for k = 1:size(linehandle.cvar,1)
          [linehandle.cvar(k,1).Color, patchhandle.cvar(k,1).FaceColor] = deal(colourlist(k,:));
          [linehandle.cvar(k,2:end).Color, patchhandle.cvar(k,2:end).FaceColor] = deal(colourlist2(k,:));
        end
    end
    if rateflag
      [patchhandle.rate(:).FaceAlpha] = deal(.5);
      for k = 1:size(linehandle.rate,1)
        [linehandle.rate(k,1).Color, patchhandle.rate(k,1).FaceColor] = deal(colourlist(k,:));
        [linehandle.rate(k,2:end).Color, patchhandle.rate(k,2:end).FaceColor] = deal(colourlist2(k,:));
      end
    end

    figure(hd.perfnsd.fg)
    axes(hd.perfnsd.ax(1))

    textxpos = max(uni.nsd)*1.1;

    for covindex = 1:numel(uni.cov)
      textstring = 'C';
      text(textxpos,textypos(covindex),[textstring,' = ',num2str(uni.cov(covindex))],'VerticalAlignment','middle','HorizontalAlignment','right')
    end

    xlabel('noise standard deviation')
    [hd.perfnsd.ax(:).XLim] = deal([0, max([max(hd.perfnsd.ax(1).XLim), max(textxpos)])]);
    
    % bring to front and arrange
    figure(hd.perfnsd.fg);
    if cvarflag, figure(hd.cvarnsd.fg), end
    if rateflag, figure(hd.ratensd.fg), end
    hd.perfnsd.fg.Position = [100 450 560 420];
    if cvarflag, hd.cvarnsd.fg.Position = [670 450 560 420]; end
    if rateflag, hd.ratensd.fg.Position = [1240 450 560 420]; end
    
  end

  result.perfquantile = perfquantilecell;
  if cvarflag, result.cvarquantile = cvarquantilecell; end
  if rateflag, result.ratequantile = ratequantilecell; end
  result.uni = uni;
  
end