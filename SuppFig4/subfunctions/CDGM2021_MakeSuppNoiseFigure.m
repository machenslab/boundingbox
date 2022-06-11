
clearvars -except result data
fontsize = 12;

hd.noisesupp.fg = figure();
figurelist = findobj(allchild(groot), 'flat', 'type', 'figure');

hd.noisesupp.fg.Color = [1 1 1];
hd.noisesupp.fg.Position = [0 80 715 775];


% facilitate logarithmic plotting (so "x=0" can be shown at 1e-9, manually labelled 0)
for k = 1:numel(result.perfquantile)
  result.perfquantile{k}(result.perfquantile{k}<0) = 1e-9;
  result.cvarquantile{k}(result.cvarquantile{k}<0) = 1e-9;
  result.ratequantile{k}(result.ratequantile{k}<0) = 1e-9;
end

panel.xorigin = 45;
panel.yorigin = 620;
panel.width = 170;
panel.height = 110;
panel.xoffset = 70;
panel.yoffset = 30;

for row = 1:5
  for column = 1:3
    
    % create axes object
    hd.noisesupp.ax(row,column) = axes();
    hd.noisesupp.ax(row,column).Units = 'pixels';
    hd.noisesupp.ax(row,column).FontSize = fontsize;
    hd.noisesupp.ax(row,column).Position = ...
      [panel.xorigin(1) + (column-1)*panel.width(1) + (column-1)*panel.xoffset(1), ...
       panel.yorigin(1) - (row-1)*panel.height(1) - (row-1)*panel.yoffset(1), ...
       panel.width(1), panel.height(1)];
     
    % create boundedline plots of results
    for boxindex = 1:size(result.perfquantile,2)
      for covindex = 1:size(result.perfquantile,1)
        for dimindex = row
          hold on
          if column == 1
            [linehandle.perf(covindex,boxindex,dimindex), ...
             patchhandle.perf(covindex,boxindex,dimindex)] = ...
              boundedline(result.uni.nsd, result.perfquantile{covindex,boxindex,dimindex}(2,:),...
                          diff(result.perfquantile{covindex,boxindex,dimindex},1)');
          elseif column == 3
            [linehandle.cvar(covindex,boxindex,dimindex), ...
             patchhandle.cvar(covindex,boxindex,dimindex)] = ...
                boundedline(result.uni.nsd, result.cvarquantile{covindex,boxindex,dimindex}(2,:), ...
                            diff(result.cvarquantile{covindex,boxindex,dimindex},1)');
          elseif column == 2
    %         relativespikenum = result.ratequantile{covindex,boxindex,dimindex}(2,:); % = spikenum / networksize
    %         relativespikenumpersec = relativespikenum / data(1).par.duration;
    %         actualexcess = (relativespikenumpersec-relativespikenumpersec(1))/relativespikenumpersec(1);
    %         [linehandle.rate(covindex,boxindex,dimindex), patchhandle.rate(covindex,boxindex,dimindex)] = ...
    %             boundedline(result.uni.nsd, actualexcess, ...
    %                                         diff(result.ratequantile{covindex,boxindex,dimindex},1)');
            [linehandle.rate(covindex,boxindex,dimindex), ...
             patchhandle.rate(covindex,boxindex,dimindex)] = ...
                boundedline(result.uni.nsd, result.ratequantile{covindex,boxindex,dimindex}(2,:), ...
                            diff(result.ratequantile{covindex,boxindex,dimindex},1)');
          end
        end
        hold off
      end
    end
    if column == 2
  %     hd.noisesupp.ax(row,column).YLabel.String = '# excess spikes/neuron';
      hd.noisesupp.ax(row,column).YLabel.String = 'population firing rate';
      hd.noisesupp.ax(row,column).YScale = 'log';
%       hd.noisesupp.ax(row,column).YLim(1) = 40;
    elseif column == 3
      hd.noisesupp.ax(row,column).YLabel.String = 'coefficient of variation';
      hd.noisesupp.ax(row,column).YLim(1) = 0;
    elseif column == 1
      hd.noisesupp.ax(row,column).YLabel.String = 'performance';
      hd.noisesupp.ax(row,column).YScale = 'lin';
      hd.noisesupp.ax(row,column).YLim = [0 1.1];
%       hd.noisesupp.ax(row,column).YLim = [.9 1.02];
%       hd.noisesupp.ax(row,column).YLim = [.9 1.02];
    end  
    hd.noisesupp.ax(row,column).Box = 'off';
    hd.noisesupp.ax(row,column).YGrid = 'off';
    hd.noisesupp.ax(row,column).XScale = 'log';
    if row == 5
      hd.noisesupp.ax(row,column).XLabel.String = 'noise standard deviation';
    else
      hd.noisesupp.ax(row,column).XLabel.Visible = 'off';
      hd.noisesupp.ax(row,column).XTickLabel = [];
    end
    if row == 1
      hd.noisesupp.ax(row,column).Title.String = hd.noisesupp.ax(row,column).YLabel.String;
      hd.noisesupp.ax(row,column).Title.FontSize = hd.noisesupp.ax(row,column).YLabel.FontSize;
      hd.noisesupp.ax(row,column).Title.FontWeight = hd.noisesupp.ax(row,column).YLabel.FontWeight;
      hd.noisesupp.ax(row,column).Title.Units = 'pixels';
      hd.noisesupp.ax(row,column).Title.Position(2) = ...
        hd.noisesupp.ax(row,column).Title.Position(2) + 10;
    end
    if column == 1
      hd.noisesupp.ax(row,column).YLabel.String = [num2str(result.uni.dim(row)),'D'];
      hd.noisesupp.ax(row,column).YLabel.String = ['M = ',num2str(result.uni.dim(row))];
    else
      hd.noisesupp.ax(row,column).YLabel.Visible = 'off';
    end
    
    
    
    % add text label to each line:
    textypos = NaN(numel(result.uni.cov),numel(result.uni.box),5);
    for boxindex = 1:numel(result.uni.box)
      dimindex = row;
      switch column
        case 3
          for covindex = 1:size(result.cvarquantile,1)
            if ~(boxindex == 2 && covindex < 3)
              textypos(covindex,boxindex,dimindex) = result.cvarquantile{covindex,boxindex,dimindex}(2,end);
            end
          end
        case 2
          for covindex = 1:numel(result.ratequantile)
    %         textypos(k) = result.ratequantile{k}(2,end);
            try
              textypos(covindex,boxindex,dimindex) = linehandle.rate(covindex,boxindex,dimindex).YData(end);
            catch
            end
          end
        case 1
          for covindex = 1:size(result.perfquantile,1)
            if ~(boxindex == 2 && covindex < 3)
              textypos(covindex,boxindex,dimindex) = result.perfquantile{covindex,boxindex,dimindex}(2,end);
            end
          end
      end
    end
    textxpos = max(result.uni.nsd)*1.1;
    for boxindex = 1:numel(result.uni.box)
      for covindex = 1:numel(result.uni.cov)
        dimindex = row;
%         textstring = 'C';
        textstring = '\rho';
        if boxindex > 1
%           boxstring = 'w';
          boxstring = '';
        else
          boxstring = '';
        end
        try
          linelabel(covindex,boxindex,dimindex) = text(textxpos,textypos(covindex,boxindex,dimindex), ...
            [textstring,'=',num2str(result.uni.cov(covindex)),boxstring], ...
            'VerticalAlignment','middle','HorizontalAlignment','left');
        catch
        end
      end
    end

    hd.noisesupp.ax(row,column).XLim = [0, max([max(hd.noisesupp.ax(row,column).XLim), max(textxpos)*3])];

    % move 0 to small positive value for logarithmic plotting
    zerolocation = 1e-2;
    kids = hd.noisesupp.ax(row,column).Children;
    for k = 7:2:numel(kids)
      kids(k).XData(1) = zerolocation;
    end
    for k = 8:2:numel(kids)
      kids(k).Vertices([1 end],1) = zerolocation*[1;1];
  %     kids(k).XData(end) = kids(k).XData(1);
    end
    
  end
end



% set plot colours
colourlist = [1 .5 .2; ...
              .2 .8 .4; ...
              .4 .6 1];
colourlist = flipud(colourlist);
% colourlist = flipud(cbrewer('seq','Greens',ceil(1.5*size(linehandle.perf,1))));
colourlist2 = flipud(cbrewer('seq','Greys',ceil(1.5*size(linehandle.perf,1))));
% colourlist2 = flipud(colourlist2);
% darkening = 1.0;  % 1 <= darkening < Inf
% for covindex = 1:size(colourlist,1)
%   [maxvalue,maxposition] = max(colourlist(covindex,:));
%   colourlist(covindex,:) = colourlist(covindex,:) / darkening;
%   colourlist(covindex,maxposition) = mean([1 maxvalue]);
% end
[patchhandle.perf(:).FaceAlpha] = deal(.3);
% dimindex = row;
for covindex = 1:size(linehandle.perf,1)
  [linehandle.perf(covindex,1,:).Color, ...
   patchhandle.perf(covindex,1,:).FaceColor] = deal(colourlist(covindex,:));
  [linehandle.perf(covindex,2:end,:).Color, ...
   patchhandle.perf(covindex,2:end,:).FaceColor] = deal(colourlist2(covindex,:));
%   try linelabel(covindex,boxindex,3).Color = colourlist2(covindex,:); end
end
[patchhandle.cvar(:).FaceAlpha] = deal(.3);
for covindex = 1:size(linehandle.cvar,1)
  [linehandle.cvar(covindex,1,:).Color, ...
   patchhandle.cvar(covindex,1,:).FaceColor] = deal(colourlist(covindex,:));
  [linehandle.cvar(covindex,2:end,:).Color, ...
   patchhandle.cvar(covindex,2:end,:).FaceColor] = deal(colourlist2(covindex,:));
%   try linelabel(covindex,boxindex,1).Color = colourlist2(covindex,:); end
end
[patchhandle.rate(:).FaceAlpha] = deal(.3);
for covindex = 1:size(linehandle.rate,1)
  [linehandle.rate(covindex,1,:).Color, ...
   patchhandle.rate(covindex,1,:).FaceColor] = deal(colourlist(covindex,:));
  [linehandle.rate(covindex,2:end,:).Color, ...
   patchhandle.rate(covindex,2:end,:).FaceColor] = deal(colourlist2(covindex,:));
%   try linelabel(covindex,boxindex,2).Color = colourlist2(covindex,:); end
end




% add an explicit legend onto one (or several) panel(s)
for row = 1:5
  for column = 1:3
    hd.noisesupp.ax(row,column).XLim(1) = min(hd.noisesupp.ax(row,column).Children(end).XData);
    hd.noisesupp.ax(row,column).XLim(2) = max(hd.noisesupp.ax(row,column).Children(end).XData);
  end
end
% ylist = .5:-.1:.1;
ylist = [.65 .50 .35 .20];
% clist = [3,2,1,4,5,6];
clist = [3,2,1,4];
count = 0;
for k = [1 4 5 6]
  count = count+1;
%   legendcolour(k,:) = hd.noisesupp.ax(1,1).Children(5-k).Color;
  legendstring{count} = hd.noisesupp.ax(1,1).Children(k).String;
  if k < 4
%   if k > 3
    legendstring{count} = [legendstring{count},' (wide)'];
  end
end
axes(hd.noisesupp.ax(1,1))
% hd.noisesupp.ax(1,1).YLim = [0 6.5];
% hd.noisesupp.ax(1,1).YLim = [0 ];
hold on
% for k = 1:numel(legendstring)
count = 0;
% for k = [1, 4:6]
for k = 1:4
  count = count+1;
  t = text(0,0,legendstring{k});
  t.FontSize = fontsize;
  t.Units = 'normalized';
  switch hd.noisesupp.ax(1,1).XScale
    case {'lin','linear'}
      t.Position = [.3, ylist(count), 0];
    case {'log','logarithmic'}
      t.Position = [.2, ylist(count), 0];
  end
%   t.Color = legendcolour(k,:);
  t.VerticalAlignment = 'middle';
  
  xlimit = hd.noisesupp.ax(1,1).XLim;
  ylimit = hd.noisesupp.ax(1,1).YLim;
  if ~isempty(zerolocation)
    xlimit(1) = zerolocation;
  end
  t.Units = 'Data';
  switch hd.noisesupp.ax(1,1).XScale
    case {'lin','linear'}
      p = plot(diff(xlimit)*[.1 .25], diff(ylimit)*ylist(count)*[1 1]);
    case {'log','logarithmic'}
      p = plot([1.7*xlimit(1) .9*t.Position(1)], ylimit(1) + diff(ylimit)*ylist(count)*[1 1]);
%       p = plot([.105 .15], diff(ylimit)*ylist(count)*[1 1]);
% %       p = plot([.12 .18], diff(ylimit)*ylist(count)*[1 1]);
  end
  p.LineWidth = 2;
  p.LineStyle = '-';
%   p.Color = t.Color;
end
hold off



% remove classic labels next to each line (comment out for debugging)
for row = 1:5
  for column = 1:3
    if column == 1 && row == 1
      [hd.noisesupp.ax(row,column).Children(9:14).Visible] = deal('off');
    else
      [hd.noisesupp.ax(row,column).Children(1:6).Visible] = deal('off');
    end
  end
end

% colour legend as desired
count = 0;
% for k = [6, 4,2,1]
for k = [1:3 6]
  count = count+1;
  [hd.noisesupp.ax(1,1).Children(2*count-1).Color, ...
   hd.noisesupp.ax(1,1).Children(2*count).Color] = deal(hd.noisesupp.ax(1,1).Children(end+1-2*k).Color);
%   [hd.noisesupp.ax(1,1).Children(2*count-1).Color, ...
%    hd.noisesupp.ax(1,1).Children(2*count).Color] = deal(hd.noisesupp.ax(1,1).Children(end-12+2*k).FaceColor);
%   hd.noisesupp.ax(1,1).Children(2*k-1).Color = hd.noisesupp.ax(1,1).Children(end-9+2*k).Color;
end




% add markers for split x axis
if ~isempty(zerolocation)
  x1 = .045;
  dx = .02;
  x2 = .07;
  dy = .3;
  for row = 1:5
    for column = 1:3
      hd.noisesupp.ax(row,column,2) = axes;
      hd.noisesupp.ax(row,column,2).Units = hd.noisesupp.ax(row,column).Units;
      hd.noisesupp.ax(row,column,2).Position([1 3]) = hd.noisesupp.ax(row,column).Position([1 3]);
      hd.noisesupp.ax(row,column,2).Position(4) = 30;
      hd.noisesupp.ax(row,column,2).Position(2) = hd.noisesupp.ax(row,column).Position(2) ...
                                          - hd.noisesupp.ax(row,column,2).Position(4)/2;
      hd.noisesupp.ax(row,column,2).FontSize = hd.noisesupp.ax(row,column).FontSize;
      hd.noisesupp.ax(row,column,2).Color = 'none';
      [hd.noisesupp.ax(row,column,2).XAxis.Visible, ...
       hd.noisesupp.ax(row,column,2).YAxis.Visible] = deal('off');
      hold on
      plot([x1+dx/2 x2+dx/2],[0 0],'-w','Linewidth',2)
      plot([x1,x1+dx],[-dy,dy],'-k')
      plot([x2,x2+dx],[-dy,dy],'-k')
      if row == 5
        text(0,-2.32*dy,'0','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fontsize)
      end
      hold off
      hd.noisesupp.ax(row,column,2).XLim = [0 1];
      hd.noisesupp.ax(row,column,2).YLim = [-1 1];
      hd.noisesupp.ax(row,column).TickLength = [.04 0];
      hd.noisesupp.ax(row,column).XTick = [.1 1];
    end
  end
end



% impose same Y axis labels on all panels of a column
for column = 1:3
  columnlimit = [min([hd.noisesupp.ax(:,column,1).YLim]), max([hd.noisesupp.ax(:,column,1).YLim])];
%   columnkids = findobj(hd.noisesupp.ax(:,column,1),'-not','Type','Text','-not','Type','Axes');
%   columndata = [];
%   for kid = 1:numel(columnkids)
%     columndata = [columndata; columnkids(kid).YData(:)];
%   end
%   columnlimit = [min(columndata), max(columndata)];
  for row = 1:5
    hd.noisesupp.ax(row,column).YLim = columnlimit;
  end
end

for row = 1:5
  hd.noisesupp.ax(row,1,1).YTick = [hd.noisesupp.ax(row,1).YTick(1), 1];
end



% % HACK: manual correction
% % Oddly, the data for "low coverage narrow box" and "high coverage wide box"
% % appear to have traded places. This is likely a bug from the version
% % analyseServerFiles being used. For now, swapping them manually (fix actual
% % problem eventually!)
% for row = 1:5
%   for column = 1:3
% %     if row == 1 && column == 1
% %       
% %     else
%     indexA = 1;
%     indexB = 11;
%     colourA = hd.noisesupp.ax(row,column,1).Children(end-indexA).Color;
%     colourB = hd.noisesupp.ax(row,column,1).Children(end-indexB).Color;
%     [hd.noisesupp.ax(row,column,1).Children(end-indexA+1).FaceColor, ...
%      hd.noisesupp.ax(row,column,1).Children(end-indexA).Color] = deal(colourB);
%     [hd.noisesupp.ax(row,column,1).Children(end-indexB+1).FaceColor, ...
%      hd.noisesupp.ax(row,column,1).Children(end-indexB).Color] = deal(colourA);
%   end
% end

compact = false;

if compact

  % upon request by CM, cut out two rows (no. 2 and no. 4) for the paper supplement
  for row = [2 4]
    for column = 1:3
      for layer = 1:2
        hd.noisesupp.ax(row,column,layer).Visible = 'off';
        kids = hd.noisesupp.ax(row,column,layer).Children;
        for kid = 1:numel(kids)
          kids(kid).Visible = 'off';
        end
      end
    end
  end
  for column = 1:3
    for layer = 1:2
      hd.noisesupp.ax(1,column,layer).Position(2) = hd.noisesupp.ax(1,column).Position(2) ...
                                                    - 2*panel.height(1) - 2*panel.yoffset(1);
      hd.noisesupp.ax(3,column,layer).Position(2) = hd.noisesupp.ax(3,column).Position(2) ...
                                                    - 1*panel.height(1) - 1*panel.yoffset(1);
    end
  end
  hd.noisesupp.fg.Position(4) = hd.noisesupp.fg.Position(4) - 2*panel.height(1) - 2*panel.yoffset(1);

  
  % FD20200510: OPTION A, remove first row
  firstrow  = hd.noisesupp.ax(1,:,:);
  secondrow = hd.noisesupp.ax(3,:,:);
  % replicate titles of first row on second row
  for a = 1:numel(firstrow)
    secondrow(a).Title = firstrow(a).Title;
    secondrow(a).Title.Visible = 'on';
    secondrow(a).Title.Position(2) = ceil(secondrow(a).Title.Position(2) * 1.02); % minor adjustment
  end
  % replicate line labels of first row on second row
  copyobj(firstrow(1).Children(1:8), secondrow(1));
  % hide first row
  [firstrow.Visible] = deal('off');
  for a = 1:numel(firstrow)
    kids = firstrow(a).Children;
    [kids.Visible] = deal('off');
  end
  % reduce figure height
  hd.noisesupp.fg.Position(4) = 340;
  
else

  % FD20200510: OPTION B, remove first row IFF rows 2 and 4 were not removed after all
  firstrow  = hd.noisesupp.ax(1,:,:);
  secondrow = hd.noisesupp.ax(2,:,:);
  % replicate titles of first row on second row
  for a = 1:numel(firstrow)
    secondrow(a).Title = firstrow(a).Title;
    secondrow(a).Title.Visible = 'on';
    secondrow(a).Title.Position(2) = ceil(secondrow(a).Title.Position(2) * 1.04); % minor adjustment
  end
  % replicate line labels of first row on second row
  copyobj(firstrow(1).Children(1:8), secondrow(1));
  % hide first row
  [firstrow.Visible] = deal('off');
  for a = 1:numel(firstrow)
    kids = firstrow(a).Children;
    [kids.Visible] = deal('off');
  end
  % reduce figure height
  hd.noisesupp.fg.Position(4) = 620;
  
end


figurelist = findobj(allchild(groot), 'flat', 'type', 'figure');
for k = 1:numel(figurelist)-1
  close(k)
end