partialfigurelist = findobj(allchild(groot), 'flat', 'type', 'figure');
hd.noise.fg = figure();
figurelist = findobj(allchild(groot), 'flat', 'type', 'figure');

hd.noise.fg.Color = [1 1 1];
hd.noise.fg.Position = [50 80 770 250];
hd.noise.fg.Position = [50 80 515 250];

cmap = parula(22);
colormap(cmap(1:2:19,:))

panel.width = 173;
panel.height = 150;
panel.xorigin = 75;
panel.yorigin = 62;
panel.xoffset = 75;

fontsize = 11;

for letter = 1:3
  hd.noise.lt(letter) = axes();
  hd.noise.lt(letter).Units = 'pixels';
  axis off
  switch letter
    case 1
      letterstring = 'E';
    case 2
      letterstring = 'F';
    case 3
      letterstring = 'X';
  end
  text(1,1,letterstring,'FontWeight','normal','FontSize',fontsize*1.2)
%   text(1,1,char(-32+96+letter),'FontWeight','normal','FontSize',fontsize*1.2)
end

% to facilitate logarithmic plotting
for k = 1:numel(result.perfquantile)
  result.perfquantile{k}(result.perfquantile{k}<0) = 1e-9;
  result.cvarquantile{k}(result.cvarquantile{k}<0) = 1e-9;
  result.ratequantile{k}(result.ratequantile{k}<0) = 1e-9;
end
            
for column = 1:3
  % panels bottom row (numbered "6" to allow insertion of extra rows): evaluation results
  hd.noise.ax(6,column) = axes();
  hd.noise.ax(6,column).Units = 'pixels';
  hd.noise.ax(6,column).FontSize = fontsize;
  hd.noise.ax(6,column).Position = [panel.xorigin + (column-1)*panel.width + (column-1)*panel.xoffset, ...
                                    panel.yorigin, panel.width, panel.height];
  for covindex = 1:size(result.perfquantile,1)
    for boxindex = 1:size(result.perfquantile,2)
      hold on
      if column == 3
        [linehandle.perf(covindex,boxindex), patchhandle.perf(covindex,boxindex)] = ...
          boundedline(result.uni.nsd, result.perfquantile{covindex,boxindex}(2,:),...
                                      diff(result.perfquantile{covindex,boxindex},1)');
      elseif column == 1
        [linehandle.cvar(covindex,boxindex), patchhandle.cvar(covindex,boxindex)] = ...
            boundedline(result.uni.nsd, result.cvarquantile{covindex,boxindex}(2,:), ...
                                        diff(result.cvarquantile{covindex,boxindex},1)');
      elseif column == 2
%         relativespikenum = result.ratequantile{covindex,boxindex}(2,:); % = spikenum / networksize
%         relativespikenumpersec = relativespikenum / data(1).par.duration;
%         actualexcess = (relativespikenumpersec-relativespikenumpersec(1))/relativespikenumpersec(1);
%         [linehandle.rate(covindex,boxindex), patchhandle.rate(covindex,boxindex)] = ...
%             boundedline(result.uni.nsd, actualexcess, ...
%                                         diff(result.ratequantile{covindex,boxindex},1)');
        [linehandle.rate(covindex,boxindex), patchhandle.rate(covindex,boxindex)] = ...
            boundedline(result.uni.nsd, result.ratequantile{covindex,boxindex}(2,:), ...
                                        diff(result.ratequantile{covindex,boxindex},1)');

      end
      hold off
    end
  end
  if column == 2
%     hd.noise.ax(6,column).YLabel.String = '# excess spikes/neuron';
    hd.noise.ax(6,column).YLabel.String = 'population firing rate';
    hd.noise.ax(6,column).YScale = 'log';
    hd.noise.ax(6,column).YLim(1) = 40;
    
  elseif column == 1
    hd.noise.ax(6,column).YLabel.String = 'coefficient of variation';
    hd.noise.ax(6,column).YLim(1) = 0;
    
  elseif column == 3
    hd.noise.ax(6,column).YLabel.String = 'performance';
    hd.noise.ax(6,column).YScale = 'lin';
    hd.noise.ax(6,column).YLim = [0 1.1];
    
  end  
  hd.noise.ax(6,column).XLabel.String = 'noise standard deviation';
  hd.noise.ax(6,column).Box = 'off';
  hd.noise.ax(6,column).YGrid = 'off';
  hd.noise.ax(6,column).XScale = 'log';
  
  xlabel('noise standard deviation')

  textypos = NaN(numel(result.uni.cov),numel(result.uni.box));
  for boxindex = 1:numel(result.uni.box)
    switch column
      case 1
        for covindex = 1:size(result.cvarquantile,1)
          if ~(boxindex == 2 && covindex < 3)
            textypos(covindex,boxindex) = result.cvarquantile{covindex,boxindex}(2,end);
          end
        end
      case 2
        for covindex = 1:numel(result.ratequantile)
          try
            textypos(covindex,boxindex) = linehandle.rate(covindex,boxindex).YData(end);
          catch
          end
        end
      case 3
        for covindex = 1:size(result.perfquantile,1)
          if ~(boxindex == 2 && covindex < 3)
            textypos(covindex,boxindex) = result.perfquantile{covindex,boxindex}(2,end);
          end
        end
    end
  end
    textxpos = max(result.uni.nsd)*1.1;

     for boxindex = 1:numel(result.uni.box)
      for covindex = 1:numel(result.uni.cov)
        textstring = 'C';
        try
          linelabel(covindex,boxindex,column) = text(textxpos,textypos(covindex,boxindex), ...
            [num2str(result.uni.cov(covindex)*2),' neurons'], ...
            'VerticalAlignment','middle','HorizontalAlignment','left');
        catch
        end
      end
    end

    hd.noise.ax(6,column).XLim = [0, max([max(hd.noise.ax(6,column).XLim), max(textxpos)*3])];
    
    % move 0 to small positive value for logarithmic plotting
    zerolocation = 1e-2;
    kids = hd.noise.ax(6,column).Children;
    for k = 7:2:numel(kids)
      kids(k).XData(1) = zerolocation;
    end
    for k = 8:2:numel(kids)
      kids(k).Vertices([1 end],1) = zerolocation*[1;1];
    end

end

colourlist2 = bone(ceil(1.5*size(linehandle.perf,1)));

colourlist = [1 .5 .2; ...
              .2 .8 .4; ...
              .4 .6 1];
colourlist = flipud(colourlist);


[patchhandle.perf(:).FaceAlpha] = deal(.3);
for covindex = 1:size(linehandle.perf,1)
  [linehandle.perf(covindex,1).Color, patchhandle.perf(covindex,1).FaceColor] = deal(colourlist(covindex,:));
  [linehandle.perf(covindex,2:end).Color, patchhandle.perf(covindex,2:end).FaceColor] = deal(colourlist2(covindex,:));
end
[patchhandle.cvar(:).FaceAlpha] = deal(.3);
for covindex = 1:size(linehandle.cvar,1)
  [linehandle.cvar(covindex,1).Color, patchhandle.cvar(covindex,1).FaceColor] = deal(colourlist(covindex,:));
  [linehandle.cvar(covindex,2:end).Color, patchhandle.cvar(covindex,2:end).FaceColor] = deal(colourlist2(covindex,:));
end
[patchhandle.rate(:).FaceAlpha] = deal(.3);
for covindex = 1:size(linehandle.rate,1)
  [linehandle.rate(covindex,1).Color, patchhandle.rate(covindex,1).FaceColor] = deal(colourlist(covindex,:));
  [linehandle.rate(covindex,2:end).Color, patchhandle.rate(covindex,2:end).FaceColor] = deal(colourlist2(covindex,:));
end

for column = 1:3
  hd.noise.ax(6,column).XLim(1) = min(hd.noise.ax(6,column).Children(end).XData);
  hd.noise.ax(6,column).XLim(2) = max(hd.noise.ax(6,column).Children(end).XData);
end

ylist = 1:-.1:.1;
clist = [3,2,1,4,5,6];
for k = 1:6
  legendcolour(k,:) = hd.noise.ax(6,1).Children(7-k).Color;
  legendstring{k} = hd.noise.ax(6,1).Children(k).String;
  if k < 4
    legendstring{k} = [legendstring{k},' (wide)'];
  end
end
axes(hd.noise.ax(6,1))
hd.noise.ax(6,1).YLim = [0 6.5];
hold on
count = 0;
for k = [1, 4:6]
  count = count+1;
  t(count) = text(0,0,legendstring{k});
  t(count).Units = 'normalized';
  switch hd.noise.ax(6,1).XScale
    case {'lin','linear'}
      t(count).Position = [.3, ylist(count), 0];
    case {'log','logarithmic'}
      t(count).Position = [.2, ylist(count), 0];
  end
  t(count).Color = legendcolour(k,:);
  t(count).VerticalAlignment = 'middle';
  
  xlimit = hd.noise.ax(6,1).XLim;
  ylimit = hd.noise.ax(6,1).YLim;
  if ~isempty(zerolocation)
    xlimit(1) = zerolocation;
  end
  t(count).Units = 'Data';
  switch hd.noise.ax(6,1).XScale
    case {'lin','linear'}
      p(count) = plot(diff(xlimit)*[.1 .25], diff(ylimit)*ylist(count)*[1 1]);
    case {'log','logarithmic'}
      p(count) = plot([1.7*xlimit(1) .9*t(count).Position(1)], diff(ylimit)*ylist(count)*[1 1]);
  end
  p(count).LineWidth = 2;
  p(count).LineStyle = '-';
  p(count).Color = t(count).Color;
end
hold off



% remove classic labels next to each line
for column = 1:3
  switch column
    case 1
      [hd.noise.ax(6,column).Children(9:14).Visible] = deal('off');
    otherwise
      [hd.noise.ax(6,column).Children(1:6).Visible] = deal('off');
  end
end

% colour legend as desired
count = 0;
for k = [6, 4,2,1]
  count = count+1;
  [hd.noise.ax(6,1).Children(2*count-1).Color, ...
   hd.noise.ax(6,1).Children(2*count).Color] = deal(hd.noise.ax(6,1).Children(end-12+2*k).FaceColor);
end

[t1.Units, t2.Units, t3.Units] = deal('data');

t2.Position(2) = 3.9;
t3.Position(2) = 2.6;

% add markers for split x axis
if ~isempty(zerolocation)
  x1 = .045;
  dx = .02;
  x2 = .07;
  dy = .3;
  for column = 1:3
    hd.noise.ax(7,column) = axes;
    hd.noise.ax(7,column).Units = hd.noise.ax(6,column).Units;
    hd.noise.ax(7,column).Position([1 3]) = hd.noise.ax(6,column).Position([1 3]);
    hd.noise.ax(7,column).Position(4) = 30;
    hd.noise.ax(7,column).Position(2) = hd.noise.ax(6,column).Position(2) ...
                                        - hd.noise.ax(7,column).Position(4)/2;
    hd.noise.ax(7,column).FontSize = hd.noise.ax(6,column).FontSize;
    hd.noise.ax(7,column).Color = 'none';
    [hd.noise.ax(7,column).XAxis.Visible, ...
     hd.noise.ax(7,column).YAxis.Visible] = deal('off');
    hold on
    plot([x1+dx/2 x2+dx/2],[0 0],'-w','Linewidth',2)
    plot([x1,x1+dx],[-dy,dy],'-k')
    plot([x2,x2+dx],[-dy,dy],'-k')
    text(0,-2.32*dy,'0','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fontsize)
    hold off
    hd.noise.ax(7,column).XLim = [0 1];
    hd.noise.ax(7,column).YLim = [-1 1];
    hd.noise.ax(6,column).TickLength = [.04 0];
    hd.noise.ax(6,column).XTick = [.1 1];
  end
end


hd.noise.ax(3,1)

for letter = 1:3
  hd.noise.lt(letter).Position = [hd.noise.ax(6,letter).Position(1) - 74, ...
                                  hd.noise.ax(6,letter).Position(2) + panel.height + 2, 20, 20];
end

% close(rawhd.noise.fg)
hd.noise.fg.Name = 'Noise figure';

% swap bottom right and bottom left panels, but keep legend on bottom left
bottomleft = hd.noise.ax(6,1).Position;
bottomright = hd.noise.ax(6,3).Position;
hd.noise.ax(6,1).Position = bottomright;
hd.noise.ax(6,3).Position = bottomleft;
oldlegend = hd.noise.ax(6,1).Children(1:8);
newlegend = copyobj(oldlegend, hd.noise.ax(6,3));
[oldlegend.Visible] = deal('off');
for k = 1:2:numel(newlegend)
  newlegend(k).YData = oldlegend(k).YData / 6 - .55;
end
for k = 2:2:numel(newlegend)
  newlegend(k).Position(2) = oldlegend(k).Position(2) / 6 - .55;
end


% move wide box data to bottom layer
uistack(hd.noise.ax(6,1).Children(15:16),'bottom')

uistack(hd.noise.ax(6,3).Children(15:16),'bottom')



% % change colour of line plots as discussed for the paper
changecolour = true;
if changecolour
  
  lineobject  = findobj([hd.noise.ax(6,:)],'Type','Line');
  patchobject = findobj([hd.noise.ax(6,:)],'Type','Patch');
  textobject =  findobj([hd.noise.ax(6,:)],'Type','Text');

  linecolour = NaN(numel(lineobject),3);
  for k = 1:numel(lineobject)
    linecolour(k,:) = lineobject(k).Color;
  end
  patchcolour = NaN(numel(patchobject),3);
  for k = 1:numel(patchobject)
    patchcolour(k,:) = patchobject(k).FaceColor;
  end
  textcolour = NaN(numel(textobject),3);
  for k = 1:numel(textobject)
    textcolour(k,:) = textobject(k).Color;
  end
  uniquecolour = unique([linecolour;patchcolour;textcolour],'rows');

  linecolourtype  = NaN(size(lineobject));
  patchcolourtype = NaN(size(lineobject));
  textcolourtype  = NaN(size(lineobject));
  for c = 1:size(uniquecolour,1)
    for k = 1:numel(lineobject)
      if logical(prod(lineobject(k).Color == uniquecolour(c,:)))
        linecolourtype(k) = c;
      end
    end
    for k = 1:numel(patchobject)
      if logical(prod(patchobject(k).FaceColor == uniquecolour(c,:)))
        patchcolourtype(k) = c;
      end
    end
    for k = 1:numel(textobject)
      if logical(prod(textobject(k).Color == uniquecolour(c,:)))
        textcolourtype(k) = c;
      end
    end
  end

  greencolour = summer(10);
  greycolour = .65*[1 1 1];

  newcolour = repmat(greycolour,[10 1]);             
	newcolour([7 3 5],:) = greencolour([1 3 5],:); % replacing orange, green, blue

  for k = 1:numel(lineobject)
    lineobject(k).Color = newcolour(linecolourtype(k),:);
  end
  for k = 1:numel(patchobject)
    patchobject(k).FaceColor = newcolour(patchcolourtype(k),:);
  end
  for k = 1:numel(textobject)
    textobject(k).Color = newcolour(textcolourtype(k),:);
  end

end


% reset y axis limits for prettiness (previously shifted because of invisible labels)
hd.noise.ax(6,1).YLim = [0 4.2];
% hd.noise.ax(6,2).YLim = [4e2 2.13e4];
hd.noise.ax(6,2).YLim = [1e1 2.5e4];
% hd.noise.ax(6,3).YLim = [.8 1.002];
hd.noise.ax(6,3).YLim = [.75 1.1];
offset  = newlegend(4).Position(2);
spacing = newlegend(4).Position(2) - newlegend(4).Position(2);
axislim = hd.noise.ax(6,3).YLim;
relativespread = [.2 .5];
for k = 1:numel(newlegend)
  newposition(k) = axislim(1) + relativespread(1)*diff(axislim) ...
                              + relativespread(2)*diff(axislim)*2*floor((k-1)/2)/numel(newlegend);
end
for k = 1:2:numel(newlegend)
  newlegend(k).YData = newposition(k)*ones(size(newlegend(k).YData));
end
for k = 2:2:numel(newlegend)
  newlegend(k).Position(2) = newposition(k);
end


for k = 1:numel(partialfigurelist)
  close(k)
end
figure(hd.noise.fg)