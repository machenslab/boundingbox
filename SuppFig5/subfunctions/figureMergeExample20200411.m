function figureMergeExample(varargin)
  % arguments are an arbitrary number of figure indices to be merged, e.g. (1,2,3,4)
  
  numrawfg = numel(varargin);
  for k = 1:numrawfg
    hd.raw(k).fg = figure(varargin{k});
  end
  numrawax = numel(hd.raw(1).fg.Children);
  
  hd.merged.fg = figure();
  hd.merged.fg.Color = [1 1 1];
  hd.merged.fg.Position = [50 50 hd.raw(1).fg.Position(3) numrawfg/2*hd.raw(1).fg.Position(4) + 20];
  
  for k = 1:numrawfg
    hd.merged.ax(k,1:numrawax/2) = copyobj(hd.raw(k).fg.Children(1:floor(end/2)), hd.merged.fg)';
    for a = 1:numrawax/2
      hd.merged.ax(k,a).Position(2) = hd.merged.ax(k,a).Position(2) + ...
                                      (numrawfg-k)/2*hd.raw(1).fg.Position(4);
    end
  end
  
  % remove superfluos axis ticks and axis labels
  for k = 1:numrawfg-1
    [hd.merged.ax(k,3:4).XTickLabel] = deal([]);
    for a = 1:numrawax/2
      hd.merged.ax(k,a).XLabel.Visible = 'off';
    end
  end
  
  % add column titles at the top of the figure
  for a = 1:numrawax/2
    hd.merged.ax(1,a).Title = hd.raw(1).fg.Children(a).Title;
    hd.merged.ax(1,a).Title.Visible = 'on';
    hd.merged.ax(1,a).Title.Units = 'pixels';
    hd.merged.ax(1,a).Title.Position(2) = hd.merged.ax(1,a).Title.Position(2) + 5;
  end
  
  % add panel lettering
  fontsize = 11;
  for letter = 1:2*numrawfg
    hd.merged.lt(letter) = axes();
    hd.merged.lt(letter).Units = 'pixels';
    axis off
    text(1,1,char(-32+96+letter),'FontWeight','normal','FontSize',fontsize*1.2)
    row = ceil(letter/2);
    column = mod(letter,2)+3;
    reference = hd.merged.ax(row,column).Position;
    hd.merged.lt(letter).Position = [reference(1) - 60, reference(2) + reference(4) - 7, 20, 20];
  end
  
  % add new row label on the left
  for k = 1:numrawfg
    networksize = floor(hd.merged.ax(k,numrawax/2-1).YLim(2));
    dimension = numel(hd.merged.ax(k,numrawax/2).Children);
    redundancy = networksize/dimension;
%     hd.merged.ax(k,numrawax/2).YLabel.String = ['\rho = ',num2str(redundancy)];
%     if k>1 && strcmp(hd.merged.ax(k,numrawax/2).YLabel.String, ...
%                      hd.merged.ax(k-1,numrawax/2).YLabel.String)
%       hd.merged.ax(k,numrawax/2).YLabel.String = [hd.merged.ax(k,numrawax/2).YLabel.String,' (wide)'];
%     end
    switch hd.merged.ax(k,numrawax/2).YLabel.String
      case 'full'
        colour = hd.merged.ax(k,numrawax/2).Children(end).Color;
        if colour(2) > colour(3) % more green than blue
%           hd.merged.ax(k,numrawax/2).YLabel.String = 'undelayed (narrow)';
          hd.merged.ax(k,numrawax/2).YLabel.String = 'default undelayed';
        else
%           hd.merged.ax(k,numrawax/2).YLabel.String = 'full connectivity';
          hd.merged.ax(k,numrawax/2).YLabel.String = 'default delayed';
        end
      case 'wide'
        hd.merged.ax(k,numrawax/2).YLabel.String = 'wider box';
      case 'reduced'
        hd.merged.ax(k,numrawax/2).YLabel.String = 'less excitation';
      case 'none'
        hd.merged.ax(k,numrawax/2).YLabel.String = 'inhibition only';
    end
    hd.merged.ax(k,numrawax/2).YLabel.FontSize = fontsize;
    hd.merged.ax(k,numrawax/2).YLabel.FontWeight = 'bold';
    hd.merged.ax(k,numrawax/2).YLabel.Units = 'pixels';
    hd.merged.ax(k,numrawax/2).YLabel.Position(1) = hd.merged.ax(k,numrawax/2).YLabel.Position(1) - 3;
    hd.merged.ax(k,numrawax/2).YLabel.Visible = 'on';
  end
  xshift = 15;
  hd.merged.fg.Position(3) = hd.merged.fg.Position(3) + xshift;
  kids = hd.merged.fg.Children;
  for k = 1:numel(kids)
    kids(k).Position(1) = kids(k).Position(1) + xshift;
  end
  
end