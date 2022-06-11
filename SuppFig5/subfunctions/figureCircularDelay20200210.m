function figureCircularDelay(figureindex1,figureindex2)

  hd.raw.fg(1) = figure(figureindex1);
  hd.raw.fg(2) = figure(figureindex2);
  hd.composite.fg = figure();
  hd.composite.fg.Position = [50 50 760 390];
  hd.composite.fg.Color = [1 1 1];

  for source = 1:2
    
    offset = numel(hd.raw.fg(source).Children) - 12;

    hd.composite.ax(1+4*(source-1)) = copyobj(hd.raw.fg(source).Children(8+offset), hd.composite.fg);
    hd.composite.ax(2+4*(source-1)) = copyobj(hd.raw.fg(source).Children(11+offset), hd.composite.fg);
    hd.composite.ax(3+4*(source-1)) = copyobj(hd.raw.fg(source).Children(10+offset), hd.composite.fg);
    hd.composite.ax(4+4*(source-1)) = copyobj(hd.raw.fg(source).Children(9+offset), hd.composite.fg);

  end

  [hd.composite.ax.Units] = deal('pixels');

  % y1 = 260;
  % y2 = 40;
  y1 = 218;
  y2 = 40;

  x1 = 30;
  w1 = 220;
  h1 = 150;
  w2 = 70;
  d1 = 40;
  hd.composite.ax(1).Position = [x1 y1 w1 h1];
  hd.composite.ax(2).Position = [x1+w1+1.7*d1 y1 w1 h1];
  hd.composite.ax(3).Position = [x1+2*w1+2.4*d1 y1 w2 h1];
  hd.composite.ax(4).Position = [x1+2*w1+w2+3.2*d1 y1 w2 h1];

  for k = 5:8
    hd.composite.ax(k).Position = [hd.composite.ax(k-4).Position(1), y2, ...
                                   hd.composite.ax(k-4).Position(3:4)];
  end

  for k = 1:4, hd.composite.ax(k).Title.Units = 'pixels'; end
  for k = 1:4, hd.composite.ax(k).Title.Position(2) = 157; end
  for k = 5:8, hd.composite.ax(k).Title.Visible = 'off'; end
  for k = 1:4, hd.composite.ax(k).XLabel.Visible = 'off'; end
  for k = 1:4, hd.composite.ax(k).XTickLabel = []; end
  for k = [2 6], hd.composite.ax(k).YLabel.String = 'neuron ID'; end
  for k = [3 4 7 8], hd.composite.ax(k).YTickLabel = []; end
  
%   for k = [1 5], hd.composite.ax(k).YLim = 4*[-1 1]; end
  ymax = max(abs([hd.composite.ax(1).Children(2:end).YData, ...
                  hd.composite.ax(5).Children(2:end).YData]));
  for k = [1 5], hd.composite.ax(k).YLim = ymax*[-1 1]; end
  
  % newer style, without honeymoon highlight and without most axes and grids
  for k = [1 5]
    hd.composite.ax(k).Children(1).Visible = 'off';
  end
  for k = [2 6]
    hd.composite.ax(k).Children(end).Visible = 'off';
  end
  for k = 1:4
    hd.composite.ax(k).XAxis.Visible = 'off';
  end
  for k = [3 4 7 8]
    hd.composite.ax(k).XGrid = 'off';
  end
  for k = 1:8
    hd.composite.ax(k).Box = 'off';
  end
  
  
end