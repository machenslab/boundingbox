function updateMergeExample(figureindex)
  
  fg = figure(figureindex);
  ax = fg.Children((end/3+1+3):4:end);
  
  ymin = [];
  ymax = [];
  
  for k = 1:numel(ax)
    
    ymin = min([ymin, ax(k).YLim(1)]);
    ymax = max([ymax, ax(k).YLim(2)]);
    
  end
  
  for k = 1:numel(ax)
    
    ax(k).YLim = [ymin ymax];
    
  end
  
end