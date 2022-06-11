% script with functions to reproduce results from Calaim*, Dehmelt*, Goncalves* & Machens 2021
% (panels from the main figure "Synaptic transmission delays cause uninformed spikes, 
% but networks with high-dimensional inputs are less affected.")


fgfull = figure(1);
axfull = fgfull.Children;
axfull.XLim = [2 axfull.XLim(2)];
axfull.Units = 'pixels';
fgfull.Position = [50 330 550 450];
axfull.Position = [70 60 400 330];
axes(axfull)
hold on
plot([2 50], [1 1], '--k')
plot([2 50], [0 0], '--k')
hold off
funkyscale = true;
axfull = scaleYAxis(axfull,funkyscale,fontsize);
axfull = formatAxis(axfull,fontsize);
axfull.Title.String = [axfull.Title.String, ' (default box, full connectivity)'];

fgwide = figure(2);
axwide = fgwide.Children;
axwide.XLim = axfull.XLim;
axwide.Units = 'pixels';
fgwide.Position = [650 330 550 450];
axwide.Position = axfull.Position;
axes(axwide)
hold on
plot([2 50], [1 1], '--k')
plot([2 50], [0 0], '--k')
hold off
funkyscale = true;
axwide = scaleYAxis(axwide,funkyscale,fontsize);
axwide = formatAxis(axwide,fontsize);
axwide.Title.String = [axwide.Title.String, ' (wide box, full connectivity)'];

fgred = figure(3);
axred = fgred.Children;
axred.Units = 'pixels';
fgred.Position = [1250 330 550 450];
axred.Position = axfull.Position;
axred.XLim = axfull.XLim;
axes(axred)
hold on
plot([2 50], [1 1], '--k')
plot([2 50], [0 0], '--k')
hold off
funkyscale = true;
axred = scaleYAxis(axred,funkyscale,fontsize);
axred = formatAxis(axred,fontsize);
axred.Title.String = [axred.Title.String, ' (default box, reduced connectivity)'];



function axishandle = scaleYAxis(axishandle,funkyscaleflag,fontsize)

  if funkyscaleflag
    
    funkyoffset = 0;
    funkyticklist = -3.0:.5:1;

    childlist = axishandle.Children;
    textlist  = findobj(childlist,'type','Text');
    linelist  = findobj(childlist,'type','Line');
    patchlist = findobj(childlist,'type','Patch');
    numpatch = numel(patchlist);
    otherlinelist = linelist(1:end-numpatch);
    datalinelist  = linelist(end-numpatch+1:end);

    for k = 1:numel(patchlist)
      patchlist(k).YData = exp(patchlist(k).YData + funkyoffset);
    end
    for k = 1:numel(datalinelist)
      datalinelist(k).YData = exp(datalinelist(k).YData + funkyoffset);
    end
    for k = 1:numel(otherlinelist)
      otherlinelist(k).YData = exp(otherlinelist(k).YData + funkyoffset);
    end
    for k = 1:numel(textlist)
      textlist(k).Position(2) = exp(textlist(k).Position(2) + funkyoffset);
    end

    axishandle.YTick = funkyticklist;
    axishandle.YTick = exp(axishandle.YTick + funkyoffset);
    axishandle.YTickLabel = {};
    for k = 1:numel(axishandle.YTick)
      if ismember(k, numel(axishandle.YTick)-[4 2 1 0])
        axishandle.YTickLabel{k} = num2str(funkyticklist(k));
      else
        axishandle.YTickLabel{k} = '';
      end
    end
    
  end
  
  for k = 1:numel(textlist)
    textlist(k).FontSize = fontsize;
  end
  
end



function axishandle = formatAxis(axishandle,fontsize)

  axishandle.FontSize = fontsize;
  axishandle.TitleFontSizeMultiplier = 1;
  axishandle.LabelFontSizeMultiplier = 1;

  axishandle.XLabel.String = 'dimensions';
  axishandle.YLabel.String = 'performance';
  axishandle.YLabel.Visible = 'on';  
  
  axishandle.Title.Units = 'normalized';
  axishandle.Title.Position = [.5 1.08 0];
  
  childlist = axishandle.Children;
  textlist  = findobj(childlist,'type','Text');
  
  % shift labels if the axis limits are not the default [2 50]:
  for k = 1:numel(textlist)
    textlist(k).Position(1) = textlist(k).Position(1) * axishandle.XLim(2)/50;
  end

end