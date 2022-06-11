function CDGM2021_ProcessFile(datafolder)

  filelist = dir([datafolder,'\*.mat']);
  isfile = ~[filelist.isdir];
  filelist = filelist(isfile);

  badfile = [];
  badfilestring = [];
  
  for file = 1:numel(filelist)

    showProgress(file,numel(filelist))    % subfunction included below

    badfilestring = cell(1);
    try
      data(file) = load([datafolder,filelist(file).name], 'par','performance','hon'); %#ok<AGROW>
    catch
      badfile = [badfile; file]; %#ok<AGROW>
      badfilestring{end+1,1} = filelist(file).name; %#ok<AGROW>
    end
    
  end
  
  if isempty(badfilestring{:})
    disp('No files contained unexpected data. That''s good.')
  else
    disp({'The following files contained unexpected data (empty : ',badfilestring{:}}') %#ok<CCAT>
  end

  if ~isempty(badfile)
    if badfile(end) == numel(filelist)
      data(badfile(1:end-1)) = [];
    else
      data(badfile) = [];
    end
  end

  tic
  rawdata = data;
  data = rawdata;
  for k = 1:numel(data)
    data(k).par.recurrent = [];
  end
  if exist([datafolder,'summary'],'dir') == 0
    mkdir([datafolder,'summary'])
  end
  save([datafolder,'summary\datasummary.mat'],'data','-v7.3')
  toc

end



function showProgress(progress,target)
  
  countstring = [num2str(progress),'/',num2str(target),'...'];
  
  if progress == 1
    fprintf(['\n','Processing .mat file ',countstring])
  else
    for numdigit = 0:log10(max(0,progress-1)) + numel(num2str(target)) + 4
      fprintf('\b'); % delete previous counter display
    end
    fprintf(countstring);
  end
  
  if progress == target
    fprintf(' done.\n')
  end
    
end