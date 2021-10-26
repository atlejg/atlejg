function r = configread(filename)
% function configread(filename)
% Parses a config file (ini file) and stores the result in a struct (f)
% Uses # and % as comment style
% Author: Geir Arne Evjen
% Date:   15/06/2005
% Modifications:
%   Removed dependency of eval function (GAE:22/06/2005)
  
  fid   = fopen(filename);
  lines = [];
  while 1
    line    = fgetl(fid);
    if ~ischar(line), break; end
    [s,f,t] = regexp(line, '^(.*)(\%|\#).*$');    % Find comments
    if ~isempty(s) 
      if any(t{1}(1,:))
        sidx = t{1}(1,:); eidx = t{1}(2,:);
        line    = line(sidx(1):eidx(1)-1);
      else, continue; end
    end
    
    
    sline = strip(line);   % kill spaces from line
    if ~isempty(sline),
      lines{end+1} = line;   
    end
  end  
  
  fclose(fid);
 
  % Go over the values in lines and make cells 
  
  msec      = [];
  subsec    = []; subsecval = [];
  for i=1:length(lines)
    % i want it to be possible to give values as '[2 4 7]' (which is easy to eval)
    % so must make sure they are not interpreted as sections. -atle,
    %[s,f,t]     = regexp(lines{i}, '\[(.*)\]'); 
    [s,f,t]     = regexp(lines{i}, '^\s*\[(.*)\]'); 
    if ~isempty(s),
      msec{end+1}      = lines{i}(t{1}(1):t{1}(2));            
      subsec{end+1}    = {};
      subsecval{end+1} = {};
      continue;
    end
    [s, f, t] = regexp(lines{i}, '^(.*?)\s*=\s*(.*)$');
    if ~isempty(s) && length(t{1}(:,1)) == 2,
      nidx = t{1}(1,:); vidx = t{1}(2,:);
      subsec{end}{end+1}    = lines{i}(nidx(1):nidx(2));
      subsecval{end}{end+1} =  lines{i}(vidx(1):vidx(2));
      continue;  
    end
  end
  
  % Build structure
  lms = length(msec);
  r   = cell2struct(cell(1,lms), msec, 2);
  for i=1:lms
    if ~isempty(subsec{i})  
      r.(msec{i}) = cell2struct(subsecval{i},subsec{i},2);
    end  
  end
  
  
  
  
  
  
  
  
  
