function [t, y] = read_frft_file_02(file, varnm, n_tsteps)

nvals_pr_column = 4;

fid = fopen(file);
if fid < 0 error(['could not open file: ',file]);, end

% init
y = [];
t = [];
tstep  = 0;
lineno = 0;
nvals  = 0;
j      = 0; % for debug
counter = 0; % tricky...
while 1
   j = j + 1;
   line = fgetl(fid);
   if ~isstr(line) break, end
   
   % REQUESTED PROFILE
   if strfind(line,varnm)
       % append data for previous timestep - if any
      if tstep > 0 y(tstep,:) = vals; , end

      % reset
      nvals = str2double(line(22:23));  % done redunantly
      vals = [];
      lineno = 1;
      tstep = tstep + 1;
      disp(sprintf('file: %s - tstep = %i', file, tstep));
      if (tstep > n_tsteps) break , end
      continue;
   end
   
   % TIME
   if strfind(line,'TIME') , counter = 0; , end
   if counter == 1 , t = [t ;sscanf(line,'%g')]; , end
   counter = counter + 1;
       
   % are we done for this timestep?
   if lineno == 0 || (lineno-1) * nvals_pr_column > nvals , continue; , end

   vals = [vals;sscanf(line,'%g')];
   lineno = lineno + 1;
end

% dont forget last timestep (unless stopped by n_tsteps)
if vals , y(tstep,:) = vals; , end

fclose(fid);

