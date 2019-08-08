function xy = read_xyfile(fname,skip,debug_)
%based on lesdata by Are Mjaavatten, april 1999
% not effective, should use textread

if (debug_)
   disp(sprintf('read_xyfile: fname= %s',fname));
end

fid = fopen(fname);
if fid < 0
   error(['Kunne ikke opne ',fname]);
end

% skip 4 header lines + user request
for (i=1 : 4 + skip)
   fgetl(fid);
end

j = 0;
xy = [];
while 1
   j = j+1;
   line = fgetl(fid);
   if ~isstr(line), break, end
   if (line == ')') break, end
   try
      xy = [xy;sscanf(line,'%g')'];
   catch
      disp(['Kunne ikke lese linje ',num2str(j),' i ',fname]);
   end
end
fclose(fid);

% make sure time is monotonically increasing
xy = sortrows(xy,1);
