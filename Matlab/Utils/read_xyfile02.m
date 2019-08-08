function xy = read_xyfile02(fname,skip,debug_)
%effective, uses textread

if (debug_)
   disp(sprintf('read_xyfile: fname= %s',fname));
end

% skip 4 header lines + user request
xy = textread(fname,'',-1,'headerlines',4 + skip);

% make sure time is monotonically increasing
xy = sortrows(xy,1);
