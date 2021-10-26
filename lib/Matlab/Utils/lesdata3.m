function data = lesdata3(fname,skip,max_n_data,ncols)
%lesdata3
% -> fname  : filename
% -> skip   : how many lines to skip (header)
% -> max_n_data   : max data lines to read (<= 0 means all)
% -> ncols  : number of columns to read

fid = fopen(fname);
if (fid < 0) error(['lesdata3 : can not open ',fname]), end

% build format string
fs = '%g';
for (i=1 : ncols-1)
   fs = sprintf('%s %%g',fs);
end

tmp  = zeros(ncols,1);
data = [];
head_count = 0;
data_count = 0;
while (1)
   head_count = head_count + 1;
   line = fgetl(fid);

   % skip header
   if (head_count <= skip) continue, end

   % have passed header. count data-lines
   data_count = data_count + 1;
   if ~isstr(line), break, end
   tmp = sscanf(line,fs,ncols);
   data(data_count,:) = tmp(1:ncols);

   if (max_n_data > 0 && data_count >= max_n_data) break, end
end

fclose(fid);
