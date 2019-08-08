function plot_xyfiles(directory,pattern,decimate,skip,debug_)
%useful for checking convergence of solution based on xy-files monitoring 
% input: directory, pattern, decimate, skip,debug_

% history
%
% Mon Sep  1 11:14:06 NST 2008


if (decimate ==0) error('decimate should be > 0'); , end

files = dir(sprintf('%s/%s',directory,pattern));

n = length(files);
files = files(1:decimate:n);
n = length(files);

if (n ==0) error('no files found'); , end

leg = [];
figure;
hold on;
for (i=1 : n)
   monit = read_xyfile(sprintf('%s/%s',directory,files(i).name),skip,debug_);
   % make sure time is monotonically increasing
   monit = sortrows(monit,1);
   plot(monit(:,1),monit(:,2),symbols02(5,i,2));
end
hold off;

lh = legend(files.name);
set(lh,'interpreter','none');
title(directory,'interpreter','none');

