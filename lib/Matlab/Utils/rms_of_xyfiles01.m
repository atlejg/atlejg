function rms = rms_of_xyfiles01(directory,pattern,decimate,skip,debug_)
%useful for checking convergence of solution based on xy-files
% input: directory,pattern,decimate,debug_

files = dir(sprintf('%s/%s',directory,pattern));

files = files(1:decimate:end);

n = length(files);
rms = zeros(n-1,1);

if (debug_) n , end

if (n <= 1)
   warning 'I NEED AT LEAST 2 FILES';
   return
end

for (i=1 : n-1)
   if (i == 1)  xy1 = read_xyfile(sprintf('%s/%s',directory,files(1).name),skip,debug_);,end
   xy2 = read_xyfile(sprintf('%s/%s',directory,files(i+1).name),skip,debug_);
   if (debug_) length(xy2);, end

   rms(i) = norm(xy1-xy2)/length(xy1);

   xy1 = xy2;
end
