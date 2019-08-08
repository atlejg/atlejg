function disp_profiles01(directory, pattern, decim, dt, disp_mode, xlab, ylab, zlab)

files     = dir(sprintf('%s/%s',directory,pattern));
nfiles = length(files);
disp(sprintf('nfiles= %i',nfiles));

for (i=1 : decim : nfiles)
   data = read_xyfile(sprintf('%s/%s',directory,files(i).name),0);

   % gotta do some initial operations when we know length of files 
   % (ie npoints)
   if (i==1)
      npoints = size(data,1);
      all = zeros(npoints,nfiles);
      l = zeros(1,npoints);
      dl = diff(data(:,1:3));
      for (j=2 : npoints)
         l(j) = l(j-1) + norm(dl(j-1));
      end
   end
   all(:,i) = data(:,4);
end

t = dt*[0:nfiles-1];
figure;
if (disp_mode == 0)
   contourf(t,l,all);
else
   surf(t,l,all);
   zlabel(zlab);
end
xlabel(xlab);
ylabel(ylab);
colorbar
