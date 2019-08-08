function eclipse_contourplot(directory, pattern, varnm, n_tsteps, well_length, well_start)
%plots a variable in time and space as a contour

files = dir(sprintf('%s/%s',directory,pattern));
nfiles = size(files,1);

data = [];

for (fileno = 1 : nfiles)
   [t, tmp] = read_frft_file_02(sprintf('%s/%s',directory,files(fileno).name),varnm,n_tsteps);
   % skip init data
   data(fileno,:,:) = tmp(1:end,:);
end

md = well_length/size(data,3) * [0 : size(data,3)-1] + well_start;

% calculate tot volumes
accum = zeros(nfiles, length(t)-1);
for fileno=1 : nfiles
   d(:,:) = data(fileno,:,:);
   tot = sum(d');
   for i=2 : length(t)
      accum(fileno,i-1) = trapz(t(1:i), tot(1:i));
   end
end

% max/min values
max_accum = max(max(accum));
max_rate = max(max(max(data)));
min_rate = min(min(min(data)));

d = [];
for (fileno = 1 : nfiles )
   d(:,:) = data(fileno,:,:);

   figure

   subplot('position',[0.1 0.1 0.85 0.65]);
   contourf(t,md,d');
   
   caxis([min_rate max_rate]);
   colorbar;

   xlabel('time [days]');
   ylabel('md [m]');
   zlabel(sprintf('relative %s [-]',varnm));

   subplot('position',[0.1 0.80 0.73 0.1]);
   plot(t(2:end),accum(fileno,:),'r-','LineWidth',3)
   axis([t(1) t(end) 0 max_accum]);
   set(gca,'YMinorGrid','on');
   set(gca,'xticklabel',[])
   ylabel(sprintf('accum.\n%s',varnm(4:end)));
   titl = files(fileno).name(1:3);
   title(titl, 'interpreter', 'none'); % avoid tex-interpretation
   saveas(gcf, sprintf('%s/%s_%s.jpg', directory, titl, varnm))
   disp(sprintf('file: %s max %s = %e', titl, varnm, max(max(d))));
end
