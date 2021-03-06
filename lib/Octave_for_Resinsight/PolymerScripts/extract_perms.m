grid_no = 0;                            % need a loop if using LGR's

permx = riGetGridProperty(grid_no, "PERMX");
permx(isnan(permx)) = -1.;               % replace NaN

maxp = max(permx(:,:,55:75),[], 3);   % max along z-axis around well
prof = mean(maxp(90:125,47:55), 2, 'h');   % harmonic average from producer to injector

% copy to all layers for convinience
maxperm = zeros(size(permx));
nz = size(permx)(3);
for k = [1:nz]
   maxperm(:,:,k) = maxp;
end
riSetGridProperty(maxperm/1000., grid_no, "MAX_PERMX");  % mD -> D

% and write profile to file
pid = fopen('perm_profile.txt', 'w');
ny = size(prof)(1)
for j = ny:-1:1
   fprintf(pid, '     PERMX      %.1f    1 20  %i %i  11 11 /\n', prof(j), ny-j+2, ny-j+2);
end
fclose(pid);

prof
%maxp(1,:)
%permx(90,:,55)
