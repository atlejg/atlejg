function helix(fname,height,n_rot,diam,rotation,qc)
%helix input: fname,height,n_rot,diam,rotation,qc

t = n_rot*2*pi*[0:0.01/n_rot:1];
z = height*t/(2*pi)/n_rot; % calculate this before rotating

t = t - rotation; 

x=(diam/2)*cos(t);
y=(diam/2)*sin(t);

if (qc) plot3(x,y,z);, end

fid = fopen(fname,'w');
if (fid < 0) 
   error(['can not open ',fname]);
else
   disp(sprintf('creating file "%s"',fname));
end

for (i=1 : length(t))
   fprintf(fid,'vertex create coordinates %e %e %e\n',x(i),y(i),z(i));
end

fclose(fid);

