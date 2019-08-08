d=0.2
h=1

t1 = 2*pi*[0:0.01:1];
t2 = 2*pi*[-.5:0.01:0.5];

z = h*t1/(2*pi);
y1=(d/2)*sin(t1);
x1=(d/2)*cos(t1);
y2=(d/2)*sin(t2);
x2=(d/2)*cos(t2);

plot3(x1,y1,z);
hold on
plot3(x2,y2,z);

fname1 = 'spiral1.jou'
fid1 = fopen(fname1,'w');
fname2 = 'spiral2.jou'
fid2 = fopen(fname2,'w');
if (fid1 < 0) error(['can not open ',fname1]), end
if (fid2 < 0) error(['can not open ',fname2]), end

for (i=1 : length(t1))
   fprintf(fid1,'vertex create coordinates %e %e %e\n',x1(i),y1(i),z(i));
   fprintf(fid2,'vertex create coordinates %e %e %e\n',x2(i),y2(i),z(i));
end

fclose(fid1);
fclose(fid2);

