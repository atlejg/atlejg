fig=figure;
set(fig,'Position',[100 100 150 600])
%directoryname = uigetdir('F:\CSSD\CSSD_v5', 'Pick a Directory');
directoryname = '.'
files=dir(strcat(directoryname,'/2_vof_0*0.tif'));
set(fig,'DoubleBuffer','on');
set(gca,'Position',[0 0 1 1],...
    'NextPlot','replace','Visible','off')
mov = avifile(strcat(directoryname,'/','2.avi'), 'FPS', 3);
for i=1:length(files)
    A=imread(strcat(directoryname,'/',files(i).name),'tif');
    image(A);
    F = getframe(gca);
    mov = addframe(mov,F);
end
mov = close(mov);
