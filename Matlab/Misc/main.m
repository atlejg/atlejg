function main
% reads and displays xy-plots (from fluent)

% get list of xy-files
xyFiles = dir('*.xy');

% all plots in one figure
hold on;

% read and plot files - one by one
for (fileNo = 1 : size(xyFiles, 1) )
   xy = lesdata(xyFiles(fileNo).name);
   plot(xy(:,1), xy(:,2),color_(fileNo) );   
end

% use filenames as legend
legend(xyFiles.name);


% ========= SUBS =========
function color = color_(no)
% want different curves to appear different...

colors = 'bgrcmyk';

% pick a colour
color = colors(mod(no,size(colors, 2) ) + 1);
