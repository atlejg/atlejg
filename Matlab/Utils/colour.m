function col= colour(no)
% returns a color-character (to be used in for example 'plot').
% based on the number given.

% typical usage:

% directory = '../Fluent6.2.16/im11/TMP'
% pattern   = 'R*.sort'
% files = dir(sprintf('%s/%s',directory,pattern));
% clf
% hold on
% for (no = 1 : size(files, 1) )
%    data = lesdata(sprintf('%s/%s',directory,files(no).name));
%    plot(data(:,1),data(:,2),colour(no))
% end
% lh = legend(files.name);
% set(lh,'interpreter','none');  % make underscore be rendered properly (no tex-interpretation)

colours = 'bgrcmyk';
%marks   = '.ox+*sd^<>ph'; - dropping marks : overkill...
lines   = {'-' '--' '-.' ':'};

count = 0;
for (n = 1 : length(lines) )
   %for (m = 1 : length(marks) )
      for (k = 1 : length(colours) )
         count = count + 1;
         %combined{count} = strcat(colours(k),marks(m),lines{n});
         combined{count} = strcat(colours(k),lines{n});
      end
   %end
end

%combined

%colours = {'b' 'g' 'r' 'c' 'm' 'y' 'k' 'b--' 'g--' 'r--' 'c--' 'm--' 'y--' 'k--'};

% pick a colour - making sure we dont exceed array size.
%col = colours{mod(no,size(colours, 2) ) + 1};
col = combined{mod(no,size(combined, 2) ) + 1};
