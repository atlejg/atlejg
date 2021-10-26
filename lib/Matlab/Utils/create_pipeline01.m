function create_pipeline01(journal,coord,r,verbosity)
%creates a Gambit journal file for straight connected 2D pipeline
% (ala OLGA).
% assumes coord is (x0,y0) ,(x1,y1) ... ,(xn,yn) and radius is r1,r2, ... ,rn
% NOTE: It adds a first horizontal section so that its easy (or possibel)
% to rotate the line into its proper position afterwards 
% (this could be done in the journal - i guess - but its not implemented yet ...
% atle.j.gyllensten@hydro.com
% Fri May  4 10:15:12 NST 2007

overlap = 4 * max(r);

% adding a horizontal first section so that it is easy to rotate
% the geometry later
tmp = zeros(length(coord) + 1,2);
tmp(1,:) = [coord(1,1) - 2*overlap, coord(1,2)];
tmp(2:end,:) = coord(:,1:2);
coord = tmp;

tmp = zeros(length(r) + 1,1);
tmp(1) = r(1);
tmp(2:end) = r(:,1);
r = tmp;

ps = zeros(length(r),2); % ps = pipesections
l  = zeros(length(r),1);
% finding pipe-sections and lengths
for (i = 1:length(ps))
   ps(i,:) = coord(i+1,:) - coord(i,:);
   l(i) = norm(ps(i,:));
end

% finding angles between pipe-sections
for (i = 1:length(ps)-1)
   a(i) = 0.5*real(acosd(dot(-ps(i,:),ps(i+1,:))/(l(i)*l(i+1))));
   % need cross product to decide whether we go 'left' or 'right' from here
   v1 = [ps(i,:) 0];
   v2 = [ps(i+1,:) 0];
   c = cross(v1,v2);
   if (c(3) >= 0)
      a(i) = 90 - a(i);
   else
      a(i) = a(i) - 90;
   end
end

pid = fopen(journal,'w');

if (pid < 0 ) error 'cannot open file', end

fprintf(pid,'reset\n');
fprintf(pid,'face create "for_split" width 1 height 1 yzplane rectangle\n');

% first two pieces must be handled separetely
fprintf(pid,'volume create height %.3f radius1 %.3f radius3 %.3f offset %.3f 0 0 xaxis frustum\n',l(1)+overlap,r(1),r(1),(l(1)+overlap)/2);
fprintf(pid,'/\n');
fprintf(pid,'volume move "volume.1" offset %.3f 0 0 connected\n',-l(1));
fprintf(pid,'/\n');

fprintf(pid,'volume create height %.3f radius1 %.3f radius3 %.3f offset %.3f 0 0 xaxis frustum\n',l(2)+2*overlap,r(2),r(2),(l(2)+2*overlap)/2);
fprintf(pid,'volume move "volume.2" offset %.3f 0 0 connected\n',-overlap);
fprintf(pid,'volume move "volume.2" dangle %.3f vector 0 0 1 origin 0 0 0 connected\n',-a(1));
fprintf(pid,'volume move "volume.1" dangle %.3f vector 0 0 1 origin 0 0 0 connected\n',a(1));
fprintf(pid,'volume split "volume.1" faces "for_split" connected keeptool\n');
fprintf(pid,'volume split "volume.2" faces "for_split" connected keeptool\n');
fprintf(pid,'volume delete "volume.2" "volume.3" lowertopology\n');
fprintf(pid,'volume split "volume.4" volumes "volume.1" connected bientity\n');
fprintf(pid,'volume move "volume.4" dangle %.3f vector 0 0 1 origin 0 0 0 connected\n',a(1));
fprintf(pid,'volume move "volume.4" offset %.3f 0 0 connected\n',-l(2));
fprintf(pid,'/\n');

% then we can generalize
for (i = 3 : length(ps))
   fprintf(pid,'volume create height %.3f radius1 %.3f radius3 %.3f offset %.3f 0 0 xaxis frustum\n',l(i)+2*overlap,r(i),r(i),(l(i)+2*overlap)/2);
   fprintf(pid,'volume move "volume.%d" offset %.3f 0 0 connected\n',4*(i-2)+2,-overlap);
   fprintf(pid,'volume move "volume.%d" dangle %.3f vector 0 0 1 origin 0 0 0 connected\n',4*(i-2)+2,-a(i-1));
   fprintf(pid,'volume move "volume.%d" dangle %.3f vector 0 0 1 origin 0 0 0 connected\n',4*(i-2),a(i-1));
   fprintf(pid,'volume split "volume.%d" faces "for_split" connected keeptool\n',4*(i-2));
   fprintf(pid,'volume split "volume.%d" faces "for_split" connected keeptool\n',4*(i-2)+2);
   fprintf(pid,'volume delete "volume.%d" "volume.%d" lowertopology\n',4*(i-2)+3,4*(i-2)+2);
   fprintf(pid,'volume split "volume.%d" volumes "volume.%d" connected bientity\n',4*(i-2)+4,4*(i-2));
   fprintf(pid,'volume move "volume.%d" dangle %.3f vector 0 0 1 origin 0 0 0 connected\n',4*(i-2)+4,a(i-1));
   fprintf(pid,'volume move "volume.%d" offset %.3f 0 0 connected\n',4*(i-2)+4,-l(i));
   fprintf(pid,'/i= %i\n',i);
end

fclose(pid);

% reporting
if (verbosity)
   l
   a
   plot(coord(:,1),coord(:,2),'-*');
end


