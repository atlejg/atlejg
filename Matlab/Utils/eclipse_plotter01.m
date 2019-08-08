% plots eclipse timeseries from RSM-files.
% see 'vars' for which variables that are plotted (against TIME), and 'files' for which files 
% are used.
% havent made it configurable (yet), so make a copy of this one and modify it!
% 2008-08-08
% atle j. gyllensten


files = {'simu_01/QFIVE.RSM', 'simu_02/QFIVE.RSM', 'simu_03/QFIVE.RSM','simu_04/QFIVE.RSM'};
casenames = files; % in some situations, we could have casenames different from filenames

n_headerlines = 6; % will be skipped

% could be changed ...

% variables = columns in RSM-file
TIME     = 1;
YEARS    = 2;
BPR      = 3;
FOPR     = 4;
WBHP_I   = 5;
WBHP_P   = 6;
FWCT     = 7;
FPR      = 8;

% could be changed
vars = [FOPR WBHP_I FWCT WBHP_P FPR BPR];

% SOMETIMES ITS USEFUL TO HAVE NAME OF VARIABLES
VARNAMES{TIME}     = 'TIME';    
VARNAMES{YEARS}    = 'YEARS';    
VARNAMES{BPR}      = 'BPR';    
VARNAMES{FOPR}     = 'FOPR';    
VARNAMES{WBHP_I}   = 'WBHP_I';    
VARNAMES{WBHP_P}   = 'WBHP_P';    
VARNAMES{FWCT}     = 'FWCT';    
VARNAMES{FPR}      = 'FPR';    

% ---- END OF INPUT ----

% technicality...
files = char(files);
casenames = char(casenames);

figlist = [];

nfiles = size(files,1);
legs = {};  % legend are tricky, as always

% read files one by one and plot variables for each file. plots for same variables are
% collected in one figure
for (fileno = 1 : nfiles)
   file = files(fileno,:)
   casenm = casenames(fileno,:)
   d=lesdata2(file,n_headerlines);

   varno  = 1;
   for (var = vars)
      if (fileno == 1)
         figlist(varno) = figure;
         lh = title(VARNAMES(var));
         set(lh,'interpreter','none');  % make underscore be rendered properly (no tex-interpretation)
         hold on;
         legs{varno} = [];
      end

      figure(figlist(varno));
      plot(d(:,TIME),d(:,var),sprintf('%s*',colour(fileno))); % using colour.m to have different marking for each curve
      legs{varno} = [legs{varno};struct('txt',casenm)];

      varno = varno + 1;
   end
end

% attach legends
for (varno = 1 : length(vars))
   figure(figlist(varno));
   lh = legend(legs{varno}.txt);
   set(lh,'interpreter','none');  % make underscore be rendered properly (no tex-interpretation)
end
