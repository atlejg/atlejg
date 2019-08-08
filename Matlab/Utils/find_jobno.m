function selected = find_jobno(jobs,expr,indexes)
% FIND_JOBNO. Finds selected Flacs job numbers using regular expressions.
% -> jobs: vector of job numbers
% -> expr: regular expression to match
% -> indexes: just return indexes. Boolean. Optional
% <-     : selected jobs in a vector
% Example: find_jobno([010100:010199],'0101\d9') will return all jobs
%          starting with 0101 and ending with a 9 since '\d' 
%          represents any digit.
%          find_jobno([010100:010199],'0101\d9',1) does the same, but
%          will give you the indexes in stead.

error(nargchk(2,3,nargin) );
if (nargin == 2)
   indexes = 0;
end

% want them as a (row) vector cause they're eaxeier to loop
if (size(jobs,1) > 1) jobs = jobs'; end
   
selected = [];
index = 0;
for (job = jobs)
   index = index + 1;

   % wanna treat'em as strings
   job_str = num2str(job);
   
   if (length(job_str) == 5)
      % typically, the leading zeros are removed...
      job_str = strcat('0',job_str);
   end
   
   % check if this job is one we're looking for
   if (regexp(job_str,expr) )
      % yes it is, so we keep it, according to what user wants
      if (indexes)
         selected = [selected, index];
      else
         selected = [selected, job]; % we'll return numbers...
      end
   end
end

