function s1 = strip(s)
% function s1 = strip(s)
% Removes spaces from the front and the back of a string  
% Modification of Matlabs deblank function
%
% Author: Geir Arne Evjen
% Date:   14/6/05
  if isempty(s)
    s1 = s([]);
  else
    if ~isstr(s),
      warning('MATLAB:deblank:NonStringInput','Input must be a string.')
    end
    
    [r,c] = find( ~isspace(s) );
 
    if isempty(c),
      s1 = s([]);
   else
      s1 = s(:,min(c):max(c));
    end
  end  
  