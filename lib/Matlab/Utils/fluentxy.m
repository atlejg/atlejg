function fluentxy(x ,y, fileNm, title, label1, label2)
%fluentxy. Creates a xy file for fluent based on x and y vectors

error(nargchk(2, 6, nargin) );

if (nargin < 3) fileNm = 'fromMatlab.xy'; end
if (nargin < 4) title  = 'notitle'; end
if (nargin < 5) label1 = 'Position'; end
if (nargin < 6) label2 = 'none'; end

fid = fopen(fileNm, 'w');

if (fid < 0)
   error('cannot open file');
end

n=size(x);

fprintf(fid, '(title "%s")\n', title);
fprintf(fid, '(labels "%s" "%s")\n', label1, label2);
fprintf(fid, '\n');
fprintf(fid, '((xy/key/label "%s")\n', title);

for (k=1 : n)
   fprintf(fid,'%1.6f\t%1.6f\n',x(k),y(k));
end

fprintf(fid, ')\n');

fclose(fid);
