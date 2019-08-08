function descr = describe01(vals,formatStr,prefix,unit)

descr = cell(size(vals));

for (k = 1 : length(vals))
   s = sprintf('%%s= %s%%s',formatStr);
   descr(k) = {sprintf(s,prefix,vals(k),unit)};
end
