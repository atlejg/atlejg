function s = symbols02(sym_indx,col_indx,ln_indx)
%input 
% sym_indx: symbol
% col_indx: colour
% ln_indx : line

syms    = 'h*.<>+oxsd^vp';
colours = 'mkrcgyb';
lines   = '-:';

s = '';

if (sym_indx > 0)
   i = mod(sym_indx,length(syms))+1;
   s = [s syms(i)];
end

if (col_indx > 0)
   i = mod(col_indx,length(colours))+1;
   s = [s colours(i)];
end

if (ln_indx > 0)
   i = mod(ln_indx,length(lines))+1;
   s = [s lines(i)];
end
