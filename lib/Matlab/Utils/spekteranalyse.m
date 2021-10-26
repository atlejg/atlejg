function spekteranalyse(fname,skip,plot_shifted,plot_semilogy,freq_limits,tscaler)
%input: fname,skip,plot_shifted,plot_semilogy,freq_limits,tscale
% main functionality done by Magne W. Mathiesen (as far as i know ... -atle)

[t,y] = textread(fname,'%n%n%*[^\n]','headerlines',skip);

t = t*tscaler;

subplot(211);

if (~plot_shifted) plot(t,y);,end
y = y - mean(y);
if (plot_shifted) plot(t,y);,end

N = 2*floor(length(y)/2);
sampf = 1/mean(diff(t));
Y = fft(y(1:end),N);
%Pyy = Y.* conj(Y)/N/2000;
Pyy = real(Y).*real(Y)/N/1000;
f = sampf*(0:N/2)/N;

lh = title(sprintf('%s (tscaler= %f)',fname,tscaler));
set(lh,'interpreter','none');  % make underscore be rendered properly (no tex-interpretation)
subplot(212);
if (plot_semilogy) 
   semilogy(f,Pyy(1:N/2+1))
else
   plot(f,Pyy(1:N/2+1))
end
title('Frequency content of y')
xlb = 'frequency (Hz)';
xlabel(xlb)
if (freq_limits == 0)
   axis([min(f) max(f) 0 max(Pyy)*1.05]);
   axis([0 max(f) 0 max(Pyy)*1.05]);
else
   axis(freq_limits);
end

saveas(gcf,[fname '.jpg']);
