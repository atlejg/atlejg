function data = lesdata2(filnavn,skip)
%  lesdata:   Leser ASCII datafiler.  Som LOAD men tolererer tekstlinjer og feil.
%  Fila ma inneholde kolonner med tall.  
%  Første linje med tall ma ha rett antall kolonner.  
%  Linjer hvor antall tolkbare kolonner avviker fra dette tolkes som feil.
%  Tall etter (til hoyre for) en tekstkolonne tolkes ikke.
%  BRUK: data = lesdata2('filanvn',skip);
%  
%  Se ogso: LOAD, DLMREAD, TEXTREAD(Matlab 5.3+)

% Are Mjaavatten, april 1999

% lesdata2: skip header lines. 

fid = fopen(filnavn);
if fid < 0
   error(['Kunne ikke opne ',filnavn]);
end
j = 0;
data = [];
while 1
   j = j+1;

   line = fgetl(fid);

   if (j <= skip), continue, end

   if ~isstr(line), break, end
   try
      data = [data;sscanf(line,'%g')'];
   catch
      disp(['Kunne ikke lese linje ',num2str(j),' i ',filnavn]);
   end
end
fclose(fid);
