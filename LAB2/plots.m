filocaldo = importfile(file_name, [1, Inf]);
data = filocaldo.VarName2; % vettore dalla II colonna (letture HWA)
  
data = ((3.65+(5-data).^2)/7.83).^2; % curva di taratura corretta HWA
data = data - mean(data); % pulizia del segnale dal valor medio

figure(1);
time = 0:ts:ts*(length(data)-1);
plot(freq,abs(dft),'k-','LineWidth',0.5);
title('Spettro del segnale');
xlabel('f [Hz]');
ylabel("U(f)");