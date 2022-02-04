file_name = "filocaldo28cm";
hwa_data = importfile(file_name, [1, Inf]);
data = hwa_data.VarName2; % vettore dalla II colonna (letture HWA)
data = ((3.65+(5-data).^2)/7.83).^2;

a=0.5;
w = a - (1-a)*cos(2*pi*(1:1:length(data))/length(data))';
[R,tau] = xcorr(data .* w);

fs = 2500; %frequenza di campionamento
ts = 1/fs; %tempo di campionamento

tau = tau*ts;
fnyq = fs/2; %frequenza di nyquist
T = ts*(length(data)-1);
df = 1/T(end);

R_positive = R(floor(0.5*length(R)):end);
psd = fft(R_positive-mean(R_positive));
freq = df*(0:1:length(data)-1);

psd = psd(freq<fnyq);
freq = freq(freq<fnyq);


freq_oct = log2(freq(freq>380 & freq < 1000)');
db_psd = 10*log10(abs(psd(freq>380 & freq < 1000)));
x0 = [1 1]; 
fitfun = fittype(@(a,b,x) a*x+b);
[fitted_curve,gof] = fit(freq_oct,db_psd,fitfun,'StartPoint',x0);

figure(1); hold on; grid on;
plot(log2(freq),10*log10(abs(psd)),'b-');
plot(freq_oct,fitted_curve(freq_oct),'r-','LineWidth',5)
title('Power Spectral Density');
xlabel('log2(f)')
ylabel('10*log10(|psd|)');
legend('PSD','Fitted Psd');

hold off;