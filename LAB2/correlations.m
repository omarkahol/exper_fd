clear variables;
clc;

file_name = "filocaldo28cm";
hwa_data = importfile(file_name, [1, Inf]);
data = hwa_data.VarName2; % vettore dalla II colonna (letture HWA)
data = ((3.65+(5-data).^2)/7.83).^2;

a=0.5;
w = a - (1-a)*cos(2*pi*(1:1:length(data))/length(data))';
[R,tau] = xcorr(data .* w);
[R_max, i_max] = max(R);

fs = 2500; %frequenza di campionamento
ts = 1/fs; %tempo di campionamento

tau = tau*ts;
fnyq = fs/2; %frequenza di nyquist
T = ts*(length(data)-1);
df = 1/T(end);

R_vary = xcorr(data-mean(data));
R_positive = R_vary(floor(0.5*length(R_vary)):end);
psd = fft(R_positive-mean(R_positive));
freq = df*(0:1:length(data)-1);

psd = psd(freq<fnyq);
freq = freq(freq<fnyq);
[psd_max, i_psd_max] = max(abs(psd));
phi = angle(psd(i_psd_max));
f = freq(i_psd_max);
dt = 1/f;
skip = floor(dt/ts);
search = 200;
x_data = zeros(1,2*search-1);
y_data = zeros(1,2*search-1);
for i=1:search
    [mm,ii]=max(R(i_max+i*skip:i_max+(i+1)*skip));
    x_data(i) = tau(ii+i_max+i*skip-1);
    y_data(i) = mm;
    [mm,ii]=max(R(i_max-(i+1)*skip:i_max-i*skip));
    x_data(2*search-i) = tau(i_max-(i+1)*skip + 1 +ii);
    y_data(2*search-i) = mm;
end

p = polyfit(x_data,y_data,4);
AR = polyval(p,0);
figure(1); hold on; grid on;
plot(x_data,y_data,'r*');
plot(tau,R,'b-','LineWidth',2);
plot(tau(i_max-5000:i_max+5000),polyval(p,tau(i_max-5000:i_max+5000)),'g-','LineWidth',2);
xlabel('s');
ylabel('R');
title('Correlazione');
legend('Inviluppo','R(s)','Fit')
hold off;

Sn2 = R_max-AR;
SNR = 10*log10(AR/Sn2)

% data = data - mean(data);
% psd = fft(data .* w);
% psd = psd(freq<fnyq);

freq_oct = log2(freq(freq>380 & freq < 1000)');
db_psd = 10*log10(abs(psd(freq>380 & freq < 1000)));
x0 = [1 1]; 
fitfun = fittype(@(a,b,x) a*x+b);
[fitted_curve,gof] = fit(freq_oct,db_psd,fitfun,'StartPoint',x0);

figure(2); hold on; grid on;
plot(log2(freq),10*log10(abs(psd)),'b-');
plot(freq_oct,fitted_curve(freq_oct),'r-','LineWidth',5)
title('Power Spectral Density');
xlabel('log2(f)')
ylabel('10*log10(|psd|)');
legend('PSD','Fitted Psd');

hold off;





