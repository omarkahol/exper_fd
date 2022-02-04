clear variables;

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 1);

% Specify range and delimiter
opts.DataLines = [15, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = "Value";
opts.VariableTypes = "double";

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
tbl = readtable("/home/omarkahol/MATLAB/SPERIMENTALE/ESERCITAZIONE1/Strouhal150ks.csv", opts);

%% Convert to output type
samples = tbl.Value;

%% Clear temporary variables
clear opts tbl

ts = 4.0e-5;
h = 5.5*10^4;
k = .94;
U = (samples.^2- k).^2 ./ h^2;
t = ts * (0:1:length(samples)-1)';

a = 0.5;
i = (0:1:length(samples)-1)';

w = a - (1-a).*cos(2.*pi.*i./length(samples));

wu = w.* U;

figure(1);
plot(t, wu, 'k-', 'LineWidth',1);
title('velocity');
xlabel('t');
ylabel('U');

T = t(end);
df = 1/T;

dft = abs(fft(wu));

f_nyq = 12.5e+3;
freq = df*i;
freq = freq(freq < f_nyq);
dft = dft(freq < f_nyq);

figure(2);
plot(freq, dft, 'k-', 'LineWidth',1);
title('dft u');
xlabel('f');
ylabel('dft');


[R,tau] = xcorr(U);
R = R(floor(0.5*length(R)):end);
R = R-mean(R);

psd = abs(fft(R));
psd = psd(freq < f_nyq);

freq_oct = log2(freq(freq>380));
db_psd = 10*log10(abs(psd(freq>380)));
x0 = [1 1]; 
fitfun = fittype(@(a,b,x) a*x+b);
[fitted_curve,gof] = fit(freq_oct,db_psd,fitfun,'StartPoint',x0);

figure(3); hold on; grid on;
plot(log2(freq),10*log10(abs(psd)),'b-');
plot(freq_oct,fitted_curve(freq_oct),'r-','LineWidth',5)
title('Power Spectral Density');
xlabel('log2(f)')
ylabel('10*log10(|psd|)');
legend('PSD','Fitted Psd');












