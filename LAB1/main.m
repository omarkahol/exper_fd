clear variables;
data = importfile("/home/omarkahol/MATLAB/SPERIMENTALE/LAB1/Misure Speimentali Lab 1.xlsx", "Foglio1", [2, 35]);

deg = rmmissing(data.GRADI);
e1 = rmmissing(data.E1);
e2 = rmmissing(data.E2);

K1 = 497.68; %Pa/V
K2 = 400; %Pa/V

dp1 = 0.01*K1*5;
dp2 = 0.01*K2*5;

syms p1(k1)
syms p2(k2)

p1=k1.*e1;
p2=k2.*e2;
k_zerodeg = (K2*e2(deg==0))/(K1*e1(deg==0));

f = p2 ./(p1.*k_zerodeg);
err = sqrt(dp1.*diff(f,k1).^2 + dp2.*diff(f,k2).^2);

err = eval(subs(err,[k1 k2],[K1 K2]));
f = eval(subs(f,[k1 k2],[K1 K2]));

figure(1); hold on; grid on;
errorbar(deg,f,err,'k-s','MarkerSize',2,...
    'MarkerEdgeColor','red','MarkerFaceColor','red')
title('ERROR PLOT');
xlabel('alfa [°]');
ylabel('curva di taratura [-]');

dft = fft(f);
deg_length = max(deg)-min(deg);
freqs = (1/deg_length) .* (0:1:length(dft)-1)';
w = 2*pi/deg_length;
sigma = 4*w/(2*pi);
filter = exp(-0.5*freqs.^2 ./ sigma.^2);
dft = dft .* filter;
w = 2*pi/deg_length;
an = (2/length(dft)).*real(dft(1:floor(0.5*length(dft))));
an(1) = an(1) / 2;

x_space = linspace(min(deg),max(deg),1000);
y_space = zeros(1,1000);
for i = 1:length(an)
    y_space = y_space + an(i).*cos((i-1)*w.*(x_space+min(deg)));
end

plot(x_space,y_space,'r-','LineWidth',2);












