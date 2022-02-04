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

x0 = [1 1 1 1]; 
fitfun = fittype(@(a,b,c,d,x) a*x.^8 + b*x.^6 + c*x.^4 + d*x.^2 + 1.0);
[fitted_curve,gof] = fit(deg,f,fitfun,'StartPoint',x0);

a = linspace(min(deg),max(deg),1000);
plot(a,fitted_curve(a),'g-','LineWidth',2);











