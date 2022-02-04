clear variables; 

%FREQUENZE DI CAMPIONAMENTO E DI NYQUIST
fs = 2500; %frequenza di campionamento
ts = 1/fs; %tempo di campionamento
fnyq = fs/2; %frequenza di nyquist

%DATASET
files = ["filocaldo16cm","filocaldo17cm","filocaldo18cm","filocaldo24cm",...
    "filocaldo25cm","filocaldo26.5cm","filocaldo26cm","filocaldo27.5cm","filocaldo27cm",...
    "filocaldo28.5cm","filocaldo28cm","filocaldo29cm"];
N_samples = length(files); %numero di campioni ottenibili
wake_frequencies = zeros(1,N_samples);%vettore con i campioni di frequenza

%CICLO SU OGNI FILE PER CALCOLARE LE FREQUENZE MEDIE
count = 1;
for file_name = files
    
    filocaldo = importfile(file_name, [1, Inf]);
    data = filocaldo.VarName2;
  
    k = 3.65;
    h = 7.83;
    data = (data.^2 + k).^2 / h^2; %da cambiare con la curva di taratura corretta
    data = data - mean(data);
    
    a = 0.5;
    w = a - (1-a)*cos(2*pi*(1:1:length(data))/length(data))'; %windowing function
    data = data .* w; %windowed signal

    T = ts*(length(data)-1); %tempo
    dft = abs(fft(data)); %trasformata di fourier discreta del segnale

    df = 1/T(end); %spaceing dello spazio delle frequenze
    freq = df*(0:1:length(data)-1); %spazio delle frequenze


    freq = freq(freq<fnyq); %seleziono solo le frequenze inferiori a quelle di nyquist
    dft = dft(freq<fnyq); %seleziono i dati alle frequenze inferiori a quella di nyquist
    [~,i] = max(dft); %cerco la frequenza massima
    wake_frequencies(count) = i*df;
    count = count + 1;
end

%ANALISI STATISTICA DELLE FREQUENZE 
nu = N_samples -1; %gradi di libertà della students't
f_mean = mean(wake_frequencies); %frequenza media
f_std = std(wake_frequencies)/sqrt(N_samples); %deviazione standard

freq_space = linspace(f_mean-5,f_mean+5,1000); 
t_var = (freq_space - f_mean)./(f_std);

students = @(t) gamma((nu+1)/2)./(sqrt(nu*pi)...
    *gamma(0.5*nu).*(1+t.^2 ./ nu).^(0.5*(nu+1)));

pdf_t = students(t_var);
figure(1); hold on; grid on;
plot(freq_space,pdf_t./f_std,'k-','LineWidth',2);
plot([f_mean,f_mean],[0,max(pdf_t)/f_std],'r--','LineWidth',2);
legend('PDF','<f>=53.8792 Hz');
title('Distribuzione di probabilità per la media delle frequenze');
xlabel('f [Hz]');
ylabel('p(f)');
hold off;

figure(2);
histogram(wake_frequencies,8);
title('istogramma delle frequenze rilevate');
xlabel('f [Hz]');


conf_interval = 0.95; %una deviazione standard
t_conf = fzero(@(z) integral(@(x) students(x), -z, z) - conf_interval,1);
f_err = t_conf * f_std;


%ANALISI STATISTICA DELLA PRESSIONE E DELL'ERRORE SULLA VELOCITÀ
Ep_all=[]; % inizializzo vettore con tutte le Ep registrate
i=1;
for file_name = files
    
    % Import delle pressioni e calcolo media
    pres = importpressure(file_name, [1, Inf]);
    
    % Costruisco il vettore con tutte le Ep registrate
    Ep_all=[Ep_all; pres.VarName1]; 
    
    % Aggiorno contatore
    i=i+1;
end
N_press = length(Ep_all);
ep_mean = mean(Ep_all);
ep_std = std(Ep_all)/sqrt(N_press);

figure(3); hold on; grid on;
N_bins = floor(1.87*(N_press-1)+1);
histogram(Ep_all,N_bins,'Normalization','pdf');
gaussian = @(z) (1.0./sqrt(2*pi)).*exp(-0.5*z.^2);
z_conf = fzero(@(x) integral(gaussian, -x, x) - conf_interval,1.96);
title('Istogramma dei valori di E_p');
xlabel('ep [V]');
ylabel('p(ep)');


figure(4); hold on; grid on;
ep_space = linspace(ep_mean-0.001,ep_mean+0.001,1000);
z_var = (ep_mean-ep_space)./ep_std;
plot(ep_space,gaussian(z_var)./ep_std,'k-','LineWidth',2);
plot([ep_mean,ep_mean],[0,max(gaussian(z_var))/ep_std],'r--','LineWidth',2);
legend('PDF','<ep>=2.709 V');
title('Distribuzione di probabilità per la media delle pressioni');
xlabel('ep [V]');
ylabel('p(ep)');
hold off;

ep_err=ep_std*z_conf;

k_tras = 494.6;
k_err = 0.05*0.95;

ep_zero = 2.533;
ep_zero_err = ep_err;
P_mean = 102000;
P_err = 500*0.95;

T_mean = 23+273.15;
T_err = 1*0.95;

R_gas = 287.05;

%DENSITÀ
syms rho(T,P)
rho = P/(R_gas*T);
rho_T = subs(diff(rho,T),[T,P],[T_mean,P_mean]);
rho_P = subs(diff(rho,P),[T,P],[T_mean,P_mean]);
rho_mean = eval(subs(rho,[T,P],[T_mean,P_mean]));
rho_err = eval(sqrt((rho_T*T_err)^2 + (rho_P*P_err)^2));

%VISCOSITÀ
mu_mean = 1.827e-5; % Viscosit� dinamica aria a 23�C e 1 atm
mu_err= (1e-8)*0.95/2;

%RSS SULLA VELOCIT�

syms U(ep,ep0,k,rho) % Definisco la velocit� vista come funzione simbolica, dipende dai parametri tra parentesi
U = sqrt(2*k*(ep-ep0)/rho); %dalla pressione dinamica

%ora scriviamo le derivate, subs valuta la funzione U_k inserendo i valori
%medi

U_k = subs(diff(U,k),[ep,ep0,k,rho],[ep_mean,ep_zero,k_tras,rho_mean]); % dU/dk
U_ep = subs(diff(U,ep),[ep,ep0,k,rho],[ep_mean,ep_zero,k_tras,rho_mean]); % dU/dep
U_ep0 = subs(diff(U,ep0),[ep,ep0,k,rho],[ep_mean,ep_zero,k_tras,rho_mean]); % dU/dep0
U_rho = subs(diff(U,rho),[ep,ep0,k,rho],[ep_mean,ep_zero,k_tras,rho_mean]); % dU/drho

% Trovo U media valutando U inserendo i parametri medi, eval fa la
% valutazione della funzione
U_mean = eval(subs(U,[ep,ep0,k,rho],[ep_mean,ep_zero,k_tras,rho_mean]));
% Trovo l'errore con RSS di U rispetto alle tre variabili ep,ep0,k
U_err = eval(sqrt((U_k*k_err)^2+(U_ep*ep_err)^2+(U_ep0*ep_zero_err)^2+(U_rho*rho_err)^2)); 

fprintf('La velocità di galleria è %f +- %f con intervallo di confidenza pari al %d%% \n', ...
         U_mean, U_err, conf_interval*100);

%% ANALISI STATISTICA DEL NUMERO DI STROUHAL (RSS)

syms St(f,U) % Definisco Strouhal come equazione simbolica

c = 0.10; % Corda del profilo
cvert = c*sin(deg2rad(25)); % Lunghezza caratteristica del profilo a 25�

St = f*cvert/U;
St_U = subs(diff(St,U),[f,U],[f_mean,U_mean]); % dSt/dU
St_f = subs(diff(St,f),[f,U],[f_mean,U_mean]); % dSt/df

% Trovo Strouhal medio inserendo i parametri medi
St_mean = eval(subs(St,[f,U],[f_mean,U_mean]));
% Trovo l'errore con RSS di St rispetto a U,f
St_err = eval(sqrt((St_U*U_err)^2+(St_f*f_err)^2));

fprintf('Il  numero di Strouhal è %f +- %f con intervallo di confidenza pari al %d%% \n', ...
         St_mean, St_err, conf_interval*100);

%% ANALISI STATISTICA DEL NUMERO DI REYNOLDS  (RSS)

syms Re(U,rho,mu) % Definisco Reynolds come equazione simbolica

Re = rho*c*U/mu;
Re_U = subs(diff(Re,U),[U,rho,mu],[U_mean,rho_mean,mu_mean]); % dRe/dU
Re_rho = subs(diff(Re,rho),[U,rho,mu],[U_mean,rho_mean,mu_mean]); % dRe/dU
Re_mu = subs(diff(Re,mu),[U,rho,mu],[U_mean,rho_mean,mu_mean]); % dRe/dU

% Trovo Reynolds medio inserendo i parametri medi
Re_mean = eval(subs(Re,[U,rho,mu],[U_mean,rho_mean,mu_mean]));
% Trovo l'errore con RSS di St rispetto a U
Re_err = eval(sqrt((Re_U*U_err)^2+(Re_rho*rho_err)^2+(Re_mu*mu_err)^2)); 

fprintf('Il numero di Reynolds è %f +- %f con intervallo di confidenza pari al %d%% \n', ...
         Re_mean, Re_err, conf_interval*100);