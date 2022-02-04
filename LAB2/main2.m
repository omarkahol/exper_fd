clear variables;

% FREQUENZE DI CAMPIONAMENTO E DI NYQUIST

fs = 2500; %frequenza di campionamento
ts = 1/fs; %tempo di campionamento
fnyq = fs/2; %frequenza di nyquist

% DATASET

files = ["filocaldo16cm","filocaldo17cm","filocaldo18cm",...
         "filocaldo24cm","filocaldo25cm","filocaldo26.5cm",...
         "filocaldo26cm","filocaldo27.5cm","filocaldo27cm",...
         "filocaldo28.5cm","filocaldo28cm","filocaldo29cm"]; %tutti i file in un vettore
     
N_samples = length(files); % numero di campioni ottenuti
wake_frequencies = zeros(1,N_samples);% inizializzo vettore delle frequenze

%% IMPORT DELLE FREQUENZE DAI FILE ACQUISITI

count = 1;
for file_name = files
    
    % Import dal file
    filocaldo = importfile(file_name, [1, Inf]);
    data = filocaldo.VarName2; % vettore dalla II colonna (letture HWA)
  
    data = ((3.65+(5-data).^2)/7.83).^2; % curva di taratura corretta HWA
    data = data - mean(data); % pulizia del segnale dal valor medio
    
    % Windowing del segnale
    a = 0.5;
    w = a-(1-a)*cos(2*pi*(1:1:length(data))/length(data))'; %Hanning window
    data = data .* w; % segnale windowed 

    % FFT, correzioni dell'output, salvataggio output
    T = ts*(length(data)-1); % Calcolo tempo di acquisizione
    dft = abs(fft(data)); % trasformata di Fourier discreta del segnale

    df = 1/T(end); % Risoluzione in frequenza 
    freq = df*(0:1:length(data)-1); % spazio delle frequenze osservabili

    freq = freq(freq<fnyq); % seleziono solo le f. inferiori a f. Nyquist
    dft = dft(freq<fnyq); % seleziono output a f. inferiori a f. Nyquist
    [~,i] = max(dft); % cerco la frequenza massima
    wake_frequencies(count) = i*df; % inserisco la max f. nel vettore
    
    % Aggiorno iterazione
    count = count + 1;
end

figure(1)
plot(w);
hold on
title('Fuzione di finestratura di Hanning');

figure(2)
plot(data);
title('Dati finestrati (da file 29cm)');

figure(3)
plot(dft);
title('FFT (da file 29cm)');

%% ANALISI STATISTICA DELLE FREQUENZE USANDO DISTRIBUZIONE T-STUDENTS

nu = N_samples -1; % Gradi di libert� della T-Students
f_mean = mean(wake_frequencies); % Frequenza media
f_sigma = std(wake_frequencies); % Deviazione standard 

conf_interval = 0.95; % Scelgo intervallo di confidenza del 95%
t_value = 1.796; % Dalle tabelle della t-students (.95 con nu=11)
f_err = t_value*f_sigma/sqrt(N_samples); % Errore sulla frequenza con distribuzione t-students

fprintf('La frequenza di shedding � %f +- %f con intervallo di confidenza pari al %d%% \n', ...
         f_mean, f_err, conf_interval*100);


% Creo la distribuzione T-Student per plottarla

% Vettori
freq_space = linspace(f_mean-10,f_mean+10,1000); 
t_var = (freq_space - f_mean)./f_sigma;

% Funzione T-Students scritta in Symbolic math 
students = @(t) gamma((nu+1)/2)./(sqrt(nu*pi)...
           *gamma(0.5*nu).*(1+t.^2 ./ nu).^(0.5*(nu+1)));

pdf_t = students(t_var); % Creazione distribuzione T-Students

figure(4);
hold on;
grid on;
plot(freq_space,pdf_t,'k-','LineWidth',2);
plot([f_mean,f_mean],[0,max(pdf_t)],'r--','LineWidth',2);
leg = sprintf('Valor medio = %.3f ',f_mean);
legend('t-Student PDF dei dati',leg);
title('PDF of the wake frequencies');
xlabel('f [Hz]');
ylabel('pdf');
hold off;



%% IMPORT DELLE LETTURE DEL TRASDUTTORE SULLE PRESSIONI DAI FILE ACQUISITI

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

%% ANALISI STATISTICA DELLA PRESSIONE E DELL'ERRORE SULLA VELOCIT� (RSS)

n = numel(Ep_all); % Dimensione di Ep_all

% Creo 100 campioni da Ep_all per poter usare una distribuzione gaussiana
Ep_resample = mat2cell(Ep_all,diff([0:floor(n/100):n-1,n])); %da matrice a cella

% Inizializzo e riempio il vettore con le medie dei 100 campioni
mean_Ep = zeros(1,100); 
for i=1:100
    mean_Ep(i) = mean(Ep_resample{i,1}); %vettore con le medie dei 100 campioni
end

ep_mean = mean(mean_Ep); % Media delle medie dei 100 campioni
ep_sigma = std(mean_Ep); % Deviazione standard delle medie dei 100 campioni

%% section ev
figure(5); hold on; grid on;
histogram(mean_Ep,10,'Normalization','pdf');
ep_space = linspace(ep_mean-0.005,ep_mean+0.005,1000);
ep_dist = (1.0./(sqrt(2*pi*ep_sigma^2))).*exp(-(ep_mean-ep_space).^2 ./ (2*ep_sigma^2));
plot(ep_space,ep_dist,'k--','LineWidth',2);
title('PDF dei valori di E_p');
xlabel('v');
ylabel('p(v)');
legend('histogram (pdf normalized)','pdf')
ep_err=ep_sigma*1.96; % Int. di conf. = 95% (2*sigma per dist. gaussiana)
%% section2
ep_zero = 2.533; % Azzeramento del trasduttore
ep_zero_err = ep_err; % Errore sull'azzeramento (= errore sulle letture) ASSUNZIONE 

k_tras = 494.6; % Costante di calibrazione del trasduttore
k_err = 0.1*0.95/2; % Errore sulla costante di calibrazione, assumendo distribuzione uniforme
%errore come ultima cifra significativa

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

fprintf('La velocit� di galleria � %f +- %f con intervallo di confidenza pari al %d%% \n', ...
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

fprintf('Il  numero di Strouhal � %f +- %f con intervallo di confidenza pari al %d%% \n', ...
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

fprintf('Il numero di Reynolds � %f +- %f con intervallo di confidenza pari al %d%% \n', ...
         Re_mean, Re_err, conf_interval*100);
