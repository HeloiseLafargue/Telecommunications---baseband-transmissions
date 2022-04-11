% ------------------------------------------------------------------------
% Télécommunications : Etudes de chaines de transmission en bande de base                
% ------------------------------------------------------------------------

clear
close all

Fe = 24000;
Te = 1 / Fe;
Rb = 3000;
Tb = 1 / Rb;

%% ____ Etude des interférences entre symbole et du crière de Nyquist _____

% _______________ Etude sans canal de propagation _________________________

%------------------------- Modulation -------------------------------------

n = 1;          
Ns = n * Tb / Te;
Ts = Ns * Te; 
nb_bits=1000;
bits=randi([0,1],1,nb_bits);
Symboles=2*bits-1;
Suite_diracs=kron(Symboles, [1 zeros(1,Ns-1)]);
h=ones(1,Ns);
x=filter(h,1,Suite_diracs);
figure(5); sgtitle('Etude sans canal de propagation');
subplot(2,1,1); plot(x); title("Signal modulé");
axis([0 nb_bits-1 -1.5 1.5]); xlabel("Temps (s)"); ylabel("signal");

%------------------------ Démodulation ------------------------------------

z = filter(h,1,x);
figure(5); subplot(2,1,2); plot(z); title("Signal démodulé");
axis([0 nb_bits-1 -1.5 1.5]); xlabel("Temps (s)"); ylabel("signal");


%---------------------- Réponse impulsionnelle ----------------------------
g = conv(h,h);
figure(6); plot(g); xlabel('Temps (s)');
title("Réponse impulsionelle globale de la chaîne de transmission");
% c'est optimal pour n0 = 8

%---------------------- Diagramme de l'oeil -------------------------------

figure(7); plot(reshape(z,Ns,length(z)/Ns));
title("Diagramme de l'oeil");

% Pour ne pas avoir d'interférences on se place là où il n'y que 2 points :
% à n0 = 8s. On retrouve la même valeur qu'avec la réponse impusionnelle.

%------------------------ Echantillonage ----------------------------------

n0 = 8;
echant = z([n0:8:length(z)]);
TEB = length(find(abs(echant/8 - Symboles) > 0)) / length(Symboles);
fprintf("4.2.7) Pour n0 = 8 (optimal), le TEB obtenu est %f \n",TEB);
% TEB = 0 : la transmission est totale.

n0_2 = 3;
echant_2 = z([n0_2:8:length(z)]);
TEB = length(find(abs(echant_2/8 - Symboles) > 0)) / length(Symboles);
fprintf("4.2.8) Pour n0 = 3, le TEB obtenu est %f \n",TEB);
% TEB = 0.501 : la transmission n'est pas parfaite à cause des
% interférences.


% ____________ Etude avec canal de propagation sans bruit _________________

%---------------------- Cas BW = 8000 Hz ----------------------------------

fc = 8000;
hc = (2*fc/Fe)*sinc(2*(fc/Fe)*[-(Ns-1)/2 : (Ns-1)/2]);

g = conv( h, conv(hc,h));
figure(8); sgtitle('Cas BW = 8000 Hz'); subplot(3,1,1);
plot(g); title("Réponse impulsionelle"); xlabel('Temps (s)');

z = filter(h,1,filter(h,1,Suite_diracs));
figure(8); subplot(3,1,2);
plot(reshape(z,Ns,length(z)/Ns)); title("Diagramme de l'oeil");

% les réponses fréquentielles
F1 = linspace(-fc, fc, length(fft(g, 2^15))); 
F2 = linspace(-fc, fc, length(fft(h, 2^15)));
F3 = linspace(-fc, fc, length(fft(conv(h,h), 2^15))); 

figure(8); subplot(3,1,3);
plot(F1, fftshift(abs(fft(g, 2^15)))); % |Hc|
hold on;
plot(F2, fftshift(abs(fft(conv(h,h), 2^15)))); % |H*Hr|
xlabel("Fréquence (Hz)"); ylabel("Y(f)");
title("Réponse fréquentielle"); legend('|Hc(f)|', '|H(f)Hr(f)|');

% échantillonage et TEB
n0 = 8;
echant = z([n0:8:length(z)]);
TEB = length(find(abs(echant/8 - Symboles) > 0)) / length(Symboles);
fprintf("4.3.1) Pour n0 = 8 (optimal), le TEB obtenu est %f \n",TEB);
% TEB = 0 : la transmission est totale.


%---------------------- Cas BW = 1000 Hz ----------------------------------

fc = 1000;
hc = (2*fc/Fe)*sinc(2*(fc/Fe)*[-(Ns-1)/2 : (Ns-1)/2]);

g = conv( h, conv(hc,h));
figure(9); sgtitle('Cas BW = 1000 Hz'); subplot(3,1,1);
plot(g); title("Réponse impulsionelle"); xlabel('Temps (s)');

z = filter(h,1,filter(h,1,Suite_diracs));
figure(9); subplot(3,1,2);
plot(reshape(z,Ns,length(z)/Ns)); title("Diagramme de l'oeil");

% les réponses fréquentielles
F1 = linspace(-fc, fc, length(fft(g, 2^15))); 
F2 = linspace(-fc, fc, length(fft(h, 2^15)));
F3 = linspace(-fc, fc, length(fft(conv(h,h), 2^15))); 

figure(9); subplot(3,1,3);
plot(F1, fftshift(abs(fft(g, 2^15)))); % |Hc|
hold on;
plot(F2, fftshift(abs(fft(conv(h,h), 2^15)))); % |H*Hr|
xlabel("Fréquence (Hz)"); ylabel("Y(f)");
title("Réponse fréquentielle"); legend('|Hc(f)|', '|H(f)Hr(f)|');

% échantillonage et TEB
n0 = 8;
echant = z([n0:8:length(z)]);
TEB = length(find(abs(echant/8 - Symboles) > 0)) / length(Symboles);
fprintf("4.3.2) Pour n0 = 8 (optimal), le TEB obtenu est %f \n",TEB);
% TEB = 0 : la transmission est totale.
