
%%initialisation des constantes
clear all;
close all;
clc;

[y,fs] = audioread('mix.wav');
NbEch=length(y);

te=1/fs;
tps=(0:NbEch-1)*te;
df=fs/NbEch;
freq =(0:NbEch-1)*df;
dt = 1/df;

Nwin = 1024; 
Nfft = Nwin;
Nhop = Nwin/2;
%%
%%Calcul de la TFCT de y
[Z,dtv,freq2] = TFCT(Nwin,Nhop, Nfft, y,fs);

Y = abs(Z);
 figure;
 imagesc(dtv,freq2,20*log10((Y)),[-10 10]);
axis xy
xlabel('Temps(s)');
ylabel('Frequence (Hz)')
title ('Spectrogramme du signal avec bruit (methode manuelle)')
ylim([0 fs/2])
colorbar 

%% Identification de la signature spectrale du bruit 
trame = Y(:,1:14);
moyen_trame = mean(abs(trame),2);

figure;
df2=fs/Nfft;
freq3 = (-Nfft/2:Nfft/2-1)*df2;
plot (freq3,moyen_trame);
xlim ([0 fs/2])
title ('Spectre du bruit')
%% Appliation de la methode de la soustraction spectrale 
matrice_bruit = zeros(size(Y));
for i=(1:(NbEch-Nwin)/Nhop)
    
    matrice_bruit(:,i) = moyen_trame;
end 
  
phase = angle(Z);
Y2 = (Y - matrice_bruit).*exp(phase*j); %reconstruction du signal y



 
 %%Le nouveau signal debruité 
figure;
imagesc (dtv,freq2,20*log10(abs(Y2)),[-10 10]);
axis xy
xlabel('Temps(s)');
ylabel('Frequence (Hz)')
title ('Spectrogramme du signal sans bruit (methode manuelle)')
ylim([0 fs/2])
colorbar 

[y2, t] = itfct(Y2, Nwin, Nhop, Nfft, fs);
y2 = real(y2);
%%negligeable = max(imag(y2));
%%rsb = db(rms(y2)./rms(y(1:length(y2))'-y2));
soundsc (y,fs)
pause(10)
soundsc (y2,fs)

 