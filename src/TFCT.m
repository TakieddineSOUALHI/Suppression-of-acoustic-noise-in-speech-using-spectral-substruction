function [X,Y,Fq] = TFCT(Nwin,Nhop, Nfft, y,fs)


Te=1/fs;%% Periode d'échantillonage
Nb_Echantillon=length(y); %% nombres d'échantillon
X=zeros(Nwin,1);

w=hamming(Nwin); % fenêtre de hamming
y_0=y(1: Nfft); % vecteur de Nfft valeur
ywin_0=y_0.*w;

Ywin_0=fft(ywin_0);
X=Ywin_0;

for i=(1:(Nb_Echantillon-Nwin)/Nhop)
    y_hop=y(1+i*Nhop : Nfft+i*Nhop);
    y_win=y_hop.*w;
    Y_win=(fft(y_win));
    X=[X Y_win];
end 


Ni=((Nb_Echantillon-Nwin)/Nhop);
dure=(Nb_Echantillon-1)*Te;

delta_t=dure/(Ni+1);
Y=(1:(Ni+1))*delta_t;
df=fs/Nfft;

Fq=(1:Nfft)*df;