function [x, t] = itfct(X_m, wlen, hop, nfft, fs)

% [x, t] = itfct(stft, wlen, hop, wlen)
% transformee de Fourier a court-terme inverse
% par overlapp-add (OLA)
%
% ENTREES
% X_m       matrice TFCT (FxT)
% wlen      taille de la fenetre d'analyse (en echantillons)
% hop       pas d'avancement (en echantillons)
% nfft      nombre de points frequentiels de la FFT
% fs        frequence d'echnatillonage (Hz)
%
% SORTIES
% x_v       signal temporel (1xT)
% t_v       temps correspondants (s)

% initialize variables
coln    = size(X_m, 2);
xlen    = nfft + (coln-1)*hop;
x       = zeros(1, xlen);
x2      = x;

% form a hamming window
win = hamming(wlen);

% perform IFFT and weighted-OLA
for b = 0:hop:(hop*(coln-1))
        
    % compute inverse FFT
    xprim = ifft(X_m(:, 1+b/hop), wlen);

    % compute weighted-OLA
    x((b+1):(b+wlen)) = x((b+1):(b+wlen)) + (xprim.*win)';
    x2((b+1):(b+wlen)) = x2((b+1):(b+wlen)) + (win.^2)';
    
end

% normalize OLA by window
x = x./x2;

% calculate the time vector
t = (0:length(x)-1)/fs;

end