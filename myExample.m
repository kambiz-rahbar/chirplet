clc
clear
close all

Fs = 1000;    % samples per second
dt = 1/Fs;    % seconds per sample
StopTime = 1; % seconds
Time = (0:dt:StopTime)';

Sig = chirp(Time, 50, StopTime, 450, 'quadratic');

N = 15;
WinLen = 100;
Chirplet = GLCT(Sig, N, Fs, WinLen);
ChirpletPowerSpect = abs(Chirplet).^2;

figure(1);

subplot(2, 2, 1); plot(Time, Sig);
axis square
xlabel('Time (Sec)')
ylabel('Amplitude');
title('Signal in time domain');

dF = Fs/length(Time); % hertz
Freq = 0:dF:Fs/2;     % hertz

[XX, YY] = meshgrid(Time, Freq);
subplot(2, 2, 2); mesh(XX, YY, ChirpletPowerSpect);
axis xy square
ylabel('Freq (Hz)');
xlabel('Time (Sec)')
title(sprintf('Chirplet PowerSpect (N = %d)', N));

subplot(2, 2, 3); imagesc(Time, Freq, ChirpletPowerSpect);
axis xy square
ylabel('Freq (Hz)');
xlabel('Time (Sec)')
title('Chirplet PowerSpect (TOP view)');
