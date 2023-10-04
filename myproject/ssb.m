clc; %clears all the text from the Command Window
close all; % closes all open MATLAB figure windows
clear; %delete data in Workspace
%1))
%Reading from audio file
[mt,Fs]=audioread('eric.wav');
Nsamps=length(mt);
sound(mt,Fs);
pause(Nsamps/Fs);
t =linspace(0,Nsamps/Fs,Nsamps); %Nsamps/Fs=411248/48000 =8.5668=endpoint      
f =linspace(-Fs/2,Fs/2,Nsamps);
figure;
subplot(2, 1, 1)
plot(t, mt);
grid on;
title('msg in time domain');
subplot(2, 1, 2);
plot(f, abs(fftshift(fft(mt))));
grid on;
title('msg in freq domain');
%--------------------------------------------------------------------------
%2))
% remove all freq above 4000
%2*BW=8000
d = designfilt('lowpassfir', 'FilterOrder', 8, 'CutoffFrequency',4000, 'SampleRate', Fs);
mtnew = filter(d, mt);
%--------------------------------------------------------------------------
%3))

figure;
subplot(2, 1, 1);
plot(t, mtnew);
grid on;
title('Filtered msg in time domain');

subplot(2, 1, 2);
plot(f, abs(fftshift(fft(mtnew))));
grid on;
title('Filtered msg in freq domain');

%sound(mtnew, Fs);
%--------------------------------------------------------------------------
%4))
%Generate a DSB-SC modulated signal
Fc = 100000;
re_mt = resample(mtnew, 5 * Fc, Fs); 
Fs = 5*Fc;
t = linspace(0, length(re_mt)/Fs, length(re_mt));
f = linspace(-Fs/2, Fs/2, length(re_mt));
ct = cos(2*pi*Fc*t);
DSBSC = re_mt .* transpose(ct);

figure; 
subplot(2,1,1);
plot(t, DSBSC);
grid on;
title('DSB-SC Time Domain')

subplot(2,1,2); 
plot(f, abs(fftshift(fft(DSBSC))));
grid on;
title('DSB-SC Frequency Domain')
%--------------------------------------------------------------------------
%5))
%SSB SC From DSBSC Modulation(LSB)

d = designfilt('lowpassfir', 'FilterOrder', 8000, 'CutoffFrequency', Fc, 'SampleRate', Fs);
DSBSC = filter(d, DSBSC);

figure; 
subplot(2,1,1);
plot(t, DSBSC);
grid on;
title('SSB-SC Time Domain')
subplot(2,1,2); 
plot(f, abs(fftshift(fft(DSBSC))));
grid on;
title('SSB-SC Frequency Domain')
%--------------------------------------------------------------------------
%6))
% Coherent Detection

ct = cos(2*pi*Fc*t);
deSignal = DSBSC .* transpose(ct);
d = designfilt('lowpassfir', 'FilterOrder', 8000, 'CutoffFrequency', 4000, 'SampleRate', 5 * Fc);
coherentSC  = filter(d, deSignal);

%sound(coherentSC, Fs);

figure;
subplot(2, 1, 1);
plot(t, coherentSC);
title('SC Coherent time domain');
grid on;
subplot(2, 1, 2);
plot(f, abs(fftshift(fft(coherentSC))));
xlim([-5000,5000])
grid on;
title('SC Coherent freq domain');
%--------------------------------------------------------------------------
%6))

% Butterworth filter

ct = cos(2*pi*Fc*t);
deSignal = DSBSC .* transpose(ct);
[b, a] = butter(3, Fc/(Fc * 5 / 2));
butteredSig = filter(b, a, deSignal);
butteredSigmod = filter(b, a, DSBSC);

figure;
subplot(2, 1, 1);
plot(t, butteredSigmod);
grid on;
title('SC modulated signal (Butterworth) in time domain');
subplot(2, 1, 2);
plot(f, abs(fftshift(fft(butteredSigmod))));
grid on;
title('SC modulated signal Spectrum (Butterworth) in freq domain');


figure;
subplot(2, 1, 1);
plot(t, butteredSig);
grid on;
title('SC Coherent (Butterworth) in time domain');
subplot(2, 1, 2);
plot(f, abs(fftshift(fft(butteredSig))));
xlim([-5000,5000])
grid on;
title('SC Coherent Spectrum (Butterworth) in freq domain');
%--------------------------------------------------------------------------
%8))

coherentSNR0 = awgn(coherentSC, 0);
coherentSNR10 = awgn(coherentSC, 10);
coherentSNR30 = awgn(coherentSC, 30);

%sound(coherentSC0SNR, Fs);
%sound(coherentSC10SNR, Fs);
%sound(coherentSC30SNR, Fs);

figure;
subplot(2, 1, 1);
plot(t, coherentSNR0);
grid on;
title('SC Coherent 0 SNR in time domain');
subplot(2, 1, 2);
plot(f, abs(fftshift(fft(coherentSNR0))));
grid on;
title('SC Coherent 0 SNR in freq domain');

figure;
subplot(2, 1, 1);
plot(t, coherentSNR10);
grid on;
title('SC Coherent 10 SNR in time domain');
subplot(2, 1, 2);
plot(f, abs(fftshift(fft(coherentSNR10))));
grid on;
title('SC Coherent 10 SNR in freq domain');

figure;
subplot(2, 1, 1);
plot(t, coherentSNR30);
grid on;
title('SC Coherent 30 SNR in time domain');
subplot(2, 1, 2);
plot(f, abs(fftshift(fft(coherentSNR30))));
grid on;
xlim([-5000,5000])
title('SC Coherent 30 SNR in freq domain');

% Transmitted Carrier Modulation

ct = cos(2*pi*Fc*t);
re_mt = re_mt + (max(re_mt) * 2);
TC_Signal = re_mt .* transpose(ct); 

d = designfilt('lowpassfir', 'FilterOrder', 8000, 'CutoffFrequency', Fc, 'SampleRate', Fs);
TC_Signal = filter(d, TC_Signal);

figure;
subplot(2, 1, 1);
plot(t, TC_Signal);
grid on;
title('DSB-TC Time Domain');
subplot(2, 1, 2);
plot(f, abs(fftshift(fft(TC_Signal))));
grid on;
title('DSB-TC Frequency Domain');
%--------------------------------------------------------------------------
%9))
% Envelope Detection

envelopeTC = abs(hilbert(TC_Signal)); % DSB-TC envelope detection

%sound(envelopeTC, Fs);

figure;
subplot(2, 1, 1);
plot(t, envelopeTC);
ylim([0 0.3])
grid on;
title('TC Envelope in time domain');
subplot(2, 1, 2);
plot(f, abs(fftshift(fft(envelopeTC))));
grid on;
title('TC Envelope in freq domain');