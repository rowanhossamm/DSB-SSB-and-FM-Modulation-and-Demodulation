clc; %clears all the text from the Command Window
close all; % closes all open MATLAB figure windows
clear; %delete data in Workspace
%1))
%Reading from audio file
[mt,Fs]=audioread('eric.wav');
Nsamps=length(mt);
clear sound;
sound(mt,Fs);
pause(Nsamps/Fs);
%samplefreq(Fs)in audio.wav=48000
%time
t = linspace(0,Nsamps/Fs,Nsamps); %Nsamps/Fs=411248/48000 =8.5668=endpoint      
f=linspace(-Fs/2,Fs/2,Nsamps);
%fftshift make midpoint=0
mf=fftshift(fft(mt));

%Plot audio File in Time Domain
subplot(2,1,1);
plot(t, mt)
xlabel('Time(s)')
ylabel('Amplitude')
title('audiowave in Time Domain')

%Plot Sound File in Frequency Domain
subplot(2,1,2);
%abs takes real part only
plot(f,abs(mf))
xlabel('Frequency(Hz)')
ylabel('Amplitude')
title('audiowave in frequency Domain')
%--------------------------------------------------------------------------
%2))
% remove all freq above 4000
d = designfilt('lowpassfir', 'FilterOrder', 8000, 'CutoffFrequency',4000, 'SampleRate', Fs);
mtnew = filter(d, mt);
% fc=4000;
% mf(f>= fc | f<=-fc) = 0;
% mtnew = ifft(ifftshift(mf));
Nsamp2 = length(mtnew);
mf = fftshift(fft(mtnew));
f = Fs/2*linspace(-1,1,Nsamp2);

%--------------------------------------------------------------------------
%3))
% plotting Time Domain
t = linspace(0,Nsamp2/Fs,Nsamp2);
figure;
subplot(2,1,1);
plot(t,real(mtnew));
title('Time domain of filtered signal');
ylim([-0.3,0.2])

% plotting freq Domain
subplot(2,1,2);
plot(f,abs(mf)/Nsamp2);
title('Spectrum of filtered signal');
%--------------------------------------------------------------------------
%4))
sound(real(mtnew),Fs);
pause(Nsamp2/Fs);
%--------------------------------------------------------------------------
%5))
%modulation
Fc=100000;
mod=0.5;
c_Fs=5*Fc;
%to resample mtnew from Fs to c_Fs which is 5*Fc
re_mt=resample(mtnew,c_Fs,Fs);
re_Nsamps=length(re_mt);
re_t = linspace(0,(re_Nsamps/c_Fs), re_Nsamps);
re_t=re_t';
ct=cos(2*pi*Fc*re_t);
DSBSC= ct.*re_mt;
DC = 2 * max(abs(re_mt));
DSBTC= DC*(1 + mod*re_mt).*ct;

%plotting time domain
figure;
subplot(2,1,1);
plot(re_t,real(DSBSC))
title ('DSB-SC: modulated  resampled filtered signal time domain')
xlabel('Time');
ylabel('Amplitude');
subplot(2,1,2);
plot(re_t,real(DSBTC));
title ('DSB-TC modulated resampled filtered signal time domain');
xlabel('Time');
ylabel('Amplitude');

%plotting frequency domain
freq_DSBSC=(fftshift(fft(DSBSC))/c_Fs);
SC_f=linspace(-c_Fs/2,c_Fs/2,length(freq_DSBSC));
figure;
subplot(2,1,1);
plot(SC_f,abs(freq_DSBSC))
title ('DSB-SC modulated resampled filtered signal freq. domain');
xlabel('Frequency');
ylabel('Amplitude');
freq_DSBTC = (fftshift(fft(DSBTC))/c_Fs);
TC_f=linspace(-c_Fs/2,c_Fs/2,length(freq_DSBTC));
subplot(2,1,2);
plot(TC_f,abs(freq_DSBTC))
title ('DSB-TC modulated resampled filtered signal freq. domain');
%--------------------------------------------------------------------------
%6))
DSBSC_envelope = abs(hilbert(DSBSC));
DSBTC_envelope = abs(hilbert(DSBTC));
figure;
subplot(2,1,1);
plot(re_t,DSBSC_envelope);
title('DSB-SC envelope detector');
xlabel('Time');
ylabel('Amplitude');

subplot(2,1,2);
plot(re_t,DSBTC_envelope);
title('DSB-TC envelope detector');
xlabel('Time');
ylabel('Amplitude');
%--------------------------------------------------------------------------
%7))
%to resample mtnew from c_Fs to Fs
re_DSBTC = resample(DSBTC_envelope, Fs, c_Fs);
re_DSBSC = resample(DSBSC_envelope, Fs, c_Fs);
sound(real(re_DSBTC),Fs);  
pause(length(re_DSBTC)/Fs);
sound(real(re_DSBSC),Fs);
pause(length(re_DSBSC)/Fs);
%observation:DSBTC modulated signal is clearer,envelope can be used with DSB_TC
%--------------------------------------------------------------------------
%8))
%demodulation and coherent detector DSB TC
SNR1 = 0;
SNR2 = 10;
SNR3 = 30;
y1 = awgn(DSBTC,SNR1); %awgn=add white Gaussian noise
y2 = awgn(DSBTC,SNR2);
y3 = awgn(DSBTC,SNR3);
envelopeSNR1 = abs(hilbert(y1));
envelopeSNR2 = abs(hilbert(y2));
envelopeSNR3 = abs(hilbert(y3));

%plotting time domain
figure;
x1 = y1.*ct;
subplot(3,1,1);
plot(re_t,[x1 DSBTC_envelope]);
title('DSBTC, SNR=0db in time domain');
xlabel('Time');
ylabel('Amplitude');
x2 = y2.*ct;
subplot(3,1,2);
plot(re_t,[x2 DSBTC_envelope]);
title('DSBTC, SNR=10db in time domain');
xlabel('Time');
ylabel('Amplitude');
x3 = y3.*ct;
subplot(3,1,3);
plot(re_t,[x3 DSBTC_envelope]);
title('DSBTC, SNR=30db in time domain');
xlabel('Time');
ylabel('Amplitude');

%plotting frequency domain
figure;
subplot(3,1,1);
freq_x1=fftshift(fft(x1));
plot(TC_f,[abs(freq_DSBTC),abs(freq_x1)]);
title('DSBTC, SNR=0db in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');
subplot(3,1,2);
freq_x2=fftshift(fft(x2));
plot(TC_f,[abs(freq_DSBTC),abs(freq_x2)]);
title('DSBTC, SNR=10db in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');
subplot(3,1,3);
freq_x3=fftshift(fft(x3));
plot(TC_f,[abs(freq_DSBTC),abs(freq_x3)]);
title('DSBTC, SNR=30 db in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');

soundsc(envelopeSNR1);
pause(1);
clear sound;
soundsc(envelopeSNR2);
pause(1);
clear sound;
soundsc(envelopeSNR3);
pause(1);
clear sound;

% demodulation and coherent detector DSB SC
fs = 1000000;
y1 = awgn(DSBSC,SNR1);
y2 = awgn(DSBSC,SNR2);
y3 = awgn(DSBSC,SNR3);

%plotting time domain
x1 = y1.*ct;
[b,a] = butter(5,Fc/(fs/2));
%filtfilt remove the high frequency noise and has no phase delay
x1 = filtfilt(b,a,x1);
subplot(3,1,1);
plot(re_t,[x1,DSBSC_envelope]);
title('DSBSC, SNR=0db in time domain');
xlabel('Time');
ylabel('Amplitude');
x2 = y2.*ct;
[b,a] = butter(5,Fc/(fs/2));
x2 = filtfilt(b,a,x2);
subplot(3,1,2);
plot(re_t,[x2,DSBSC_envelope]);
title('DSBSC, SNR=10db in time domain');
xlabel('Time');
ylabel('Amplitude');
x3 = y3.*ct;
[b,a] = butter(5,Fc/(fs/2));
x3 = filtfilt(b,a,x3);
subplot(3,1,3);
plot(re_t,[x3,DSBSC_envelope]);
title('DSBSC, SNR=30db in time domain');
xlabel('Time');
ylabel('Amplitude');

%plotting frequency domain
figure;
subplot(3,1,1);
freq_x1=fftshift(fft(x1));
plot(SC_f,[abs(freq_DSBSC),abs(freq_x1)]);
title('DSBSC, SNR=0db in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');
subplot(3,1,2);
freq_x2=fftshift(fft(x2));
plot(SC_f,[abs(freq_DSBSC),abs(freq_x2)]);
title('DSBSC, SNR=10db in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');
subplot(3,1,3);
freq_x3=fftshift(fft(x3));
plot(SC_f,[abs(freq_DSBSC),abs(freq_x3)]);
title('DSBSC, SNR=30db in frequency domain');
xlabel('Frequency');
ylabel('Amplitude');

soundsc(real(double(x1)));
pause(1);
clear sound;
soundsc(real(double(x2)));
pause(1);
clear sound;
soundsc(real(double(x3)));
pause(1);
clear sound;

%--------------------------------------------------------------------------
%9))
%frequency error
fc=100100;
a1 = y1.*cos(2*pi*fc*re_t);
a2 = y2.*cos(2*pi*fc*re_t);
a3 = y3.*cos(2*pi*fc*re_t);
figure;
[b,a] = butter(5,fc/(fs/2));
%filtfilt remove the high frequency noise and has no phase delay
a1 = filtfilt(b,a,a1);

%frequency error plotting in time domain
subplot(3,1,1);
plot(re_t,[x1,a1]);
title('SNR=0db in time domain with frequency error');
xlabel('Time');
ylabel('Amplitude');
[b,a] = butter(5,fc/(fs/2));
a2 = filtfilt(b,a,a2);
subplot(3,1,2);
plot(re_t,[x2,a2]);
title('SNR=10db in time domain with frequency error');
xlabel('Time');
ylabel('Amplitude');
subplot(3,1,3);
[b,a] = butter(5,fc/(fs/2));
a3 = filtfilt(b,a,a3);
plot(re_t,[x3,a3]);
title('SNR=30 db in time domain with frequency error');
xlabel('Time');
ylabel('Amplitude');

%frequency error plotting in time domain
figure;
subplot(3,1,1);
freq_a1=fftshift(fft(a1));
plot(SC_f,[abs(freq_x1),abs(freq_a1)]);
title('SNR=0db in frequency domain with frequency error');
xlabel('Frequency');
ylabel('Amplitude');
subplot(3,1,2);
freq_a2=fftshift(fft(x2));
plot(SC_f,[abs(freq_x2),abs(freq_a2)]);
title('SNR=10db in frequency domain with frequency error');
xlabel('Frequency');
ylabel('Amplitude');
subplot(3,1,3);
freq_a3=fftshift(fft(a3));
plot(SC_f,[abs(freq_x3),abs(freq_a3)]);
title('SNR=30db in frequency domain with frequency error');
xlabel('Frequency');
ylabel('Amplitude');

soundsc(real(double(a1)));
pause(1);
clear sound;
soundsc(real(double(a2)));
pause(1);
clear sound;
soundsc(real(double(a3)));
pause(1);
clear sound;

%The name of this phenomenon is Doppler effect.

%--------------------------------------------------------------------------
%10))
%phase error

fc=100000;
%phase error plotting in time domain
figure;
a1 = y1.*cos(2*pi*fc*re_t+20);
[b,a] = butter(5,fc/(fs/2));
a1 = filtfilt(b,a,a1);
subplot(3,1,1);
plot(re_t,[x1,a1]);
title('SNR=0db in time domain with phase error');
xlabel('Time');
ylabel('Amplitude');
a2 = y2.*cos(2*pi*fc*re_t+20);
[b,a] = butter(5,fc/(fs/2));
a2 = filtfilt(b,a,a2);
subplot(3,1,2);
plot(re_t,[x2,a2]);
title('SNR = 10 db in time domain with phase error');
xlabel('Time');
ylabel('Amplitude');
a3 = y3.*cos(2*pi*fc*re_t+20);
[b,a] = butter(5,fc/(fs/2));
a3 = filtfilt(b,a,a3);
subplot(3,1,3);
plot(re_t,[x3,a3]);
title('SNR=30 db in time domain with phase error');
xlabel('Time');
ylabel('Amplitude');

%phase error plotting in time domain
figure;
subplot(3,1,1);
freq_a1=fftshift(fft(a1));
plot(SC_f,[abs(freq_x1),abs(freq_a1)]);
title('SNR=0 db in frequency domain with phase error');
xlabel('Frequency');
ylabel('Amplitude');
subplot(3,1,2);
freq_a2=fftshift(fft(x2));
plot(SC_f,[abs(freq_x2),abs(freq_a2)]);
title('SNR=10 db in frequency domain with phase error');
xlabel('Frequency');
ylabel('Amplitude');
subplot(3,1,3);
freq_a3=fftshift(fft(a3));
plot(SC_f,[abs(freq_x3),abs(freq_a3)]);
title('SNR=30db in frequency domain with phase error');
xlabel('Frequency');
ylabel('Amplitude');

soundsc(real(double(a1)));
pause(1);
clear sound;
soundsc(real(double(a2)));
pause(1);
clear sound;
soundsc(real(double(a3)));
pause(1);
clear sound;


