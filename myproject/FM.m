clc; %clears all the text from the Command Window
close all; % closes all open MATLAB figure windows
clear; %delete data in Workspace
%1))
%Reading from audio file
[mt,Fs]=audioread('eric.wav');
Nsamps=length(mt);
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
% remove all freq above 4000
d = designfilt('lowpassfir', 'FilterOrder', 8000, 'CutoffFrequency',4000, 'SampleRate', Fs);
mtnew = filter(d, mt);
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
sound(real(mtnew),Fs);
pause(Nsamp2/Fs);
%--------------------------------------------------------------------------
%2))
%Generating NBFM
kf =0.00003*pi;
Fc=100000;
c_Fs=5*Fc;
%to resample mtnew from Fs to c_Fs which is 5*Fc
re_mt=resample(mtnew,c_Fs,Fs);
re_Nsamps=length(re_mt);
re_t = linspace(0,(re_Nsamps/c_Fs), re_Nsamps);%time of modulation
%NBFM=Ac*cos(2?fct)cos(2?kf?m(t)dt)? Ac*sin(2?fct)sin(2?kf?m(t)dt)=
%=cos(2?fct)?(2?kf?m(t)dt)sin(2?fct)
%let Q=(2?kf?m(t)dt)
Q= 2*pi*kf.*cumsum(re_mt)';  
NBFM= cos(2*Fc*pi*re_t)-(Q.*sin(2*Fc*pi*re_t)); 
F_NBFM = fftshift(fft(NBFM));
f_NBFM = linspace((-c_Fs/2)-Fc,c_Fs/2+Fc,length(NBFM)); 
figure;
plot(f_NBFM,abs(F_NBFM));
title('NBFM in Frequency Domain');
%--------------------------------------------------------------------------
%3))
%condition to generate NBFM that modulation index < 1
%--------------------------------------------------------------------------
%4))
%Demodulation
NBFM_diff = diff(NBFM) ;  %Convert FM to AM by differentiation 
env =abs(hilbert(NBFM_diff))- mean(abs(hilbert(NBFM_diff))); %remove effect of AM from envelop detector                                               
figure(4);
plot(env);                                 
title('NBFM after demodulatin');



