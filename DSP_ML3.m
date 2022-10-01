close all; clear;clc; load('projIB.mat');
set(0, 'defaultFigureUnits', 'centimeters', 'defaultFigurePosition', [0 0 20 30]);

wp = 2500*2/fs;
ws = 4000*2/fs;
rp = 3;
rs = 95;
IIR_gain = 10^(40/20);
FIR_gain = 10^(38/20);
%% Butterworth

close all;
[n, Wn] = buttord(wp,ws,rp,rs);
fprintf("Order of Butterworth filter: %d\n", n);

% Convert to SOS
[z, p, k] = butter(n,Wn);
[SOS, G] = zp2sos(z,p,k);
hd = dfilt.df2sos(SOS,G);
fprintf("Number of multiplication operations: %d per sample\n", 2*n + 1);

hd2 = dfilt.scalar(G);
final_Butt = cascade(hd2,hd); % Apply amplification


filtered = filter(final_Butt,noisy);
soundsc(filtered,fs);
pause(4);
getInfo(final_Butt,'Butterworth Filter')

%% Chebyshev Type I

close all;
[n, Wn] = cheb1ord(wp,ws,rp,rs);
fprintf("Order of Chebyshev Type I filter: %d\n", n)

[z, p, k] = cheby1(n,rp,wp);
[SOS, G] = zp2sos(z,p,k);
hd = dfilt.df2sos(SOS,G);

fprintf("Number of multiplication operations: %d per sample\n", 2*n + 1);

hd2 = dfilt.scalar(IIR_gain);
final_Cheby1 = cascade(hd2,hd); % Apply amplification


filtered = filter(final_Cheby1,noisy);
soundsc(filtered,fs);
pause(4);
getInfo(final_Cheby1,'Chebyshev Type I Filter');
%% Chebyshev Type II

close all;
[n, Wn] = cheb2ord(wp,ws,rp,rs);
fprintf("Order of Chebyshev Type II filter: %d\n", n)

[z, p, k] = cheby2(n,rs,ws);
[SOS, G] = zp2sos(z,p,k);
hd = dfilt.df2sos(SOS,G);

fprintf("Number of multiplication operations: %d per sample\n", 2*n + 1);

hd2 = dfilt.scalar(IIR_gain);
final_Cheby2 = cascade(hd2,hd); % Apply amplification

filtered = filter(final_Cheby2,noisy);
soundsc(filtered,fs);
pause(4);
getInfo(final_Cheby2,'Chebyshev Type II Filter');
%% Elliptic

close all;
[n, Wn] = ellipord(wp,ws,rp,rs);
fprintf("Order of Elliptic filter: %d\n", n)

% Convert to SOS
[z, p, k] = ellip(n,rp,rs,wp);
[SOS, G] = zp2sos(z,p,k);
hd = dfilt.df2sos(SOS, G);

fprintf("Number of multiplication operations: %d per sample\n", 2*n + 1);

hd2 = dfilt.scalar(IIR_gain);
final_Ellip = cascade(hd2, hd); % Apply amplification


filtered = filter(final_Ellip,noisy);
soundsc(filtered,fs);
pause(4);
getInfo(final_Ellip,'Elliptic Filter');
%% Parks-McClellan

close all;
[n, Fo, Ao, W] = firpmord([2500 4000],[1 0], [(10^(rp/40)-1)/(10^(rp/40)+1)  10^(-rs/20)], fs);
fprintf("Order of Parks McLellan filter: %d\n", n)

% Convert to DF1
b = firpm(n,Fo,Ao,W);
hd = dfilt.df1(b);

fprintf("Number of multiplication operations: %d per sample\n", n + 1);

hd2 = dfilt.scalar(FIR_gain);
final_parkmc = cascade(hd2, hd); % Apply amplification

filtered = filter(final_parkmc,noisy);
soundsc(filtered,fs);
pause(4);
getInfo(final_parkmc,'Parks McLellan Filter');
%% Kaiser

close all;
[n, Wn, beta, ftype] = kaiserord([2500 4000], [1 0], [(10^(rp/40)-1)/(10^(rp/40)+1)  10^(-rs/20)], fs);
fprintf("Order of Kaiser filter: %d\n", n)

fprintf("Number of multiplication operations: %d per sample\n", n + 1);

final_Kaiser = designfilt('lowpassfir','PassbandFrequency',wp, ...
         'StopbandFrequency',ws,'PassbandRipple',rp/2, ...
         'StopbandAttenuation',rs,'DesignMethod','kaiserwin');
     
filtered = filter(final_Kaiser, noisy);
soundsc(filtered,fs);
getInfo(final_Kaiser,'Kaiser Filter');
%%
function getInfo(filt, Title)
    [H,W] = freqz(filt);
    [g,w] = grpdelay(filt);
    
    figure;
    subplot(7,1,1);
    plot(W,20*log10(abs(H)));
    title("Frequency Response of " + Title);
    xlim([0 pi]);
    xticks([0 pi/4 pi/2 3*pi/4 pi]);
    xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
    xlabel('Normalized Frequency');
    ylabel('Magnitude(dB)');
    
    
    subplot(7,1,2);
    plot(W,abs(H));
    title("Frequency Response of " + Title);
    xlim([0 pi]);
    xticks([0 pi/4 pi/2 3*pi/4 pi]);
    xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
    xlabel('Normalized Frequency');
    ylabel('Magnitude');

    subplot(7,1,3);
    plot(w,g);
    title("Group Delay of " + Title);
    xlim([0 pi]);
    xticks([0 pi/4 pi/2 3*pi/4 pi]);
    xticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
    xlabel('Normalized Frequency');
    ylabel('Group Delay');
    
    subplot(7,1,4:5);
    [z,p] = filt.zpk;
    zplane(z,p);
    title("Pole-Zero Plot of  " + Title);
    
    subplot(7,1,6:7);
    imp=zeros(1,100);
    imp(1)=1;
    impresponse = filter(filt,imp);
    samps = 0:99;
    stem(samps, impresponse);
    xlabel('Samples');
    ylabel('Amplitude');
    title("Impulse Response of "+Title);
end

%% Comments
% In this project I put a sound sample with alot of noise through various
% LPF's to get a coherent message. Each filter was cascaded and the
% resulting filters each did a satisfactory job of filtering out the noise.
% Furthermore, the order of each filter was calculated along with the
% number of multiplication operations per input sample required to
% implement the filter. Also, the mag responce (in dB and focused on
% passband ripple) were plot as well as the group delay. Also, the
% pole/zero plot and impulse responses (of first 100 samples) were plot as
% well.