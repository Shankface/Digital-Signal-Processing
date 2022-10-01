%% ECE310 Project 4 - Ayden Shankman

clc; clear; close all;

load s1; load s5; load vowels;

%% Question 1

u = 4e9;
fs = 5e6;
T = 2e-4;
t = 0:1/fs:T;
x = cos(2*pi*u*t.^2);

figure;
spectrogram(x, triang(256), 255, 256, fs, 'yaxis');
title("Spectrogram of Linear FM Chirp, \mu = 4.0x10^9");
%% Question 2

f1 = u*t;
phi = 2*pi*u*t.^2;
f2 = (1/(2*pi))*(diff(phi)./diff(t));

figure;
plot(t, f1, t(2:end), f2);
title('Instantaneous Frequency');
xlabel('Time (sec)');
ylabel('Frequency (Hz)');
xlim([0 max(t)]);
legend('f_1(t)', 'f_2(t)');
% f2 corresponds to the slope of the ridge in the spectrogram
%% Question 3

u2 = 1e10;
x2 = cos(2*pi*u2*t.^2);

figure;
spectrogram(x2, triang(256), 255, 256, fs, 'yaxis');
title("Spectrogram of Linear FM Chirp, \mu = 1.0x10^{10}");
% Increasing the chirp rate causes the slope to increase and
% the frequency exceeds the Nyquist Frequency which is 2.5MHz and causes
% aliasing
%% Question 4

fs = 8e3;
length_narrow = 1024;

figure;
subplot(2, 1, 1);
spectrogram(s1, triang(length_narrow), length_narrow-1, 4096, fs, 'yaxis');
title("Narrowband Spectrogram of s1");
subplot(2, 1, 2);
spectrogram(s5, triang(length_narrow), length_narrow-1, 4096, fs, 'yaxis');
title("Narrowband Spectrogram of s5");
%% Question 5

length_wide =64

figure;
subplot(2, 1, 1);
spectrogram(s1, triang(length_wide), length_wide/2, 1024, fs, 'yaxis');
title("Wideband Spectrogram of s1");
subplot(2, 1, 2);
spectrogram(s5, triang(length_wide), length_wide-1,4096, fs, 'yaxis');
title("Wideband Spectrogram of s5");
%% Question 6

vpad = padarray(vowels, 2^nextpow2(length(vowels)) - length(vowels), 0, 'post');
stft = spectrogram(vpad, rectwin(256), 128, 1024, fs, 'yaxis');
inv = invSTFT(stft, 1024);

figure;
subplot(3, 1, 1);
plot(vowels);
title('Vowels');
xlabel('n');
ylabel('Amplitude');

subplot(3, 1, 2);
plot(inv);
title('Inverse STFT of Vowels');
xlabel('n');
ylabel('Amplitude');
xlim([0 length(vowels)]);

subplot(3, 1, 3);
plot(vowels' - inv(1:length(vowels)));
title('Difference Between Vowels & Inverse STFT of Vowels');
xlabel('n');
ylabel('Amplitude');
xlim([0 length(vowels)]);
ylim([-1e-10 1e-10]);
%% Question 7

stftDown = downsample(stft', 2);
inv = invSTFT(stftDown', 1024);

figure;
subplot(2, 1, 1);
plot(vowels);
title('Vowels');
xlabel('n');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(inv);
title('Downsampled Vowels');
xlabel('n');
ylabel('Amplitude');
xlim([0 length(vowels)/2]);

soundsc(vowels, fs);
pause(1.5);
soundsc(inv, fs);
%%
function out = invSTFT(stft, n)
    stft = [stft(1:end - 1, :); flipud(stft)];
    col_len = size(stft, 2);
    out = zeros(1, 128*(col_len + 1));
    for i = 1:col_len
        index = 128*(i-1)+1;
        inv = ifft(stft(:,i), n, 'symmetric')';
        out(index : index + 255) = out(index : index + 255) + real(inv(1:256));
    end
    out(129:length(out) - 128) = out(129:length(out) - 128)./2;
end

%% Comments
% This project consisted of  Spectrogram Analysis and Applications. In
% question 1, a chirp signal was generated as well as a spectrogram for the
% chirp. In question 2, two different definitions of "instantaneous
% frequency" were calculated and it was determined that the second
% definition corresponded to the slope of the ridge in the spectrogram of
% the chirp. In question  3, mu was increased and its spectrogram was
% generated. From the spectrogram it can be seen that this increases the
% slope significantly  and casues the frequency to exceed the Nyquist
% Frequency. In question 4 and 5, the spectrogram parameters were altered
% to obtain narrowband and wideband spectrograms for s1 and s2. In question
% 6, a function was created to get the Inverse STFT of a signal. In
% question 7, that function was used to obtain a faster version of the
% 'vowels' speech signal with the same frequency content.