clear all;

% load lpc coefficient
path='C:\Users\cat97\Documents\UniofSurrey\sem1\EEEM030-speech\assignment1\';
fname = 'heed_m_lpc20_1200.mat';
fullpath=strcat(path,fname)
load(fullpath,'lpc_coefficients');

% extract num_samples from filename
num_samples = split(fname,"_");
num_samples = split(num_samples(4),".");
num_samples = str2num(num_samples{1});

% fundemental frequency and period
f1 = 94.1176;
Fs = 24000;
f1_p = Fs/f1

% generate excitation signal
pulse_train = zeros(1, Fs); % Define the duration
pulse_train(1:f1_p:end) = 1;

% add all-pole filter to pulse train with lpc coefficients
synthesized_signal = filter(1, lpc_coefficients, pulse_train);
% synthesized_signal = conv(pulse_train, lpc_coefficients);
% pwr = 1.0193
% synthesized_signal = synthesized_signal/ max(abs(synthesized_signal)) * pwr;

% plot output signal
figure(1);
t = (0:length(synthesized_signal) - 1) / Fs;
plot(t, synthesized_signal, 'b-');
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Sound Waveform');

N = length(synthesized_signal);
frequencies = (0:N-1) * Fs / N;
amplitude_spectrum = abs(fft(synthesized_signal));
amplitude_spectrum_db = 20 * log10(amplitude_spectrum);

figure(3);
plot(frequencies, amplitude_spectrum_db);
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
title('Amplitude Spectrum of Synthesised Signal');


% [H,W]=freqz(1,[1 -A],512)


% save output
sound(synthesized_signal, Fs);
outfile = split(fname,".");
outfile = strcat(outfile{1}, ".wav");
audiowrite(outfile, synthesized_signal, Fs);