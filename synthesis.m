clear all;

% load lpc coefficient
path='C:\Users\cat97\Documents\UniofSurrey\sem1\EEEM030-speech\assignment1\';
fname = 'heed_m_lpc10_2400.mat';
fullpath=strcat(path,fname)
load(fullpath,'lpc_coefficients');

% extract num_samples from filename
num_samples = split(fname,"_");
num_samples = split(num_samples(4),".");
num_samples = str2num(num_samples{1});

% fundemental frequency and period
f1 = 41.0620;
Fs = 24000;
f1_p = Fs/f1

% generate excitation signal
pulse_train = zeros(1, num_samples*10); % Define the duration
pulse_train(1:f1_p:end) = 1;

% filter pulse train with lpc coefficients
synthesized_signal = filter(1, lpc_coefficients, pulse_train);
% pwr = 1.0193
% synthesized_signal = synthesized_signal/ max(abs(synthesized_signal)) * pwr;

% output
sound(synthesized_signal, Fs);
outfile = split(fname,".");
outfile = strcat(outfile{1}, ".wav");
audiowrite(outfile, synthesized_signal, Fs);