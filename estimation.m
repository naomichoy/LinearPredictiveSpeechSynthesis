path='C:\Users\cat97\Documents\UniofSurrey\sem1\EEEM030-speech\assignment1\speech-samples\speech\';
file=strcat(path,'heed_m.wav')

[x, Fs] = audioread(file);
num_samples = Fs * 0.1;      % 100ms samples
% num_samples = Fs * 0.06;      % pure tone samples
start = 0.05*Fs;
cropped_x = x(start:start+num_samples);
size_vector = size(cropped_x);
L = size_vector(1);
t = (0:length(cropped_x) - 1) / Fs;     % for plotting 

%% Noise filtering
% calculate autocorrolation
[autocorr_values,lags] = xcorr(cropped_x, 'coeff'); % Calculate normalized autocorrelation
figure;
plot(lags/Fs,autocorr_values);
xlabel('Lag');
ylabel('Autocorrelation');

% shortest dist peaks
[pksh,lcsh] = findpeaks(autocorr_values);
short = mean(diff(lcsh))/Fs    % noise  freq

% longest dist peaks (max lag)
noise_peak_thres = 0.25;
[pklg,lclg] = findpeaks(autocorr_values,'MinPeakDistance',ceil(short),'MinPeakheight',noise_peak_thres);
long = max(diff(lclg))/Fs   % for fundemental freq

% check values functions
% max(diff(lcsh))/Fs
% max(diff(lclg))/Fs

%% TODO calculate fundemental freq

hold on
pks = plot(lags(lcsh)/Fs,pksh,'or', lags(lclg)/Fs,pklg+0.05,'vk');
hold off
legend(pks,[repmat('Period: ',[2 1]) num2str([short;long],0)])


% Calculate the power spectral density (PSD) of the autocorrelation
N = length(autocorr_values);
frequencies = (0:N-1) * Fs / N;
psd = abs(fft(autocorr_values));

% psd graph for checking 
% plot(frequencies, psd);
% xlabel('Frequency (Hz)');
% ylabel('Power Spectral Density');
% title('PSD of Autocorrelation');

noiseFrequencies = frequencies(psd > noise_peak_thres);
% for i = 1:10
    % high pass filter
    d = designfilt('highpassiir', 'FilterOrder', 10, 'StopbandFrequency', noiseFrequencies(1), 'SampleRate', Fs);
    filtered_x = filter(d, cropped_x);
% end

% signal plot
figure;
plot(t, cropped_x, 'b-');
hold on;
plot(t, filtered_x,'r-');
hold off;
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Sound Waveform');
legend('original cropped', 'filtered')

%% plot amplitude specturm
N = length(filtered_x);
frequencies = (0:N-1) * Fs / N;
amplitude_spectrum = abs(fft(filtered_x));
amplitude_spectrum_db = 20 * log10(amplitude_spectrum);

figure;
plot(frequencies, amplitude_spectrum_db);
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
title('Amplitude Spectrum of Filtered Signal');


% % spectrogram
% segmentlen = 100;
% noverlap = 90;
% NFFT = 128;
% spectrogram(x,segmentlen,noverlap,NFFT,Fs,'yaxis');
% title('Signal Spectrogram')


% lpc
lpc_coefficients = lpc(filtered_x, 10);
[H, w] = freqz(1, lpc_coefficients, 1024);

figure;
plot(w/pi, 20*log10(abs(H)));
xlabel('Normalized Frequency (\pi radians/sample)');
ylabel('Amplitude (dB)');
title('LPC Filter Frequency Response');

% first 3 formant frequencies

