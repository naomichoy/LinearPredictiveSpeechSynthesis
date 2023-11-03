path='C:\Users\cat97\Documents\UniofSurrey\sem1\EEEM030-speech\assignment1\speech-samples\speech\';
file=strcat(path,'heed_m.wav')

[x, Fs] = audioread(file);
sample_length = 0.1;    % in seconds 
num_samples = Fs * sample_length;     
% num_samples = Fs * 0.06;      % pure tone samples
start = 0.05*Fs;
cropped_x = x(start:start+num_samples);
size_vector = size(cropped_x);
L = size_vector(1);
t = (0:length(cropped_x) - 1) / Fs;     % for plotting 

%% Noise filtering
% calculate autocorrolation
[autocorr_values,lags] = xcorr(cropped_x, 'coeff'); % Calculate normalized autocorrelation
figure(1);
plot(lags/Fs,autocorr_values);
xlabel('Lag');
ylabel('Autocorrelation');

% shortest dist peaks
[pksh,lcsh] = findpeaks(autocorr_values);
short = mean(diff(lcsh))/Fs    % noise  freq

% longest dist peaks (max lag)
[pklg,lclg] = findpeaks(autocorr_values,'MinPeakDistance',ceil(short),'MinPeakheight',0.3);
long = max(diff(lclg))/Fs  % for fundemental freq

% check values functions
% max(diff(lcsh))/Fs
% max(diff(lclg))/Fs

%% calculate fundemental freq
f1 = (long / (2*pi)) * Fs

hold on
pks = plot(lags(lcsh)/Fs,pksh,'or', lags(lclg)/Fs,pklg+0.05,'vk');
hold off
legend(pks,[repmat('Period: ',[2 1]) num2str([short;long],0)])

% Calculate the power spectral density (PSD) of the autocorrelation
N = length(autocorr_values);
frequencies = (0:N-1) * Fs / N;
psd = abs(fft(autocorr_values));

% psd graph for checking 
figure(5);
plot(frequencies, psd);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density');
title('PSD of Autocorrelation');

% filter noise frequency
noise_peak_thres = 0.25;
noiseFrequencies = frequencies(psd < noise_peak_thres);
noiseFrequencies = noiseFrequencies(noiseFrequencies > 0);
noiseFrequencies = noiseFrequencies(noiseFrequencies <= 12000);
% for i = 1:10
    % high pass filter
    d = designfilt('highpassiir', 'FilterOrder', 10, 'StopbandFrequency', noiseFrequencies(1), 'SampleRate', Fs);
    filtered_x = filter(d, cropped_x);
% end

% signal plot
figure(2);
plot(t, cropped_x, 'b-');
hold on;
plot(t, filtered_x,'r-');
hold off;
xlabel('Time (seconds)');
ylabel('Amplitude');
title('Sound Waveform');
legend('original', 'filtered')

%% plot amplitude specturm
N = length(filtered_x);
frequencies = (0:N-1) * Fs / N;
amplitude_spectrum = abs(fft(filtered_x));
amplitude_spectrum_db = 20 * log10(amplitude_spectrum);

figure(3);
plot(frequencies, amplitude_spectrum_db);
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
title('Amplitude Spectrum of Filtered Signal');


%% lpc
lpc_coefficients = lpc(filtered_x, 20);
[H, w] = freqz(1, lpc_coefficients, 1024);

%% plot freq response
figure(4);
plot(w/pi, 20*log10(abs(H)));
xlabel('Normalized Frequency (\pi radians/sample)');
ylabel('Amplitude (dB)');
title('LPC Filter Frequency Response');

%% first 3 formant frequencies
poles = roots(lpc_coefficients);
mag = abs(poles);
pole_angle = angle(poles);
formant_freq_hz = (pole_angle / (2*pi)) * Fs;
[formant_freq_hz, ind] = sort(formant_freq_hz(formant_freq_hz > 0));
first_3_formant = formant_freq_hz([1:3])
