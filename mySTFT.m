function [S, f, t] = mySTFT(signal, fs, winLen, overlap, window)
% mySTFT: Computes the Short-Time Fourier Transform (STFT) and spectrogram.
% 
% Inputs:
%   - signal: Input time-domain signal (1D vector).
%   - fs: Sampling frequency (Hz).
%   - winLen: Length of each analysis window (samples).
%   - overlap: Overlap between adjacent windows (samples).
%   - window: Window function vector (length should match winLen).
%
% Outputs:
%   - S: Magnitude spectrogram (matrix).
%   - f: Frequency vector (Hz).
%   - t: Time vector (s).

% Ensure signal is a column vector
signal = signal(:);

% Ensure the window length matches the specified winLen
if length(window) ~= winLen
    error('Window length does not match winLen. Please provide a valid window.');
end

% Define step size
step = winLen - overlap;

% Calculate the number of frames
numFrames = floor((length(signal) - overlap) / step);

% Frequency vector
f = (0:winLen/2) * (fs / winLen);

% Time vector (center of each window)
t = (0:numFrames-1) * (step / fs) + (winLen / (2 * fs));

% Preallocate spectrogram matrix
S = zeros(length(f), numFrames);

% Compute STFT
for n = 1:numFrames
    % Define start and end indices for the current window
    startIdx = (n-1) * step + 1;
    endIdx = startIdx + winLen - 1;
    
    % Check for zero-padding at the end of the signal
    if endIdx > length(signal)
        segment = signal(startIdx:end); % Truncate signal
        segment = [segment; zeros(winLen - length(segment), 1)]; % Zero-pad
    else
        segment = signal(startIdx:endIdx);
    end
    
    % Apply window function
    segment = segment .* window;
    
    % Compute FFT and retain only positive frequencies
    fftResult = fft(segment);
    S(:, n) = abs(fftResult(1:winLen/2+1));
end

end
