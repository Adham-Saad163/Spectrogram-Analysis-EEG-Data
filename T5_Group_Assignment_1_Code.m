clc;
clear;
close all;
%______________________________________________ I. Pre-processing _____________________________________________________

% 1. Import and Read the EEG Records

% Load seizure-free (non-ictal) and seizure (ictal) EEG records
NonIctalData=load('chb12_32_data.mat'); % Non-ictal (seizure-free) data
data_struct = load('chb12_29_data.mat'); % Ictal data (including both ictal and non-ictal)
raw_data_ictal = data_struct.data; % Ictal data (used in chb12_29_data.mat)
raw_data_non_ictal = NonIctalData.data; % Non-ictal data (from chb12_32_data.mat)

% Seizure intervals based on summary file (in seconds) for ictal data
seizure_intervals = [107, 146; 554, 592; 1163, 1199; 1401, 1447; 1884, 1921; 3557, 3584]; 

% Define empty channels for chb12_29.edf based on the summary file
empty_channels = [4, 10, 14, 20]; % These channels have no data

% Initialize ictal_data and non_ictal_data
fs = 256; % Sampling frequency
ictal_data = []; 
non_ictal_data = [];
num_channels = 29;

% Exclude empty channels from the ictal data
valid_channels = setdiff(1:size(raw_data_ictal, 2), empty_channels); % Indices of valid channels
filtered_data_ictal = raw_data_ictal(:, valid_channels); % Only include valid channels

% Calculate the average signal across valid channels for ictal data
average_raw_data_ictal = mean(filtered_data_ictal, 2, 'omitnan'); % Exclude NaN values during mean calculation

% Calculate the average signal across valid channels for non-ictal data
filtered_data_non_ictal = raw_data_non_ictal(:, valid_channels); % Only include valid channels
average_raw_data_non_ictal = mean(filtered_data_non_ictal, 2, 'omitnan'); % Exclude NaN values during mean calculation

% Plot both ictal and non-ictal data in the same figure
figure;
plot((1:length(average_raw_data_ictal)) / fs, average_raw_data_ictal, 'b', 'LineWidth', 1.5); % Plot ictal data
hold on; % Keep the current plot
plot((1:length(average_raw_data_non_ictal)) / fs, average_raw_data_non_ictal, 'g', 'LineWidth', 1.5); % Plot non-ictal data

% Highlight seizure intervals
for i = 1:size(seizure_intervals, 1)
    x_patch = [seizure_intervals(i, 1), seizure_intervals(i, 2), seizure_intervals(i, 2), seizure_intervals(i, 1)];
    y_patch = [min(average_raw_data_ictal), min(average_raw_data_ictal), max(average_raw_data_ictal), max(average_raw_data_ictal)];
    patch(x_patch, y_patch, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Add red highlighted region
end
title('Ictal (Blue) vs Non-Ictal (Green) Average Signal');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Ictal Data', 'Non-Ictal Data'); % Remove 'Seizure Intervals' from the legend
grid on;
hold off;



% Iterate over seizure intervals to extract seizure data
for i = 1:size(seizure_intervals, 1)
    start_sample = seizure_intervals(i, 1) * fs;
    end_sample = seizure_intervals(i, 2) * fs;
    if end_sample > size(raw_data_ictal, 1) % Check for bounds
        error('Seizure interval %d exceeds the data size.', i);
    end
   ictal_data = [ictal_data; raw_data_ictal(start_sample:end_sample, valid_channels)];
 % Concatenate ictal data
end

%________________________ Start of Plot Ictal States Stacaked Back To Back _____________________________________

% Extract seizure (ictal) data and stack back-to-back
ictal_stacked = [];
ictal_end_times = []; % To store cumulative end times of each seizure interval
cumulative_time = 0; % Initialize cumulative time

for i = 1:size(seizure_intervals, 1)
    % Get start and end sample indices
    start_sample = seizure_intervals(i, 1) * fs;
    end_sample = seizure_intervals(i, 2) * fs;
    
    % Extract ictal segment for valid channels
    ictal_segment = filtered_data_ictal(start_sample:end_sample, :);
    
    % Average across all valid channels
    ictal_average = mean(ictal_segment, 2, 'omitnan'); % Average across channels
    
    % Append data and update cumulative time
    ictal_stacked = [ictal_stacked; ictal_average];
    cumulative_time = cumulative_time + size(ictal_segment, 1) / fs; % Increment by the seizure duration
    ictal_end_times = [ictal_end_times; cumulative_time]; % Store the end time
end

% Generate time axis for the stacked data
time_axis = (0:length(ictal_stacked) - 1) / fs;

% Plot ictal (seizure) data stacked back-to-back
figure;
plot(time_axis, ictal_stacked, 'r', 'LineWidth', 1.5);
hold on;

% Add dashed lines at the end of each seizure interval
for i = 1:length(ictal_end_times)
    xline(ictal_end_times(i), '--k', sprintf('End of Ictal %d', i), ...
        'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top');
end

% Finalize plot
datacursormode on;
title('Seizure Data Stacked Back-to-Back');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Seizure Data');
grid on;
hold off;
%________________________ End of Plot Ictal States Stacaked Back To Back _____________________________________







%________________________ Start of Plot Ictal States Stacaked Back To Back Non icatal data_____________________________________


% Define empty channels for non-ictal data
empty_channels_non_ictal = [10, 13, 18, 23]; % These channels have no data

% Exclude empty channels from the non-ictal data
valid_channels_non_ictal = setdiff(1:size(raw_data_non_ictal, 2), empty_channels_non_ictal); % Valid channels
filtered_data_non_ictal = raw_data_non_ictal(:, valid_channels_non_ictal); % Filter valid channels

% Extract non-ictal data corresponding to ictal intervals
non_ictal_stacked = [];
non_ictal_end_times = []; % To store cumulative end times of each interval
cumulative_time_non_ictal = 0; % Initialize cumulative time

for i = 1:size(seizure_intervals, 1)
    % Get start and end sample indices for the ictal interval
    start_sample = seizure_intervals(i, 1) * fs;
    end_sample = seizure_intervals(i, 2) * fs;
    
    % Check for bounds
    if end_sample > size(raw_data_non_ictal, 1)
        error('Seizure interval %d exceeds the size of the non-ictal data.', i);
    end
    
    % Extract the corresponding non-ictal segment
    non_ictal_segment = filtered_data_non_ictal(start_sample:end_sample, :);
    
    % Average across all valid channels
    non_ictal_average = mean(non_ictal_segment, 2, 'omitnan'); % Average across channels
    
    % Append data and update cumulative time
    non_ictal_stacked = [non_ictal_stacked; non_ictal_average];
    cumulative_time_non_ictal = cumulative_time_non_ictal + size(non_ictal_segment, 1) / fs; % Increment by the interval duration
    non_ictal_end_times = [non_ictal_end_times; cumulative_time_non_ictal]; % Store the end time
end

% Generate time axis for the stacked non-ictal data
time_axis_non_ictal = (0:length(non_ictal_stacked) - 1) / fs;

% Plot non-ictal data stacked back-to-back
figure;
plot(time_axis_non_ictal, non_ictal_stacked, 'b', 'LineWidth', 1.5);
hold on;

% Add dashed lines at the end of each interval
for i = 1:length(non_ictal_end_times)
    xline(non_ictal_end_times(i), '--k', sprintf('End of Non-Ictal %d', i), ...
        'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top');
end

% Finalize plot
datacursormode on;
title('Non-Ictal Data Corresponding to Ictal Intervals Stacked Back-to-Back');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Non-Ictal Data');
grid on;
hold off;





%________________________ End of Plot Ictal States Stacaked Back To Back Non icatal data_____________________________________

% % Display the size of ictal_data for verification
% disp('Ictal Data:');
% disp(size(ictal_data));



%________________________ Start of Channel Selection and Data Validation _____________________________________


% Initialize channel selection loop
while true
    % Prompt user to select channels
    user_input = input('Enter the channels to select (e.g., [1, 2, 3]): ', 's');
    
    % Check if the user wants to exit
    if strcmpi(user_input, 'exit')
        disp('Exiting the program...');
        return; % Stop execution
    end
    
    % Convert input to numeric array
    try
        channels = str2num(user_input); %#ok<ST2NM> (used intentionally to evaluate string input)
        
        % Check if input is valid
        if isempty(channels) || ~isnumeric(channels)
            error('Invalid input. Please enter numeric values in the format [1, 2, 3].');
        end
        
        % Check if channels are within range
        if any(channels < 1) || any(channels > num_channels)
            error(['Channel numbers must be between 1 and ', num2str(num_channels), '.']);
        end
        
        % Remove empty channels from the selected list
        valid_channels = setdiff(channels, empty_channels);
        omitted_channels = intersect(channels, empty_channels);
        
        % Display message for omitted channels, if any
        if ~isempty(omitted_channels)
            disp(['The following channels were omitted because they are empty: ', num2str(omitted_channels)]);
        end
        
        % Check if there are valid channels left
        if isempty(valid_channels)
            error('No valid channels remain after omitting empty channels. Please try again.');
        end
        
        % If valid, display selected channels and break loop
        disp(['You selected valid channels: ', num2str(valid_channels)]);
        break; % Exit loop after successful input
        
    catch ME
        % Display error message for invalid input
        disp(['Error: ', ME.message]);
        disp('Please try again.');
    end
end

%________________________ End of Channel Selection and Data Validation _____________________________________


% Process the selected channels (example: averaging)
disp('Processing selected channels...');
average_data = mean(raw_data_ictal(:, valid_channels), 2); % Compute the average of selected non-empty channels
disp('Processing complete.');
head(average_data);

% Check if average_data contains NaN
if any(isnan(average_data))
    disp('Warning: The processed data contains NaN values. Check the selected channels and input data.');
else
    disp('Data processed successfully. No NaN values present.');
end

% Plot 6 ictals raw signal segments corresponding to seizure intervals
figure;
for i = 1:size(seizure_intervals, 1)
    % Convert seizure times (seconds) to sample indices
    start_idx = round(seizure_intervals(i, 1) * fs);
    end_idx = round(seizure_intervals(i, 2) * fs);
    
    % Extract signal segment
    segment = average_data(start_idx:end_idx);
    
    % Plot each segment
    subplot(size(seizure_intervals, 1), 1, i); % Create subplot for each interval
    plot(linspace(seizure_intervals(i, 1), seizure_intervals(i, 2), length(segment)), segment);
    title(['Seizure Interval ' num2str(i) ' (' num2str(seizure_intervals(i, 1)) '-' num2str(seizure_intervals(i, 2)) 's)']);
    xlabel('Time (s)');
    ylabel('Amplitude');
end



% Plot Raw Signal with Seizure Intervals Highlighted
figure;
plot((1:length(average_data))/fs, average_data);
hold on;
for i = 1:size(seizure_intervals, 1)
    x_patch = [seizure_intervals(i, 1), seizure_intervals(i, 2), seizure_intervals(i, 2), seizure_intervals(i, 1)];
    y_patch = [min(average_data), min(average_data), max(average_data), max(average_data)];
    patch(x_patch, y_patch, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Highlight seizure region
end
hold off;
title('Raw Signal with Seizure Intervals Highlighted');
xlabel('Time (s)');
ylabel('Amplitude');



% Plot Raw Signal with Seizure Intervals Highlighted
figure;
plot((1:length(average_data))/fs, average_data);
hold on;
for i = 1:size(seizure_intervals, 1)
    x_patch = [seizure_intervals(i, 1), seizure_intervals(i, 2), seizure_intervals(i, 2), seizure_intervals(i, 1)];
    y_patch = [min(average_data), min(average_data), max(average_data), max(average_data)];
    patch(x_patch, y_patch, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % Highlight seizure region
end
hold off;
title('Raw Signal with Seizure Intervals Highlighted');
xlabel('Time (s)');
ylabel('Amplitude');




%______________________________________________ II. Spectrogram _____________________________________________________

% 5. Spectrogram Construction Using STFT
n = 2; % Window length in seconds
window_samples = n * fs;
overlap_ratio = 0.5; %  Overlap ratio
overlap_samples = floor(window_samples * overlap_ratio);
hop_size = window_samples - overlap_samples;

% Dictionary of window types
window_dict = struct( ...
    'blackman', @blackman, ...
    'hamming', @hamming, ...
    'triangular', @triang, ...
    'rectangular', @rectwin ...
);

% Select the desired window type
window_type = 'hamming'; % Change this to any desired window type: 'hamming', 'triangular', 'rectangular'
window_func = window_dict.(window_type); % Get the corresponding window function
% Generate the window as a vector
window = window_func(window_samples); % Apply window function to get the window vector

% Divide signal into overlapping frames
signal_length = length(average_data);
disp(['average data size: ', num2str(length(average_data))]);
num_frames = ceil((signal_length - overlap_samples) / hop_size);
disp(['Number of Frames: ', num2str(num_frames), ' with overlab ratio = ', num2str(overlap_ratio)]);


% Compute the STFT using your custom function
[S, F, T] = mySTFT(ictal_stacked, fs, window_samples, overlap_samples, window);

% Convert S to decibels for better visualization
S_dB = 20 * log10(abs(S));

% Plot the spectrogram
subplot(2, 2, 1); % Custome Spectrogram of Ictal-only 
imagesc(T, F, S_dB); % T: Time vector, F: Frequency vector, S_dB: Spectrogram in dB
axis xy; 
colormap parula; 
colorbar; 
xlabel('Time (s)');
ylabel('Frequency (Hz)');
% Add window, overlap ratio, and window length annotations
title_str = sprintf('Custom STFT Ictal\nWindow: %s, Length: %d samples, Overlap: %.1f%%', ...
                    window_type, window_samples, overlap_ratio * 100);
title(title_str); % Shorter title with window type and parameters
hold off;
datacursormode on; % Enable data cursor mode

% Plot spectrogram with built-in function for reference 
[S, F, T] = spectrogram(ictal_stacked, window, overlap_samples, window_samples, fs);
% spectrogram(signal, fs, winLen, overlap)
subplot(2, 2, 2); % Built-In Spectrogram of Ictal-only 
imagesc(T, F, 20*log10(abs(S)));
axis xy;
hold off;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
% Add window, overlap ratio, and window length annotations for built-in spectrogram
title_str = sprintf('Built-In Spectrogram Ictal\nWindow: %s, Length: %d samples, Overlap: %.1f%%', ...
                    window_type, window_samples, overlap_ratio * 100);
title(title_str); % Shorter title with window type and parameters
colormap parula;
colorbar;
datacursormode on; % Enable data cursor






%________________________ End of Spectrogram Plot (Built-In STFT + Custom STFT)_____________________________________


[S, F, T] = mySTFT(non_ictal_stacked, fs, window_samples, overlap_samples, window);
S_dB = 20 * log10(abs(S));
% Plot the spectrogram for non-ictal
subplot(2, 2, 3); % Custome Spectrogram of Ictal-only 
imagesc(T, F, S_dB); % T: Time vector, F: Frequency vector, S_dB: Spectrogram in dB
axis xy; 
colormap parula; 
colorbar; 
xlabel('Time (s)');
ylabel('Frequency (Hz)');
% Add window, overlap ratio, and window length annotations
title_str = sprintf('Custom STFT Non-Ictal \nWindow: %s, Length: %d samples, Overlap: %.1f%%', ...
                    window_type, window_samples, overlap_ratio * 100);
title(title_str); % Shorter title with window type and parameters
hold off;
datacursormode on; % Enable data cursor mode

% Plot spectrogram with built-in function for reference 
[S, F, T] = spectrogram(non_ictal_stacked, window, overlap_samples, window_samples, fs);
% spectrogram(signal, fs, winLen, overlap)
subplot(2, 2, 4); % Built-In Spectrogram of Ictal-only 
imagesc(T, F, 20*log10(abs(S)));
axis xy;
hold off;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
% Add window, overlap ratio, and window length annotations for built-in spectrogram
title_str = sprintf('Built-In Spectrogram Non-Ictal\nWindow: %s, Length: %d samples, Overlap: %.1f%%', ...
                    window_type, window_samples, overlap_ratio * 100);
title(title_str); % Shorter title with window type and parameters
colormap parula;
colorbar;
datacursormode on; % Enable data cursor
