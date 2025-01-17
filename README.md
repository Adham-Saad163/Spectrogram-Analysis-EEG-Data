# EEG Spectrogram Analysis for Epileptic Seizure Detection Using STFT

## Overview
This project focuses on analyzing EEG data to detect epileptic seizures using spectrograms generated via the Short-Time Fourier Transform (STFT). It explores the impact of various parameters (window types, lengths, and overlap ratios) on the time-frequency representation of EEG signals.

## Features
- **EEG Signal Analysis:** Utilized the CHB-MIT dataset to distinguish between ictal (seizure) and non-ictal states.
- **Spectrogram Generation:** Implemented custom STFT methods and evaluated them against MATLAB's built-in functions.
- **Parameter Optimization:** Investigated effects of window types, lengths, and overlap ratios on spectrogram accuracy.

## Dataset
The EEG data is sourced from the [CHB-MIT Scalp EEG Database](https://physionet.org/content/chbmit/1.0.0/), specifically focusing on patient `chb12`.  
- **File Structure:** Raw EEG signals with seizure and non-seizure intervals annotated.  
- **Usage:** Preprocessed and divided into ictal and non-ictal states for analysis.  

## Results

### Window Type Comparison
- **Hamming Window:** Best balance between time and frequency resolution.  
- **Blackman Window:** Minimal spectral leakage but lower frequency resolution.  
- **Rectangular Window:** Significant spectral leakage, not recommended.  

### Overlap Ratio
- **85% Overlap:** Smooth transitions and high time resolution.  
- **50% Overlap:** Balanced approach with reduced computational cost.  
- **10% Overlap:** Poor time resolution with blocky artifacts.  

### Window Length
- **Short (128 samples):** High time resolution, poor frequency resolution.  
- **Medium (512 samples):** Balanced resolution for EEG analysis.  
- **Long (2560 samples):** High frequency resolution, poor time resolution.  

## Contributors
The following individuals have contributed to this project:  

- **Adham Ahmed**  
  Email: [s-adham.saad@zewailcity.edu.eg](mailto:s-adham.saad@zewailcity.edu.eg)  

- **Aya Sherif**  
  Email: [s-aya.nassef@zewailcity.edu.eg](mailto:s-aya.nassef@zewailcity.edu.eg)  

