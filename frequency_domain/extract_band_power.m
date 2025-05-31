function [object, times] = extract_band_power(object, no_objects, srate, tlimits, freq_min, winsize)
% EXTRACT_FREQ_BAND Extracts a frequency band from each object and stores it in the structure
%
%   [object, times] = extract_freq_band(object, no_objects, srate, tlimits, freq_min, freq_max)
%
%   This function performs time-frequency analysis using EEGLAB's newtimef function
%   for each object in a structure. It extracts power in a specified frequency band
%   and stores the result back in the object structure.
%
%   Inputs:
%       object      - structure with fields 'ob1', 'ob2', ..., 'obN', each containing erp data [channels x time x trials]
%       no_objects  - number of objects (fields) in the structure
%       srate       - sampling rate in Hz (e.g., 512)
%       tlimits     - time limits in ms (e.g., [-100 1000])
%       freq_min    - minimum frequency for the band (e.g., 4)
%       freq_max    - maximum frequency for the band (e.g., 6)
%       winsize     - Window size (in samples) for the Morlet wavelet convolution.
%                     Typically set to 300 samples. It must be proportional to your
%                     erp time window. bigger winsize generally means worse
%                     time resolution but better frequency resolution 
%
%   Output:
%       object      - updated structure with field 'ersp_band' added to each object
%       times       - time points corresponding to ERSP data (in ms)

eeglab;

for i = 1:no_objects

    % Access the object field dynamically and permute
    data = permute(object.(['ob' num2str(i)]), [2 3 1]);  % from [trials x channels x time] to [channels x time x trials]
    frames = size(data, 2);   % Assuming all objects have same trial size
    n_channels = size(data, 1);
    ersp_ch = [];
    ersp_ch_reverted = [];
    object.(['ob' num2str(i)]) = [];
    for ch = 1:n_channels
    
     [ersp,itc,powbase,times,freqs,erspboot,itcboot,tfdata]  = newtimef(data(ch,:,:), frames, tlimits, srate, [3 0.5], ...
                 'freqs', [4 40], 'nfreqs', 100, 'plotersp', 'off', 'baseline', 0, 'plotitc', 'off', 'winsize', winsize);
    
            freq_idx = freqs >= freq_min & freqs <= freq_max;
            ersp_z = ersp(freq_idx, :);  % [~2 freqs x times]
            ersp_ch(:,:,ch) = ersp_z;    % [freq(power_at_that_freq) x time x channels]
    end
ersp_ch_reverted = permute(ersp_ch, [1 3 2]); % [freq(power_at_that_freq) x channels x time]
object.(['ob' num2str(i)]) = ersp_ch_reverted;

end
end