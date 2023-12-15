
mean_peaks = [];
for w = 1:size(z_US_trials, 1)
    mean_peaks = horzcat(mean_peaks, max(z_US_trials(w, :)));
end

oil_peaks = mean_peaks;
grain_peaks = mean_peaks;