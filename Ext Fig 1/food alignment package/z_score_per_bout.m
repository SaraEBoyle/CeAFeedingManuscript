%% find the bout length and the peak and average z scores
% load z_scores_whole_session.mat
z_scores_whole_session;

%find the times of the events
%load('Alignment Data\white chocolate START alignment data.mat');
load('Alignment Data\white chocolate alignment data.mat');
start_times = trial_start_times + 20;
%load('Alignment Data\white chocolate STOP alignment data.mat');
load('Alignment Data\wc drops alignment data.mat');
end_times = trial_start_times + 20;
bout_lengths = end_times - start_times;
% find the frame these events happen
times = z_scores_whole_session(:, 1);
signal = z_scores_whole_session(:, 2);
avg_bout_z = [];
peak_bout_z = [];
for s_ind = 1:size(start_times, 2)
    bigs = find(times >= start_times(s_ind));
    smalls = find(times <= end_times(s_ind));
    common_inds = intersect(bigs, smalls);
    bout = mean(signal(common_inds));
    peak = max(signal(common_inds));
    avg_bout_z = horzcat(avg_bout_z, bout);
    peak_bout_z = horzcat(peak_bout_z, bout);
end

save('bout_info.mat', 'bout_lengths', 'avg_bout_z', 'peak_bout_z');