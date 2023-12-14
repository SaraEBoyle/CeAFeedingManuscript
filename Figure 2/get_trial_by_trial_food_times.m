%% First get the times of each food contact for each trial
clear variables
load('time_of_food_contact.mat');
load('z_scores_session_1.mat');
load('responsive_neurons_session_1.mat');
food_timing_mat = {};
seconds_after = 15;
first = 0;
second = 0;
missed = 0;
thresh = 3; %std multiplier for thresholding
for ind = 1:size(trial_inds, 1)
    if trial_inds(ind) == 1
        first = first + 1;
        food_timing_mat{first, 1} = relative_sorted_food_times(ind);
    elseif trial_inds(ind) == 2
        second = second + 1;
        food_timing_mat{second, 2} = relative_sorted_food_times(ind);
    elseif trial_inds(ind) == 3
        missed = missed + 1;
        food_timing_mat{missed, 3} = relative_sorted_food_times(ind);
    end
end
disp('done')
save('food_timing_mat.mat', 'food_timing_mat');

%% Next get time of food contact for finding response magnitude
% Use time of pellet dispensed to get baseline time, calculate 2 SD
% threshold from that (first 4 seconds), then find if signal goes above
% that for rest  of trial

%fat = NeuronZScores(:, 1);
fat = NeurondFF(:, 1);
fat_times = food_timing_mat(:,1);
percent_time_above_per_neuron = [];
longest_response_per_neuron = [];
avg_Z_per_neuron = [];

for neur_ind = 1:size(fat, 1)
    % for every neuron
    cur_neuron = cell2mat(fat(neur_ind));
    
    baseline = cur_neuron(:, 1:40);
    activity_thresh = mean(mean(baseline)) + thresh*std(reshape(baseline, 1, size(baseline, 1)*size(baseline, 2)));
    above_thresh = cur_neuron > activity_thresh;
    
    responses = {};
    aboves = [];
    sizes = [];
    longest_responses = [];
    mean_zs = [];
    %for every trial
    for trial_ind = 1:size(cur_neuron, 1)
        cur_delivery = floor(fat_times{trial_ind});
        cur_above_thresh = above_thresh(trial_ind, :);
        cur_trial_z = cur_neuron(trial_ind, :);
        to_delete = find(cur_trial_z < 0);
        %cur_trial_z(to_delete) = 0;
        %cur_avg_z = .1* trapz(cur_trial_z(:, cur_delivery*10:(cur_delivery*10+ seconds_after*10)));
        cur_avg_z = mean(cur_trial_z(:, cur_delivery*10:(cur_delivery*10+ seconds_after*10)));
        mean_zs = [mean_zs, cur_avg_z];
        responses = cur_above_thresh(:, cur_delivery*10:(cur_delivery*10+ seconds_after*10));
        split_os = split(num2str(responses), '0  0');
        [max_size, max_index] = max(cellfun('size', split_os, 2));
        longest_response = strfind(split_os{max_index}, '1');
        longest_responses = [longest_responses, length(longest_response)];
        % add 1 divide by 3 to get number of 1s
        aboves = [aboves, sum(cur_above_thresh(:, cur_delivery*10:(cur_delivery*10 + seconds_after*10)))];
        sizes = [sizes, size(responses, 1)*size(responses, 2)];
    end
    avg_Z_per_neuron = [avg_Z_per_neuron, mean(mean_zs)];
    percent_time_above = sum(aboves)/sum(sizes);
    longest_response_per_neuron = [longest_response_per_neuron, mean(longest_responses)];
    percent_time_above_per_neuron = [percent_time_above_per_neuron, percent_time_above];
end
save('percent_time_above_per_neuron.mat', 'percent_time_above_per_neuron', 'longest_response_per_neuron');


