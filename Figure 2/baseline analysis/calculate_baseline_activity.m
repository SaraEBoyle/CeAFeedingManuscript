%% calculates baseline from dF/F collected in sessions with no stimuli delivered
load('z_scores_session_1.mat'); load('neurons_to_delete.mat');
% delete bad neurons off the bat

baseline = 2:38;
chop_food_responses = 1;

if chop_food_responses
    disp('looking at baseline in food delivery session')
    % delete food responses to look only at baseline
    NeurondFF(neurons_to_delete, :) = [];
    new_dFF = {};
    tri_count = 1;
    for n = 1:size(NeurondFF, 1)
        all_types = [];
        for c = 1:size(NeurondFF, 2)
            current = NeurondFF{n, c};
            if isempty(current)
                continue
            end
            current = current(:, baseline);
            all_types = [all_types; current];
        end
        new_dFF{tri_count} = all_types;
        tri_count = tri_count + 1;
    end
    NeurondFF = new_dFF';
else
    NeurondFF(neurons_to_delete) = [];
end


avg_single_neuron_dFF = []; % calculate mean dFF for individual neurons- least useful
percent_above_thresh = []; % calculate 2*std threshold for individual neuron, calculate time above threshold
avg_threshold_crossings = []; % calculate 2*std threshold for individual neuron, calculate positive threshold crossings
thresholds = [];

for n = 1:size(NeurondFF, 1)
    %% first get mean dFF
    current = NeurondFF{n};
    avg_single_neuron_dFF = [avg_single_neuron_dFF, mean(mean(current))];
    flat_neuron = reshape(current, 1, size(current, 1)*size(current, 2));
    
    %% Now calculate 2*std threshold for more resilient measures
    thresh = 2*std(flat_neuron);
    thresholds = [thresholds, thresh];
    
    frames_above_thresh = find(flat_neuron > thresh);
    percent_above_thresh = [percent_above_thresh, size(frames_above_thresh, 2)/size(flat_neuron, 2)];
    
    %% find number of times crossing above threshold
    prev = frames_above_thresh(1);
    crossings = 0;
    for frame = frames_above_thresh(2:end)
        if frame > (prev + 3)
            crossings = crossings + 1;
        end
        prev = frame;
    end
    avg_threshold_crossings = [avg_threshold_crossings, crossings/size(flat_neuron, 2)];
    
end

disp('done')

save('baseline_stats.mat', 'chop_food_responses', 'baseline', 'avg_threshold_crossings', 'percent_above_thresh', 'avg_single_neuron_dFF');






