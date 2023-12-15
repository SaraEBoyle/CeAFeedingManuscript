%% Translates neuron indeces from cell registration to session specific indeces,
% Then checks responsive neurons across sessions to compare which neurons
% are active when

%% TODO
% 1. take the overlapping neuron indices from the cell reg analysis of fr
% and sated high fat and sucrose sessions. The deleted neurons need to be
% added back to adjust the indices from the cell reg in order to check
% whether the neurons are responsive to the food.

% Choose folder with stored cell reg results
clear
results_folder = uigetdir(pwd);
stim_length = 10; % check for responsiveness 5 seconds after food contact
rerun_permutations = 0; %If you've already run this before for your data, switch to 0
% Find log file and read source data folders
F = '*logFile*.txt';
S = dir(fullfile(results_folder,F));
N = S.name;
log_file = [results_folder, '\', N];
fid = fopen(log_file);
tline = fgetl(fid);
disp(tline)
in_session_list = 0;
session_folders = {};
while ischar(tline)
    tline = fgetl(fid);
    
    if or(strcmp('General data parameters:', tline), strcmp('', tline))
        % If the line is right after sessions list
        disp([num2str(size(session_folders, 2)), ' sessions read']);
        break
    end
    
    if strcmp(tline, '-----------------')
        % if the line is the one right before sessions list
        in_session_list = 1; % indicate you've reach sessions
        tline = fgetl(fid);
    end
    
    if in_session_list == 1
        disp(tline)
        session_folders{size(session_folders, 2) + 1} = tline(13:(end - 22));
    end 
end
fclose(fid);

% Find cell registered file to get the cell indeces for further analysis
F = '*cellRegistered*.mat';
S = dir(fullfile(results_folder,F));
if size(S, 2) > 1
    disp('This only can handle folders with 1 registration result. Pick your favorite and try again.')
end
N = S.name;
load([results_folder, '\', N]);
cell_indeces = cell_registered_struct.cell_to_index_map;
disp('Cell indeces loaded.');
% Find the deleted neuron indeces for each session. Responsive neurons are
% calculated before removal of bad neurons, so you will need to add in bad
% neurons to get accurate responsiveness flags for each session.

deleted_neurons = {};
first_excited = {};
first_inhibited = {};
second_excited = {};
second_inhibited = {};
dFFs = {};
ZScores = {};
trial_lengths = [];
trial_indss = {};
relative_sorted_food_timess = {};

for sess_ind = 1:size(session_folders, 2)
    data_folder = session_folders{sess_ind};
    load([data_folder, 'z_scores_session_1.mat']);
    deleted_neurons{sess_ind} = neurons_to_delete;
    NeurondFF(neurons_to_delete, :) = [];
    NeuronZScores_baseline_normed(neurons_to_delete, :) = [];
    ZScores{sess_ind} = NeuronZScores_baseline_normed;
    dFFs{sess_ind} = NeurondFF; 
    if rerun_permutations
        load([data_folder, 'responsive_neurons_session_1.mat']);
        load([data_folder, 'time_of_food_contact.mat']);
        trial_lengths = [trial_lengths, trial_length];
        trial_indss{sess_ind} = trial_inds;
        relative_sorted_food_timess{sess_ind} = relative_sorted_food_times;

        first_excited{sess_ind} = first_excited_neurons;
        first_inhibited{sess_ind} = first_inhibited_neurons;
        second_excited{sess_ind} = second_excited_neurons;
        second_inhibited{sess_ind} = second_inhibited_neurons;
    end
end












if rerun_permutations
first_excited = {};
first_inhibited = {};
second_excited = {};
second_inhibited = {};
third_excited = {};
third_inhibited = {};

for sess_ind = [1, 2]
    %% Run a permutation test of the dFFs of each neuron again
    framerate = 10;

    trial_inds = trial_indss{sess_ind};
    relative_sorted_food_times = relative_sorted_food_timess{sess_ind};
    labels = food_labels;
    %if an airpuff session
    TrialTypes = trial_inds;
    trial_n = [1, 2, 3];

    baseline = [3 38]; %3 seconds right before the stimulus
    stimulus = [60 200]; %100-125 is the CS-- tone. 126-150,160 is the US-- water/shock
    NUM_ITER = 184756; %This is 2*trials choose trials. This is for 20 trials
    N = 100;

    %FrameLength = size(DF,2)/nTrials;

    %% calculate auROC %%

    %% mean baseline per neuron in first liquid
    mean_first_baseline = {};
    mean_first_responses = {};
    mean_second_baseline = {};
    mean_second_responses = {};
    mean_third_baseline = {};
    mean_third_responses = {};
    sess_dFFs = dFFs{sess_ind};
    first_dFFs = sess_dFFs(:, 1);
    second_dFFs = sess_dFFs(:, 2);
    third_dFFs = sess_dFFs(:, 3);
    for o=1:size(first_dFFs,1)   % neuron number. Calculates average baseline of each neuron    
          %First calculate first z scores. trial num may be different so need
          %separate loops

          %Need second loop since delivery times are different trial by trial
          high_trials = find(trial_inds == 1);
          cur_neuron = first_dFFs{o,:}; %pulls DF/F for each frame for 1 neuron
          temp_baselines = [];
          temp_responses = [];
          %for each high fat trial
          for tri = 1:size(cur_neuron, 1)
              stimulus(1) = round(relative_sorted_food_times(high_trials(tri), 1) * framerate);
              stimulus(2) = stimulus(1) + stim_length*framerate;
              if stimulus(2) > size(cur_neuron, 2)
                  stimulus(2) = size(cur_neuron, 2);
              end
              temp_baselines = [temp_baselines;mean(cur_neuron(tri,baseline(1):baseline(2)))];     % calcualte average vuale within baseline period
              temp_responses = [temp_responses;mean(cur_neuron(tri,stimulus(1):stimulus(2)))];     % calcualte average vuale within baseline period
          end
          mean_first_baseline{o} = mean(temp_baselines,2);
          mean_first_responses{o} = mean(temp_responses,2); 
    end

    for o=1:size(second_dFFs,1)   % neuron number. Calculates average baseline of each neuron    
          %First calculate second z scores. trial num may be different so need
          %separate loops

          %Need second loop since delivery times are different trial by trial
          low_trials = find(trial_inds == 2);
          cur_neuron = second_dFFs{o,:}; %pulls DF/F for each frame for 1 neuron
          temp_baselines = [];
          temp_responses = [];
          %for each high fat trial
          for tri = 1:size(cur_neuron, 1)
              stimulus(1) = round(relative_sorted_food_times(low_trials(tri), 1) * framerate);
              stimulus(2) = stimulus(1) + stim_length*framerate;
              if stimulus(2) > size(cur_neuron, 2)
                  stimulus(2) = size(cur_neuron, 2);
              end
              temp_baselines = [temp_baselines;mean(cur_neuron(tri,baseline(1):baseline(2)))];     % calcualte average vuale within baseline period
              temp_responses = [temp_responses;mean(cur_neuron(tri,stimulus(1):stimulus(2)))];     % calcualte average vuale within baseline period
          end
          mean_second_baseline{o} = mean(temp_baselines,2);
          mean_second_responses{o} = mean(temp_responses,2); 
    end

    for o=1:size(third_dFFs,1)   % neuron number. Calculates average baseline of each neuron    
          %First calculate third z scores. trial num may be different so need
          %separate loops

          %Need second loop since delivery times are different trial by trial
          miss_trials = find(trial_inds == 3);
          cur_neuron = third_dFFs{o,:}; %pulls DF/F for each frame for 1 neuron
          temp_baselines = [];
          temp_responses = [];
          %for each high fat trial
          for tri = 1:size(cur_neuron, 1)
              stimulus(2) = stimulus(1) + stim_length*framerate;
              if stimulus(2) > size(cur_neuron, 2)
                  stimulus(2) = size(cur_neuron, 2);
              end
              stimulus(1) = round(relative_sorted_food_times(miss_trials(tri), 1) * framerate);
              temp_baselines = [temp_baselines;mean(cur_neuron(tri,baseline(1):baseline(2)))];     % calcualte average vuale within baseline period
              temp_responses = [temp_responses;mean(cur_neuron(tri,stimulus(1):stimulus(2)))];     % calcualte average vuale within baseline period
          end
          mean_third_baseline{o} = mean(temp_baselines,2);
          mean_third_responses{o} = mean(temp_responses,2); 
    end

    %% permutation test first liquid TODO
    % For each neuron (20 baseline 20 responses), calculate difference between
    % response and baseline. Find average difference. Next, shuffle labels and
    % calculate mean difference. Repeat 10,000 times, then see how many
    % permutations are less than the actual value, the percent is p. Might be
    % better to do area under the curve
    p_values_first = [];
    excitations = [];
    for n = 1:size(mean_first_responses, 2)
        actual_baselines = mean_first_baseline{n};
        actual_responses = mean_first_responses{n};
        pos = 0 < mean(actual_responses - actual_baselines);
        actual_dif = abs(mean(actual_responses - actual_baselines));


        to_shuffle = vertcat(actual_baselines, actual_responses);
        shuffled_difs = zeros(1, 100000);
        for g = 1:100000
            order = randperm(size(to_shuffle, 1));
            rand_baselines = to_shuffle(order(1:size(actual_baselines, 1)));
            rand_responses = to_shuffle(order((1 + size(actual_baselines, 1):end)));
            shuffled_dif = abs(mean(rand_responses - rand_baselines));
            shuffled_difs(g) = shuffled_dif;
        end

        biggers = find(shuffled_difs > actual_dif); %random difs bigger than actual
        p_value = length(biggers)/100000;
        p_values_first = horzcat(p_values_first, p_value);
        excitations = horzcat(excitations, pos);
    end

    first_responding_neurons = find(p_values_first < .05);

    disp('High fat responding neurons:');
    disp((horzcat(num2str(length(first_responding_neurons)), ' out of ', num2str(n))))
    disp(num2str(first_responding_neurons));

    first_excited_neurons = [];
    first_inhibited_neurons = [];
    for resp = first_responding_neurons
        if excitations(resp) == 1
            first_excited_neurons = horzcat(first_excited_neurons, resp);
        else
            first_inhibited_neurons = horzcat(first_inhibited_neurons, resp);
        end
    end

    disp(horzcat('Number of excited neurons: ', num2str(size(first_excited_neurons, 2))));
    disp(horzcat('Number of inhibited neurons: ', num2str(size(first_inhibited_neurons, 2))));

    p_values_second = [];
    excitations = [];
    for n = 1:size(mean_second_responses, 2)
        actual_baselines = mean_second_baseline{n};
        actual_responses = mean_second_responses{n};
        actual_dif = abs(mean(actual_responses - actual_baselines));
        pos = 0 < mean(actual_responses - actual_baselines);
        excitations = horzcat(excitations, pos);

        to_shuffle = vertcat(actual_baselines, actual_responses);
        shuffled_difs = zeros(1, 100000);
        for g = 1:100000
            order = randperm(size(to_shuffle, 1));
            rand_baselines = to_shuffle(order(1:size(actual_baselines, 1)));
            rand_responses = to_shuffle(order((1 + size(actual_baselines, 1):end)));
            shuffled_dif = abs(mean(rand_responses - rand_baselines));
            shuffled_difs(g) = shuffled_dif;
        end

        biggers = find(shuffled_difs > actual_dif); %random difs bigger than actual
        p_value = length(biggers)/100000;
        p_values_second = horzcat(p_values_second, p_value);
    end

    second_responding_neurons = find(p_values_second < .05);

    %delete the bad quality neurons you flagged earlier

    disp('Low fat responding neurons:');
    disp((horzcat(num2str(length(second_responding_neurons)), ' out of ', num2str(n))))
    disp(num2str(second_responding_neurons));

    second_excited_neurons = [];
    second_inhibited_neurons = [];
    for resp = second_responding_neurons
        if excitations(resp) == 1
            second_excited_neurons = horzcat(second_excited_neurons, resp);
        else
            second_inhibited_neurons = horzcat(second_inhibited_neurons, resp);
        end
    end

    disp(horzcat('Number of excited neurons: ', num2str(size(second_excited_neurons, 2))));
    disp(horzcat('Number of inhibited neurons: ', num2str(size(second_inhibited_neurons, 2))));



    p_values_third = [];
    excitations = [];
    for n = 1:size(mean_third_responses, 2)
        actual_baselines = mean_third_baseline{n};
        actual_responses = mean_third_responses{n};
        actual_dif = abs(mean(actual_responses - actual_baselines));
        pos = 0 < mean(actual_responses - actual_baselines);
        excitations = horzcat(excitations, pos);

        to_shuffle = vertcat(actual_baselines, actual_responses);
        shuffled_difs = zeros(1, 100000);
        for g = 1:100000
            order = randperm(size(to_shuffle, 1));
            rand_baselines = to_shuffle(order(1:size(actual_baselines, 1)));
            rand_responses = to_shuffle(order((1 + size(actual_baselines, 1):end)));
            shuffled_dif = abs(mean(rand_responses - rand_baselines));
            shuffled_difs(g) = shuffled_dif;
        end

        biggers = find(shuffled_difs > actual_dif); %random difs bigger than actual
        p_value = length(biggers)/100000;
        p_values_third = horzcat(p_values_third, p_value);
    end

    third_responding_neurons = find(p_values_third < .05);

    disp(horzcat('Missed trial responding neurons:'));
    disp((horzcat(num2str(length(third_responding_neurons)), ' out of ', num2str(n))))
    disp(num2str(third_responding_neurons));

    third_excited_neurons = [];
    third_inhibited_neurons = [];
    for resp = third_responding_neurons
        if excitations(resp) == 1
            third_excited_neurons = horzcat(third_excited_neurons, resp);
        else
            third_inhibited_neurons = horzcat(third_inhibited_neurons, resp);
        end
    end

    disp(horzcat('Number of excited neurons: ', num2str(size(third_excited_neurons, 2))));
    disp(horzcat('Number of inhibited neurons: ', num2str(size(third_inhibited_neurons, 2))));
    first_excited{sess_ind} = first_excited_neurons;
    first_inhibited{sess_ind} = first_inhibited_neurons;
    second_excited{sess_ind} = second_excited_neurons;
    second_inhibited{sess_ind} = second_inhibited_neurons;
    third_excited{sess_ind} = third_excited_neurons;
    third_inhibited{sess_ind} = second_inhibited_neurons;
end

save([results_folder, '\responsive_neurons_after_deletions.mat'], 'first_excited', 'first_inhibited', 'second_excited', 'second_inhibited', 'third_excited', 'third_inhibited')
else
    %load([results_folder, '\responsive_neurons_after_deletions.mat']);
end


cells_in_both = [];

for ind = 1:size(cell_indeces, 1)
    if cell_indeces(ind, 1)*cell_indeces(ind, 2) ~= 0
        cells_in_both = [cells_in_both, ind];
    end
end

disp(['Found ', num2str(size(cells_in_both, 2)), ' neurons overlapping']);

% Compare responses of overlapping cells
indeces_in_both_sessions = cell_indeces(cells_in_both, :);
probability_overlap = cell_registered_struct.cell_scores(cells_in_both);
save([results_folder, '\responses_in_tracked_cells.mat'], 'indeces_in_both_sessions','probability_overlap');

if ~exist('first_excited_neurons', 'var')
    disp('Baseline only session finished');
    return
end

responses_first = zeros(2, size(cells_in_both, 2)); % array n x 2. 0 is no response, 1 is positive response, -1 is inhibitory respons
responses_second = zeros(2, size(cells_in_both, 2)); % top row fat, bottom row sucrose
for c = 1:size(cells_in_both, 2)
    dual_ind = cells_in_both(c);
    disp(['Neuron ', num2str(dual_ind), ':']);
    current_neuron_inds = cell_indeces(dual_ind, :);
    
    % check first session responses
    ses_1_id = current_neuron_inds(1);
    disp(['Session 1 ID: ', num2str(ses_1_id)]);
    
    first_excited_neurons = first_excited{1};
    if ismember(ses_1_id, first_excited_neurons)
        disp('Excited by high fat in session 1');
        responses_first(1, c) = 1;
    end
    first_inhibited_neurons = first_inhibited{1};
    if ismember(ses_1_id, first_inhibited_neurons)
        disp('Inhibited by high fat in session 1');
        responses_first(1, c) = -1;
    end
    
    second_excited_neurons = second_excited{1};
    if ismember(ses_1_id, second_excited_neurons)
        disp('Excited by sucrose in session 1');
        responses_first(2, c) = 1;
    end
    second_inhibited_neurons = second_inhibited{1};
    if ismember(ses_1_id, second_inhibited_neurons)
        disp('Inhibited by sucrose in session 1');
        responses_first(2, c) = -1;
    end
    
    %check second session responses
    ses_2_id = current_neuron_inds(2);
    disp(['Session 2 ID: ', num2str(ses_2_id)]);
    
    first_excited_neurons = first_excited{2};
    if ismember(ses_2_id, first_excited_neurons)
        disp('Excited by high fat in session 2');
        responses_second(1, c) = 1;
    end
    first_inhibited_neurons = first_inhibited{2};
    if ismember(ses_2_id, first_inhibited_neurons)
        disp('Inhibited by high fat in session 2');
        responses_second(1, c) = -1;
    end
    
    second_excited_neurons = second_excited{2};
    if ismember(ses_2_id, second_excited_neurons)
        disp('Excited by sucrose in session 2');
        responses_second(2, c) = 1;
    end
    second_inhibited_neurons = second_inhibited{2};
    if ismember(ses_2_id, second_inhibited_neurons)
        disp('Inhibited by sucrose in session 2');
        responses_second(2, c) = -1;
    end
    disp(newline);
end
save([results_folder, '\responses_in_tracked_cells.mat'], 'indeces_in_both_sessions','probability_overlap', 'responses_first', 'responses_second');



%% Next check baseline activity in these cells
% if zero_negative_baseline =  1, delete negative values
zero_negative_baseline = 0;
find_baseline_activity(session_folders{1, 1}, 50, zero_negative_baseline);
load([session_folders{1, 1}, 'baseline_analysis.mat']);

tracked_cells_first_inds = indeces_in_both_sessions(:, 1);
first_session_baselines = total_session_baseline_AUCs(tracked_cells_first_inds);

find_baseline_activity(session_folders{1, 2}, 16, zero_negative_baseline);
load([session_folders{1, 2}, 'baseline_analysis.mat']);
tracked_cells_second_inds = indeces_in_both_sessions(:, 2);
second_session_baselines = total_session_baseline_AUCs(tracked_cells_second_inds);

save('baselines_tracked.mat', 'first_session_baselines', 'second_session_baselines');







































%% Here I'm trying to make up for the fact I didn't remove deleted neurons here.
%Not sure why but it doesn't work at all. Re run permutation test after
%removing bad neurons, then try on that data


%{
% Convert the cell reg indeces to what they'd be before removal so that
% they can be compared to the responsive neuron indeces

for sess_ind = 1:size(deleted_neurons, 2)
    % for each deleted neuron, adjust the index of all registered neurons
    % Ex: cell reg indeces are 2, 5, 8, 9, 10. deleted 3 previously, so to 
    % adjust, imagine 3 was added back. This means all indeces bigger or 
    %equal to 3 have to increase by 1.
    current_session_inds = cell_indeces(:, sess_ind);
    to_add_1 = zeros(size(current_session_inds, 1), 1);
    for del_neuron = deleted_neurons{sess_ind}
        to_add_1_inds = find(current_session_inds >= del_neuron);
        preserve_zeros = find(current_session_inds == 0);
        
        to_add_1(to_add_1_inds) = to_add_1(to_add_1_inds) + 1;
        to_add_1(preserve_zeros) = 0;
        
    end
    current_session_inds = current_session_inds + to_add_1;
    cell_indeces(:, sess_ind) = current_session_inds;
end

% Now find overlapping cells from each session

cells_in_both = [];

for ind = 1:size(cell_indeces, 1)
    if cell_indeces(ind, 1)*cell_indeces(ind, 2) ~= 0
        cells_in_both = [cells_in_both, ind];
    end
end

disp(['Found ', num2str(size(cells_in_both, 2)), ' neurons overlapping']);

% Compare responses of overlapping cells
for c = cells_in_both
    disp(['Neuron ', num2str(c), ':']);
    current_neuron_inds = cell_indeces(c, :);
    
    % check first session responses
    ses_1_id = current_neuron_inds(1);
    disp(['Session 1 ID: ', num2str(ses_1_id)]);
    
    first_excited_neurons = first_excited{1};
    if ismember(ses_1_id, first_excited_neurons)
        disp('Excited by high fat in session 1');
    end
    first_inhibited_neurons = first_inhibited{1};
    if ismember(ses_1_id, first_inhibited_neurons)
        disp('Inhibited by high fat in session 1');
    end
    
    second_excited_neurons = second_excited{1};
    if ismember(ses_1_id, second_excited_neurons)
        disp('Excited by sucrose in session 1');
    end
    second_inhibited_neurons = second_inhibited{1};
    if ismember(ses_1_id, second_inhibited_neurons)
        disp('Inhibited by sucrose in session 1');
    end
    
    %check second session responses
    ses_2_id = current_neuron_inds(2);
    disp(['Session 2 ID: ', num2str(ses_2_id)]);
    
    first_excited_neurons = first_excited{2};
    if ismember(ses_2_id, first_excited_neurons)
        disp('Excited by high fat in session 2');
    end
    first_inhibited_neurons = first_inhibited{2};
    if ismember(ses_2_id, first_inhibited_neurons)
        disp('Inhibited by high fat in session 2');
    end
    
    second_excited_neurons = second_excited{2};
    if ismember(ses_2_id, second_excited_neurons)
        disp('Excited by sucrose in session 2');
    end
    second_inhibited_neurons = second_inhibited{2};
    if ismember(ses_2_id, second_inhibited_neurons)
        disp('Inhibited by sucrose in session 2');
    end
    disp(newline);
end
%}




