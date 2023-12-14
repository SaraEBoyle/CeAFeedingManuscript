%% define parameter %%
clear
%% Calculates means baselines of all neurons, then compares those means to the 
% stimulus evoked response in a permutation test to check for significance
%tem = dir('*Liquids');

P = pwd;
S = dir(fullfile(P,'*.mat'));
N = {S.name};
X = ~cellfun('isempty',strfind(N,'Deliver'));

names = N(X);
session_ind = 1; %which session to analyze

%Sort by time of session
session_times = [];
for p = 1:size(names, 2)
    time_string = names{p};
    time = time_string((end - 9):(end - 4));
    session_times = [session_times;str2double(time)];
end
sorted_session_times = sort(session_times);
sorted_names = {};
for t = 1:size(sorted_session_times, 1)
    order = find(session_times == sorted_session_times(t, 1));
    sorted_names{t} = names{order};
end
names = sorted_names;

%Load bpod file
load(fullfile(P,names{session_ind}));


%FileName   = 'Sorted.mat';
FolderName = pwd;
%File       = fullfile(FolderName, FileName);
%load(File);
%FileName   = 'US_SST6_amph_water_response_pool.mat';
%FolderName = '/Users/sboyle/Documents/Li Lab/Image analysis/SST5/amph';
%File       = fullfile(FolderName, FileName);
%load(File);

%a = load('US_SST6_amph_water_response_pool.mat');
%responsive_neurons = a.Test1_response_pool;
%unresponsive_neurons = 1:size(DF_new',2);
%unresponsive_neurons(:, responsive_neurons) = []; %take out responsive neurons

%% load all sessions' z scores, consolidate neurons to delete for consistency
%load extracted z score info

X = ~cellfun('isempty',strfind(N,'z_scores_session'));
z_score_names = N(X);

all_neurons_to_delete = [];
for z = 1:size(z_score_names, 2)
    load(z_score_names{z});
    all_neurons_to_delete = [all_neurons_to_delete, neurons_to_delete];
end
all_neurons_to_delete = unique(all_neurons_to_delete);

load(horzcat('z_scores_session_', num2str(session_ind), '.mat'));
framerate = 10;
%nTrials = 66;%20 for water, 15 for shock

%% Get trial info
%trial_length = 14;
%frames = size(neuron.C,2);
%trial_IDs = zeros(frames, 1);
%frame_rate = 10;
%frame_length = trial_length * frame_rate;
%trial_num = frames/frame_length;
%for p = 1:trial_num
%    trial_IDs(((p - 1)*frame_length + 1):(p*frame_length)) = p;
%end

if isfield(SessionData, 'TrialTypes')
    %if a liquid delivery session
    TrialTypes = SessionData.TrialTypes;
    trial_n = unique(SessionData.TrialTypes);

    if ismember(-1, trial_n)
        pallet_cleansers = 1; % if 1, you included trials with no imaging
    else
        pallet_cleansers = 0;
    end

    bad_inds = find(TrialTypes == -1);
    TrialTypes(bad_inds) = [];
    first_liquid = SessionData.TrialSettings.first_liquid;
    second_liquid = SessionData.TrialSettings.second_liquid;
    third_liquid = SessionData.TrialSettings.third_liquid;
    fourth_liquid = SessionData.TrialSettings.fourth_liquid;
    fifth_liquid = SessionData.TrialSettings.fifth_liquid;
    sixth_liquid = SessionData.TrialSettings.sixth_liquid;
else
    %if an airpuff session
    TrialTypes = ones(SessionData.nTrials, 1);
    trial_n = 1;
    pallet_cleansers = 0;
end


baseline = [5 39]; %3 seconds right before the stimulus
stimulus = [40 100]; %100-125 is the CS-- tone. 126-150,160 is the US-- water/shock
NUM_ITER = 184756; %This is 2*trials choose trials. This is for 20 trials
N = 100;

%FrameLength = size(DF,2)/nTrials;

%% calculate auROC %%
%{
DF_temp = DF'; %Same thing but transposed
for j = 1:size(DF_temp,2)  % cell number
    % Run auROC
    temp = DF_temp(:,j); %Pulls fluorescence from the cell we're on
    temp2 = reshape(temp,FrameLength,nTrials); %splits the array into ntrials * 200 frame chunks
    temp3 = temp2(:,event_trials);%Chooses the trials you're interested in
    temp4 = temp3(baseline(1):baseline(2),:); %These are baselines of each frame of each trial
    negtive = median(temp4,1);  %calculate median baseline value of each column ***
    for  k = 1:FrameLength %for each frame of each trial
         positive = temp3(k,:); %whole trial
         
        x = negtive(:);  % make sure that the data is stored in a column vector 
        y = positive(:);
        % define the parameter
        lenx = length(x); %20
        leny = length(y);
        threshlo = min([x; y]);
        threshhi = max([x; y]);
        thresh = linspace(threshlo,threshhi,N);
        fa = zeros(1,N);	% allocate the false alarm vector
        hit = zeros(1,N);
        
        for i = 1:N
            fa(N-i+1) = sum(x > thresh(i));
            hit(N-i+1) = sum(y > thresh(i));
        end
        fa = fa/leny;
        hit = hit/lenx;
        fa(1) = 0;
        hit(1) = 0;
        fa(N) = 1; 
        hit(N) = 1;
        auc = trapz(fa,hit); 
        
        auROC_datat(j,k) = auc;
    end
   
end

[B I]= sortrows(mean(auROC_datat(:,stimulus(1):stimulus(2)),2));
ranking_data = auROC_datat(I,:);
ranking_data = flipud(ranking_data);
%figure, imagesc(ranking_data,[0 1]);
%}


if ~isfield(SessionData, 'TrialTypes')
    %for aipuff sessions
    stimulus = [40 80]; %4 seconds
    air_z_scores = NeuronZScores_aligned(:, 1 + pallet_cleansers);
    for o=1:size(air_z_scores,1)   % neuron number. Calculates average baseline of each neuron    
          %First calculate first z scores. trial num may be different so need
          %separate loops
          cur_neuron = air_z_scores{o,:}; %pulls DF/F for each frame for 1 neuron
          temp_baselines = cur_neuron(:,baseline(1):baseline(2));     % calcualte average vuale within baseline period
          mean_air_baseline{o} = mean(temp_baselines,2);
          temp_responses = cur_neuron(:,stimulus(1):stimulus(2));     % calcualte average vuale within baseline period
          mean_air_responses{o} = mean(temp_responses,2); 
    end
    
    p_values_air = [];
    excitations = [];
    %do a permutation test for air responses
    for n = 1:size(mean_air_responses, 2)
        actual_baselines = mean_air_baseline{n};
        actual_responses = mean_air_responses{n};
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
        p_values_air = horzcat(p_values_air, p_value);
        excitations = horzcat(excitations, pos);
    end

    air_responding_neurons = find(p_values_air < .05);
    delete = intersect(air_responding_neurons, all_neurons_to_delete);
    delete_inds = [];
    for d = delete
        delete_ind = find(air_responding_neurons == d);
        delete_inds = [delete_inds;delete_ind];
    end
    air_responding_neurons(delete_inds) = [];

    disp(horzcat('Air responding neurons:'));
    disp((horzcat(num2str(length(air_responding_neurons)), ' out of ', num2str(n - length(all_neurons_to_delete)))))
    disp(num2str(air_responding_neurons));

    air_excited_neurons = [];
    air_inhibited_neurons = [];
    for resp = air_responding_neurons
        if excitations(resp) == 1
            air_excited_neurons = horzcat(air_excited_neurons, resp);
        else
            air_inhibited_neurons = horzcat(air_inhibited_neurons, resp);
        end
    end

    disp(horzcat('Number of excited neurons: ', num2str(size(air_excited_neurons, 2))));
    disp(horzcat('Number of inhibited neurons: ', num2str(size(air_inhibited_neurons, 2))));
    save(horzcat('responsive_neurons_session_', num2str(session_ind), '.mat'), 'air_responding_neurons', ...
        'air_excited_neurons', 'air_inhibited_neurons');
    return
end

%% mean baseline per neuron in first liquid
mean_first_baseline = {};
mean_first_responses = {};
mean_second_baseline = {};
mean_second_responses = {};
mean_third_baseline = {};
mean_third_responses = {};
mean_fourth_baseline = {};
mean_fourth_responses = {};
mean_fifth_baseline = {};
mean_fifth_responses = {};
mean_sixth_baseline = {};
mean_sixth_responses = {};


normed_z = NeuronZScores;
for tri_type = 1:size(normed_z, 2)
    %for each trial type
    neuron_block = normed_z(:, tri_type);
    for neur = 1:size(neuron_block, 1)
        %for each neuron
        single_neur = neuron_block{neur, :};
        if isempty(single_neur)
            continue
        end
        new_neuron = [];
        for tri = 1:size(single_neur, 1)
            single_neur_trial = single_neur(tri, :);
            base = mean(single_neur_trial(3:(baseline(2) - 1)));
            new = single_neur_trial - base;
            new_neuron(tri, :) = new;
        end
        normed_z{neur, tri_type} = new_neuron;
    end
end

NeuronZScores = NeuronZScores_aligned;
first_z_scores = NeuronZScores(:, 1 + pallet_cleansers);
second_z_scores = NeuronZScores(:, 2 + pallet_cleansers);
third_z_scores = NeuronZScores(:, 3 + pallet_cleansers);
fourth_z_scores = NeuronZScores(:, 4 + pallet_cleansers);
fifth_z_scores = NeuronZScores(:, 5 + pallet_cleansers);
sixth_z_scores = NeuronZScores(:, 6 + pallet_cleansers);
air_z_scores = NeuronZScores(:, 7 + pallet_cleansers);

for o=1:size(first_z_scores,1)   % neuron number. Calculates average baseline of each neuron    
      %First calculate first z scores. trial num may be different so need
      %separate loops
      cur_neuron = first_z_scores{o,:}; %pulls DF/F for each frame for 1 neuron
      temp_baselines = cur_neuron(:,baseline(1):baseline(2));     % calcualte average vuale within baseline period
      mean_first_baseline{o} = mean(temp_baselines,2);
      temp_responses = cur_neuron(:,stimulus(1):stimulus(2));     % calcualte average vuale within baseline period
      mean_first_responses{o} = mean(temp_responses,2); 
end

for o=1:size(second_z_scores,1)   % neuron number. Calculates average baseline of each neuron    
      %First calculate second z scores. trial num may be different so need
      %separate loops
      cur_neuron = second_z_scores{o,:}; %pulls DF/F for each frame for 1 neuron
      temp_baselines = cur_neuron(:,baseline(1):baseline(2));     % calcualte average vuale within baseline period
      mean_second_baseline{o} = mean(temp_baselines,2);
      temp_responses = cur_neuron(:,stimulus(1):stimulus(2));     % calcualte average vuale within baseline period
      mean_second_responses{o} = mean(temp_responses,2); 
end

for o=1:size(third_z_scores,1)   % neuron number. Calculates average baseline of each neuron    
      %First calculate third z scores. trial num may be different so need
      %separate loops
      cur_neuron = third_z_scores{o,:}; %pulls DF/F for each frame for 1 neuron
      temp_baselines = cur_neuron(:,baseline(1):baseline(2));     % calcualte average vuale within baseline period
      mean_third_baseline{o} = mean(temp_baselines,2);
      temp_responses = cur_neuron(:,stimulus(1):stimulus(2));     % calcualte average vuale within baseline period
      mean_third_responses{o} = mean(temp_responses,2); 
end

for o=1:size(fourth_z_scores,1)   % neuron number. Calculates average baseline of each neuron    
      %First calculate third z scores. trial num may be different so need
      %separate loops
      cur_neuron = fourth_z_scores{o,:}; %pulls DF/F for each frame for 1 neuron
      temp_baselines = cur_neuron(:,baseline(1):baseline(2));     % calcualte average vuale within baseline period
      mean_fourth_baseline{o} = mean(temp_baselines,2);
      temp_responses = cur_neuron(:,stimulus(1):stimulus(2));     % calcualte average vuale within baseline period
      mean_fourth_responses{o} = mean(temp_responses,2); 
end

for o=1:size(fifth_z_scores,1)   % neuron number. Calculates average baseline of each neuron    
      %First calculate third z scores. trial num may be different so need
      %separate loops
      cur_neuron = fifth_z_scores{o,:}; %pulls DF/F for each frame for 1 neuron
      temp_baselines = cur_neuron(:,baseline(1):baseline(2));     % calcualte average vuale within baseline period
      mean_fifth_baseline{o} = mean(temp_baselines,2);
      temp_responses = cur_neuron(:,stimulus(1):stimulus(2));     % calcualte average vuale within baseline period
      mean_fifth_responses{o} = mean(temp_responses,2); 
end

for o=1:size(sixth_z_scores,1)   % neuron number. Calculates average baseline of each neuron    
      %First calculate third z scores. trial num may be different so need
      %separate loops
      cur_neuron = sixth_z_scores{o,:}; %pulls DF/F for each frame for 1 neuron
      temp_baselines = cur_neuron(:,baseline(1):baseline(2));     % calcualte average vuale within baseline period
      mean_sixth_baseline{o} = mean(temp_baselines,2);
      temp_responses = cur_neuron(:,stimulus(1):stimulus(2));     % calcualte average vuale within baseline period
      mean_sixth_responses{o} = mean(temp_responses,2); 
end
%% permutation test first liquid TODO
% For each neuron (20 baseline 20 responses), calculate difference between
% response and baseline. Find average difference. Next, shuffle labels and
% calculate mean difference. Repeat 10,000 times, then see how many
% permutations are less than the actual value, the percent is p. Might be
% better to do area under the curve

%% try normalizing z scores by subtracting baseline




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
delete = intersect(first_responding_neurons,all_neurons_to_delete);
delete_inds = [];
for d = delete
    delete_ind = find(first_responding_neurons == d);
    delete_inds = [delete_inds;delete_ind];
end
first_responding_neurons(delete_inds) = [];

disp(horzcat(first_liquid, ' responding neurons:'));
disp((horzcat(num2str(length(first_responding_neurons)), ' out of ', num2str(n - length(all_neurons_to_delete)))))
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
delete = intersect(second_responding_neurons,all_neurons_to_delete);
delete_inds = [];
for d = delete
    delete_ind = find(second_responding_neurons == d);
    delete_inds = [delete_inds;delete_ind];
end
second_responding_neurons(delete_inds) = [];

disp(horzcat(second_liquid, ' responding neurons:'));
disp((horzcat(num2str(length(second_responding_neurons)), ' out of ', num2str(n - length(all_neurons_to_delete)))))
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


%% Do third liquid responses
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

%delete the bad quality neurons you flagged earlier
delete = intersect(third_responding_neurons,all_neurons_to_delete);
delete_inds = [];
for d = delete
    delete_ind = find(third_responding_neurons == d);
    delete_inds = [delete_inds;delete_ind];
end
third_responding_neurons(delete_inds) = [];

disp(horzcat(third_liquid, ' responding neurons:'));
disp((horzcat(num2str(length(third_responding_neurons)), ' out of ', num2str(n - length(all_neurons_to_delete)))))
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

%% fourth liquid
p_values_fourth = [];
excitations = [];
for n = 1:size(mean_fourth_responses, 2)
    actual_baselines = mean_fourth_baseline{n};
    actual_responses = mean_fourth_responses{n};
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
    p_values_fourth = horzcat(p_values_fourth, p_value);
end

fourth_responding_neurons = find(p_values_fourth < .05);

%delete the bad quality neurons you flagged earlier
delete = intersect(fourth_responding_neurons,all_neurons_to_delete);
delete_inds = [];
for d = delete
    delete_ind = find(fourth_responding_neurons == d);
    delete_inds = [delete_inds;delete_ind];
end
fourth_responding_neurons(delete_inds) = [];

disp(horzcat(fourth_liquid, ' responding neurons:'));
disp((horzcat(num2str(length(fourth_responding_neurons)), ' out of ', num2str(n - length(all_neurons_to_delete)))))
disp(num2str(fourth_responding_neurons));

fourth_excited_neurons = [];
fourth_inhibited_neurons = [];
for resp = fourth_responding_neurons
    if excitations(resp) == 1
        fourth_excited_neurons = horzcat(fourth_excited_neurons, resp);
    else
        fourth_inhibited_neurons = horzcat(fourth_inhibited_neurons, resp);
    end
end

disp(horzcat('Number of excited neurons: ', num2str(size(fourth_excited_neurons, 2))));
disp(horzcat('Number of inhibited neurons: ', num2str(size(fourth_inhibited_neurons, 2))));

%% Fifth liquid
p_values_fifth = [];
excitations = [];
for n = 1:size(mean_fifth_responses, 2)
    actual_baselines = mean_fifth_baseline{n};
    actual_responses = mean_fifth_responses{n};
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
    p_values_fifth = horzcat(p_values_fifth, p_value);
end

fifth_responding_neurons = find(p_values_fifth < .05);

%delete the bad quality neurons you flagged earlier
delete = intersect(fifth_responding_neurons,all_neurons_to_delete);
delete_inds = [];
for d = delete
    delete_ind = find(fifth_responding_neurons == d);
    delete_inds = [delete_inds;delete_ind];
end
fifth_responding_neurons(delete_inds) = [];

disp(horzcat(fifth_liquid, ' responding neurons:'));
disp((horzcat(num2str(length(fifth_responding_neurons)), ' out of ', num2str(n - length(all_neurons_to_delete)))))
disp(num2str(fifth_responding_neurons));

fifth_excited_neurons = [];
fifth_inhibited_neurons = [];
for resp = fifth_responding_neurons
    if excitations(resp) == 1
        fifth_excited_neurons = horzcat(fifth_excited_neurons, resp);
    else
        fifth_inhibited_neurons = horzcat(fifth_inhibited_neurons, resp);
    end
end

disp(horzcat('Number of excited neurons: ', num2str(size(fifth_excited_neurons, 2))));
disp(horzcat('Number of inhibited neurons: ', num2str(size(fifth_inhibited_neurons, 2))));

%% sixth liquid
p_values_sixth = [];
excitations = [];
for n = 1:size(mean_sixth_responses, 2)
    actual_baselines = mean_sixth_baseline{n};
    actual_responses = mean_sixth_responses{n};
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
    p_values_sixth = horzcat(p_values_sixth, p_value);
end

sixth_responding_neurons = find(p_values_sixth < .05);

%delete the bad quality neurons you flagged earlier
delete = intersect(sixth_responding_neurons,all_neurons_to_delete);
delete_inds = [];
for d = delete
    delete_ind = find(sixth_responding_neurons == d);
    delete_inds = [delete_inds;delete_ind];
end
sixth_responding_neurons(delete_inds) = [];

disp(horzcat(sixth_liquid, ' responding neurons:'));
disp((horzcat(num2str(length(sixth_responding_neurons)), ' out of ', num2str(n - length(all_neurons_to_delete)))))
disp(num2str(sixth_responding_neurons));

sixth_excited_neurons = [];
sixth_inhibited_neurons = [];
for resp = sixth_responding_neurons
    if excitations(resp) == 1
        sixth_excited_neurons = horzcat(sixth_excited_neurons, resp);
    else
        sixth_inhibited_neurons = horzcat(sixth_inhibited_neurons, resp);
    end
end

disp(horzcat('Number of excited neurons: ', num2str(size(sixth_excited_neurons, 2))));
disp(horzcat('Number of inhibited neurons: ', num2str(size(sixth_inhibited_neurons, 2))));


first_and_third_responding_neurons = intersect(first_responding_neurons, third_responding_neurons);
first_and_second_responding_neurons = intersect(first_responding_neurons, second_responding_neurons);
second_and_third_responding_neurons = intersect(second_responding_neurons, third_responding_neurons);

total_neurons = n - size(all_neurons_to_delete, 2);
percent_first_responding = length(first_responding_neurons)/total_neurons;
percent_second_responding = length(second_responding_neurons)/total_neurons;
percent_third_responding = length(third_responding_neurons)/total_neurons;

percent_second_to_first = length(first_and_second_responding_neurons)/length(second_responding_neurons);
percent_second_to_third = length(second_and_third_responding_neurons)/length(second_responding_neurons);

percent_first_to_second = length(first_and_second_responding_neurons)/length(first_responding_neurons);
percent_first_to_third = length(first_and_third_responding_neurons)/length(first_responding_neurons);

%Find those responding to only one thing
copy_first_responding_neurons = first_responding_neurons;
copy_second_responding_neurons = second_responding_neurons;
copy_third_responding_neurons = third_responding_neurons;

%% Delete first responding neurons responding to multiple things

%find first and third
used_inds = [];
for k = first_and_third_responding_neurons
    rep = find(copy_first_responding_neurons == k);
    used_inds = horzcat(used_inds, rep);
end

%find first and second
for k = first_and_second_responding_neurons
    rep = find(copy_first_responding_neurons == k);
    used_inds = horzcat(used_inds, rep);
end

used_inds = unique(used_inds);
%delete first neurons that also respond to third or second
copy_first_responding_neurons(used_inds) = [];
first_only_neurons = copy_first_responding_neurons;

%% Delete second responding neurons responding to multiple things
%find first and third
used_inds = [];
for k = second_and_third_responding_neurons
    rep = find(copy_second_responding_neurons == k);
    used_inds = horzcat(used_inds, rep);
end

%find first and second
for k = first_and_second_responding_neurons
    rep = find(copy_second_responding_neurons == k);
    used_inds = horzcat(used_inds, rep);
end

used_inds = unique(used_inds);
%delete first neurons that also respond to third or second
copy_second_responding_neurons(used_inds) = [];
second_only_neurons = copy_second_responding_neurons;

%% Delete third responding neurons responding to multiple things
%find first and third
used_inds = [];
for k = second_and_third_responding_neurons
    rep = find(copy_third_responding_neurons == k);
    used_inds = horzcat(used_inds, rep);
end

%find first and second
for k = first_and_third_responding_neurons
    rep = find(copy_third_responding_neurons == k);
    used_inds = horzcat(used_inds, rep);
end

used_inds = unique(used_inds);
%delete first neurons that also respond to third or second
copy_third_responding_neurons(used_inds) = [];
third_only_neurons = copy_third_responding_neurons;

save(horzcat('responsive_neurons_session_', num2str(session_ind), '.mat'), 'first_responding_neurons', 'second_responding_neurons', ...
    'third_responding_neurons', 'first_and_third_responding_neurons', ...
    'first_and_second_responding_neurons', 'second_and_third_responding_neurons', ...
    'percent_first_responding', 'percent_second_responding', 'percent_third_responding', ...
    'percent_second_to_first', 'percent_second_to_third', 'percent_first_to_second', ...
    'percent_first_to_third', 'first_only_neurons', 'second_only_neurons', ...
    'third_only_neurons', 'first_liquid', 'second_liquid', 'third_liquid', ...
    'second_excited_neurons', 'second_inhibited_neurons', 'first_excited_neurons', ...
    'first_inhibited_neurons', 'third_inhibited_neurons', 'third_excited_neurons', 'all_neurons_to_delete');


%% TODO now it is calculating average response of each neuron, then averaging 
%those and taking the difference from baseline. You want to find average
%response of each neuron- just redo it with the z scores
%% Calculate permutation test (I'll want to do dif in baseline to stim period
