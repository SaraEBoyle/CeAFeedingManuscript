%% define parameter %%
clear
%% Calculates means baselines of all neurons in amph condition, minus mean baselines in all neurons in sal condition %%
%% Doesn't explicitly delete neurons_to_delete, but excludes those from responsive flags
%Preserves spatial info, leaving those duds in the z score data
%tem = dir('*Liquids');
first_food = 'High Fat';
second_food = 'Sucrose';
third_food = 'Missed';
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

%find trial number for each sessio
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
%load('neuron.mat');
load(horzcat('z_scores_session_', num2str(session_ind), '.mat'));
framerate = 10;
%NeuronZScores = NeuronZScores_baseline_normed;
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

%[trial_inds, relative_sorted_food_times, labels] = read_in_cam_data();
load('time_of_food_contact.mat');

%if an airpuff session
TrialTypes = trial_inds;
trial_n = [1, 2, 3];


baseline = [3 39]; %3 seconds right before the stimulus
stimulus = [60 200]; %100-125 is the CS-- tone. 126-150,160 is the US-- water/shock
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

%% mean baseline per neuron in first liquid
mean_first_baseline = {};
mean_first_responses = {};
mean_second_baseline = {};
mean_second_responses = {};
mean_third_baseline = {};
mean_third_responses = {};

first_z_scores = NeuronZScores_aligned(:, 1);
second_z_scores = NeuronZScores_aligned(:, 2);
third_z_scores = NeuronZScores_aligned(:, 3);
for o=1:size(first_z_scores,1)   % neuron number. Calculates average baseline of each neuron    
      %First calculate first z scores. trial num may be different so need
      %separate loops
      
      %Need second loop since delivery times are different trial by trial
      high_trials = find(trial_inds == 1);
      cur_neuron = first_z_scores{o,:}; %pulls DF/F for each frame for 1 neuron
      temp_baselines = [];
      temp_responses = [];
      %for each high fat trial
      for tri = 1:size(cur_neuron, 1)
          cur_trial = cur_neuron(tri, :);
          cur_trial = cur_trial - mean(cur_trial(baseline(1):baseline(2)));
          %stimulus(1) = round(relative_sorted_food_times(high_trials(tri), 1) * framerate);
          temp_baselines = [temp_baselines;mean(cur_trial(baseline(1):baseline(2)))];     % calcualte average vuale within baseline period
          temp_responses = [temp_responses;mean(cur_trial(stimulus(1):stimulus(2)))];     % calcualte average vuale within baseline period
          %stimulus(1) = round(relative_sorted_food_times(high_trials(tri), 1) * framerate);
          %temp_baselines = [temp_baselines;mean(cur_neuron(tri,baseline(1):baseline(2)))];     % calcualte average vuale within baseline period
          %temp_responses = [temp_responses;mean(cur_neuron(tri,stimulus(1):stimulus(2)))];     % calcualte average vuale within baseline period
      end
      mean_first_baseline{o} = mean(temp_baselines,2);
      mean_first_responses{o} = mean(temp_responses,2); 
end

for o=1:size(second_z_scores,1)   % neuron number. Calculates average baseline of each neuron    
      %First calculate second z scores. trial num may be different so need
      %separate loops
      
      %Need second loop since delivery times are different trial by trial
      low_trials = find(trial_inds == 2);
      cur_neuron = second_z_scores{o,:}; %pulls DF/F for each frame for 1 neuron
      temp_baselines = [];
      temp_responses = [];
      %for each high fat trial
      for tri = 1:size(cur_neuron, 1)
          cur_trial = cur_neuron(tri, :);
          cur_trial = cur_trial - mean(cur_trial(baseline(1):baseline(2)));
          %stimulus(1) = round(relative_sorted_food_times(low_trials(tri), 1) * framerate);
          temp_baselines = [temp_baselines;mean(cur_trial(baseline(1):baseline(2)))];     % calcualte average vuale within baseline period
          temp_responses = [temp_responses;mean(cur_trial(stimulus(1):stimulus(2)))];     % calcualte average vuale within baseline period
          %temp_baselines = [temp_baselines;mean(cur_neuron(tri,baseline(1):baseline(2)))];     % calcualte average vuale within baseline period
          %temp_responses = [temp_responses;mean(cur_neuron(tri,stimulus(1):stimulus(2)))];     % calcualte average vuale within baseline period
      end
      mean_second_baseline{o} = mean(temp_baselines,2);
      mean_second_responses{o} = mean(temp_responses,2); 
end

for o=1:size(third_z_scores,1)   % neuron number. Calculates average baseline of each neuron    
      %First calculate third z scores. trial num may be different so need
      %separate loops
      
      %Need second loop since delivery times are different trial by trial
      miss_trials = find(trial_inds == 3);
      cur_neuron = third_z_scores{o,:}; %pulls DF/F for each frame for 1 neuron
      temp_baselines = [];
      temp_responses = [];
      %for each high fat trial
      for tri = 1:size(cur_neuron, 1)
          cur_trial = cur_neuron(tri, :);
          cur_trial = cur_trial - mean(cur_trial(baseline(1):baseline(2)));
          %stimulus(1) = round(relative_sorted_food_times(miss_trials(tri), 1) * framerate);
          temp_baselines = [temp_baselines;mean(cur_trial(baseline(1):baseline(2)))];     % calcualte average vuale within baseline period
          temp_responses = [temp_responses;mean(cur_trial(stimulus(1):stimulus(2)))];     % calcualte average vuale within baseline period
          %temp_baselines = [temp_baselines;mean(cur_neuron(tri,baseline(1):baseline(2)))];     % calcualte average vuale within baseline period
          %temp_responses = [temp_responses;mean(cur_neuron(tri,stimulus(1):stimulus(2)))];     % calcualte average vuale within baseline period
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
delete = intersect(first_responding_neurons, neurons_to_delete);
delete_inds = [];
for d = delete
    delete_ind = find(first_responding_neurons == d);
    delete_inds = [delete_inds;delete_ind];
end
first_responding_neurons(delete_inds) = [];

disp('High fat responding neurons:');
disp((horzcat(num2str(length(first_responding_neurons)), ' out of ', num2str(n - size(neurons_to_delete, 2)))))
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
delete = intersect(second_responding_neurons, neurons_to_delete);
delete_inds = [];
for d = delete
    delete_ind = find(second_responding_neurons == d);
    delete_inds = [delete_inds;delete_ind];
end
second_responding_neurons(delete_inds) = [];

disp('Low fat responding neurons:');
disp((horzcat(num2str(length(second_responding_neurons)), ' out of ', num2str(n - size(neurons_to_delete, 2)))))
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

%delete the bad quality neurons you flagged earlier
delete = intersect(third_responding_neurons, neurons_to_delete);
delete_inds = [];
for d = delete
    delete_ind = find(third_responding_neurons == d);
    delete_inds = [delete_inds;delete_ind];
end
third_responding_neurons(delete_inds) = [];

disp(horzcat('Missed trial responding neurons:'));
disp((horzcat(num2str(length(third_responding_neurons)), ' out of ', num2str(n - size(neurons_to_delete, 2)))))
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


first_and_third_responding_neurons = intersect(first_responding_neurons, third_responding_neurons);
first_and_second_responding_neurons = intersect(first_responding_neurons, second_responding_neurons);
second_and_third_responding_neurons = intersect(second_responding_neurons, third_responding_neurons);

total_neurons = n - size(neurons_to_delete, 2);
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
    'third_only_neurons', 'first_food', 'second_food', 'third_food', ...
    'second_excited_neurons', 'second_inhibited_neurons', 'first_excited_neurons', ...
    'first_inhibited_neurons', 'third_inhibited_neurons', 'third_excited_neurons', 'TrialTypes', 'relative_sorted_food_times');


%% TODO now it is calculating average response of each neuron, then averaging 
%those and taking the difference from baseline. You want to find average
%response of each neuron- just redo it with the z scores
%% Calculate permutation test (I'll want to do dif in baseline to stim period
