function baseline_analysis_pbs_ghrelin()
%% Calculates z score using control session (after pbs injection) as baselin
% written by Sara Boyle
%% Sessions are 10 25 second trials 20 minutes after 1 x pbs subq injection, then
% 10 25 second trials after ghrelin injection. 1.5 minutes ITI

close all
clear
trial_start = 1;
trial_end = 250;
trial_num_pbs = 10;
trial_num_ghrelin = 10;
format short g
% Find the deleted neuron indeces for each session. Responsive neurons are
% calculated before removal of bad neurons, so you will need to add in bad
% neurons to get accurate responsiveness flags for each session.

deleted_neurons = {};
dFFs = {};
ZScores = {};
trial_lengths = [];
trial_indss = {};
relative_sorted_food_timess = {};

load('z_scores_session_1.mat');
deleted_neurons = neurons_to_delete;
NeurondFF(neurons_to_delete, :) = [];
NeuronZScores_baseline_normed(neurons_to_delete, :) = [];
ZScores = NeuronZScores_baseline_normed;
dFFs = NeurondFF; 


% plot whole session from tracked neurons
first_dFF = dFFs;

% means done implementing
harmonic_means_pre = [];
means_pre = [];

% stds done implementing
harmonic_stds_pre = [];
harmonic_stds_post = [];
stds_pre = [];
stds_post = [];

% avg z score done implementing
harmonic_zs_pre = [];
harmonic_zs_post = [];
zs_pbs = [];
zs_ghrelin = [];

% % time above dF/F mean thresh done implementing
harm_over_mean_pre = [];
over_mean_pre = [];
times_above_thresh_pbs = [];
times_above_thresh_ghrelin = [];

% % time above 1 std from dF/F mean thresh done implementing
harmonic_std_over_means_pre = [];
std_over_means_pre = [];

% integral above 1 std + mean done
int_above_mean_thresh_pbs = [];
int_above_mean_thresh_ghrelin = [];
h_int_above_mean_thresh_pre =  [];
p_values_first = [];
excitations = [];
p_values_time = [];
excitations_time = [];
sig_inhibited = [];
sig_excited = [];
sig_inhibited_time = [];
sig_excited_time = [];
sig_inhibited_both = [];
sig_excited_both = [];

for n = 1:size(dFFs, 1)
    first_dFF_neuron = first_dFF{n, :};
    first_dFF_neuron = first_dFF_neuron(:, trial_start:trial_end);
    pbs_trials = first_dFF_neuron(1:trial_num_pbs, :);
    ghrelin_trials = first_dFF_neuron((1 + trial_num_pbs):(trial_num_pbs + trial_num_ghrelin), :);
    
    
    flat_pbs_trials = reshape(pbs_trials', 1, size(pbs_trials, 1) * size(pbs_trials, 2));
    flat_ghrelin_trials = reshape(ghrelin_trials', 1, size(ghrelin_trials, 1) * size(ghrelin_trials, 2));
    flat_first_dFF_neuron = reshape(first_dFF_neuron', 1, size(first_dFF_neuron, 1) * size(first_dFF_neuron, 2));
    
    % make sure same length
    smallest = size(flat_first_dFF_neuron, 2);
    
    %% make sure all numbers are positive
    pooled = flat_first_dFF_neuron;
    min_both = min(pooled) - 1;
    flat_first_dFF_neuron = flat_first_dFF_neuron - min_both;
    
    % first do dF/F with no change
    %dFF_fig = figure;
    %plot(flat_first_dFF_neuron, 'LineWidth', 3);
    %hold on 
    %title('pre vs post HFD');
    
    % calculate means
    %harm_mean_first_dFF = harmmean(flat_first_dFF_neuron);
    mean_pbs_dFF = mean(flat_pbs_trials);
    mean_ghrelin_dFF = mean(flat_ghrelin_trials);
    
    % save means
    %harmonic_means_pre = [harmonic_means_pre, harm_mean_first_dFF];
    %means_pre = [means_pre, mean_first_dFF];
    
    % save time dF/F spends above mean
   
    %harm_above_mean_thresh_pre = size(find(flat_first_dFF_neuron >= harm_mean_first_dFF), 2)/smallest;
    %above_mean_thresh_pre = size(find(flat_first_dFF_neuron >= mean_first_dFF), 2)/smallest;
    
    %harm_over_mean_pre = [harm_over_mean_pre, harm_above_mean_thresh_pre];
    %over_mean_pre = [over_mean_pre, above_mean_thresh_pre];
    
    % plot means
    %plot([0, smallest], [harm_mean_first_dFF, harm_mean_first_dFF]);
    %plot([0, smallest], [mean_first_dFF, mean_first_dFF]);
    %legend({'Pre dF/F', 'h mean pre', 'mean pre'});
    
    %disp(['first arithmetic mean: ', num2str(mean_first_dFF)]);
    %disp(['first harmonic mean: ', num2str(harm_mean_first_dFF)]);
    
    % calculate standard deviations
    std_pbs_dFF = std(flat_pbs_trials);
    std_ghrelin_dFF = std(flat_ghrelin_trials);
  
    %harmonic_std_first_dFF = sqrt(sum((flat_first_dFF_neuron-harm_mean_first_dFF).^2)/(smallest-1));
  
    % save stds
    %harmonic_stds_pre = [harmonic_stds_pre, harmonic_std_first_dFF];
    %stds_pre = [stds_pre, std_first_dFF];
    
    % calculate above 1 std above mean thresh
    %h_std_thresh_pre = (harm_mean_first_dFF + harmonic_std_first_dFF);
    %harmonic_std_over_mean_pre = size(find(flat_first_dFF_neuron >= h_std_thresh_pre), 2)/smallest;
  
    std_thresh_pbs = (mean_pbs_dFF + std_pbs_dFF);
    std_thresh_ghrelin = (mean_ghrelin_dFF + std_ghrelin_dFF);
    %std_over_mean_pre = size(find(flat_first_dFF_neuron >= std_thresh_pre), 2)/smallest;
    
    % save % time over mean + std threshold
    %harmonic_std_over_means_pre = [harmonic_std_over_means_pre, 100* harmonic_std_over_mean_pre];
    %std_over_means_pre = [std_over_means_pre, 100* std_over_mean_pre];
    
    %disp(['first arithmetic std: ', num2str(std_first_dFF)]);
    %disp(['first harmonic std: ', num2str(harmonic_std_first_dFF)]);
    
    %% save integral above threshold
    %avg_harm_above_mean_thresh_pre = flat_first_dFF_neuron >= h_std_thresh_pre;
    %avg_harm_above_mean_thresh_pre = (sum(flat_first_dFF_neuron(avg_harm_above_mean_thresh_pre)) - (h_std_thresh_pre * sum(avg_harm_above_mean_thresh_pre)))/smallest;
    
    avg_above_mean_thresh_pbs = flat_pbs_trials >= std_thresh_pbs;
    time_above_pbs = sum(avg_above_mean_thresh_pbs)/(smallest/2);
    times_above_thresh_pbs = [times_above_thresh_pbs, time_above_pbs];
    avg_above_mean_thresh_pbs = (sum(flat_pbs_trials(avg_above_mean_thresh_pbs)) - (std_thresh_pbs * sum(avg_above_mean_thresh_pbs)))/(smallest/2);
    
    avg_above_mean_thresh_ghrelin = flat_ghrelin_trials >= std_thresh_ghrelin;
    time_above_ghrelin = sum(avg_above_mean_thresh_ghrelin)/(smallest/2);
    times_above_thresh_ghrelin = [times_above_thresh_ghrelin, time_above_ghrelin];
    avg_above_mean_thresh_ghrelin = (sum(flat_ghrelin_trials(avg_above_mean_thresh_ghrelin)) - (std_thresh_ghrelin * sum(avg_above_mean_thresh_ghrelin)))/(smallest/2);
    
    int_above_mean_thresh_pbs = [int_above_mean_thresh_pbs, avg_above_mean_thresh_pbs];
    int_above_mean_thresh_ghrelin = [int_above_mean_thresh_ghrelin, avg_above_mean_thresh_ghrelin];
    disp(['pbs int above: ', num2str(avg_above_mean_thresh_pbs)]);
    disp(['ghrelin int above: ', num2str(avg_above_mean_thresh_ghrelin)]);
    disp(['pbs time above: ', num2str(time_above_pbs)]);
    disp(['ghrelin time above: ', num2str(time_above_ghrelin)]);
    
    %h_int_above_mean_thresh_pre = [h_int_above_mean_thresh_pre, avg_harm_above_mean_thresh_pre];
    
    % now do arithmetic z score
    z_score_pbs = (flat_pbs_trials - mean_pbs_dFF)/std_pbs_dFF;
    z_score_ghrelin = (flat_ghrelin_trials - mean_pbs_dFF)/std_pbs_dFF;
    zs_pbs = [zs_pbs, mean(z_score_pbs)];
    zs_ghrelin = [zs_ghrelin, mean(z_score_ghrelin)];
    %plot(horzcat(z_score_pbs, z_score_ghrelin))
    
    % get trial by trial means
    pbs_tri_means = [];
    pbs_tri_above = [];
    bin_shift = 50; % amount to slide bin
    bin_size = 500; % size of bin
    % for each trial include next trial
    art_trial_num = 1 + ((smallest/2) - bin_size)/bin_shift;
    for tri = 1:art_trial_num
        cur = z_score_pbs(1, (1 + (tri - 1)*bin_shift):((1 +(tri - 1)*bin_shift) + bin_size - 1));
        dFF_cur = flat_pbs_trials(1, (1 + (tri - 1)*bin_shift):((1 +(tri - 1)*bin_shift) + bin_size - 1));
        tri_mean = mean(cur);
        avg_above_mean_thresh_pbs = dFF_cur >= std_thresh_pbs;
        tri_above = sum(avg_above_mean_thresh_pbs)/(500);
        pbs_tri_above = [pbs_tri_above, tri_above];
        pbs_tri_means = [pbs_tri_means, tri_mean];
    end
    
    ghrelin_tri_means = [];
    ghrelin_tri_above = [];
    % split into overlapping bins to increase "trial" size
    % now it is 9 trials, length 500, 250 frame dif. could be 45 trials, 50
    % dif
    for tri = 1:art_trial_num
        cur = z_score_ghrelin(1, (1 + (tri - 1)*bin_shift):((1 +(tri - 1)*bin_shift) + bin_size - 1));
        dFF_cur = flat_ghrelin_trials(1, (1 + (tri - 1)*bin_shift):((1 +(tri - 1)*bin_shift) + bin_size - 1));
        tri_mean = mean(cur);
        avg_above_mean_thresh_ghrelin = dFF_cur >= std_thresh_ghrelin;
        tri_above = sum(avg_above_mean_thresh_ghrelin)/(size(cur, 2));
        ghrelin_tri_above = [ghrelin_tri_above, tri_above];
        ghrelin_tri_means = [ghrelin_tri_means, tri_mean];
    end
    
    %% permutation test
    mean_first_baseline = pbs_tri_means;
    mean_first_responses = ghrelin_tri_means;
    mean_above_pbs = pbs_tri_above;
    mean_above_ghrelin = ghrelin_tri_above;
    %% do permutation test with avg z score values
    actual_baselines = mean_first_baseline;
    actual_responses = mean_first_responses;
    pos = 0 < mean(actual_responses - actual_baselines);
    actual_dif = abs(mean(actual_responses - actual_baselines));

    to_shuffle = vertcat(actual_baselines', actual_responses');
    shuffled_difs = zeros(1, 100000);
    for g = 1:100000
        order = randperm(size(to_shuffle, 1));
        rand_baselines = to_shuffle(order(1:size(actual_baselines, 2)));
        rand_responses = to_shuffle(order((1 + size(actual_baselines, 2):end)));
        shuffled_dif = abs(mean(rand_responses - rand_baselines));
        shuffled_difs(g) = shuffled_dif;
    end

    biggers = find(shuffled_difs > actual_dif); %random difs bigger than actual
    p_value = length(biggers)/100000;
    if p_value < 0.05
        disp(['Significantly different z score response from ', num2str(n)]);
        if pos
            disp('positive response');
            sig_excited = [sig_excited, n];
            
        else
            disp('inhibited response');
            sig_inhibited = [sig_inhibited, n];
        end
    else
        disp(['No change in z score response from ', num2str(n)])
    end

    p_values_first = horzcat(p_values_first, p_value);
    excitations = horzcat(excitations, pos);
    
    
    
    %% Do permutation test with % time above
    actual_baselines = mean_above_pbs;
    actual_responses = mean_above_ghrelin;
    pos = 0 < mean(actual_responses - actual_baselines);
    actual_dif = abs(mean(actual_responses - actual_baselines));

    to_shuffle = vertcat(actual_baselines', actual_responses');
    shuffled_difs = zeros(1, 100000);
    for g = 1:100000
        order = randperm(size(to_shuffle, 1));
        rand_baselines = to_shuffle(order(1:size(actual_baselines, 2)));
        rand_responses = to_shuffle(order((1 + size(actual_baselines, 2):end)));
        shuffled_dif = abs(mean(rand_responses - rand_baselines));
        shuffled_difs(g) = shuffled_dif;
    end

    biggers = find(shuffled_difs > actual_dif); %random difs bigger than actual
    p_value = length(biggers)/100000;
    if p_value < 0.05
        disp(['Significantly different % active time response from ', num2str(n)]);
        if pos
            disp('positive response');
            sig_excited_time = [sig_excited_time, n];
            if ismember(n, sig_excited)
                sig_excited_both = [sig_excited_both, n];
            end
        else
            disp('inhibited response');
            sig_inhibited_time = [sig_inhibited_time, n];
            if ismember(n, sig_inhibited)
                sig_inhibited_both = [sig_inhibited_both, n];
            end
        end
    else
        disp(['No change in time over response from ', num2str(n)])
    end

    p_values_time = horzcat(p_values_time, p_value);
    excitations_time = horzcat(excitations_time, pos);
  
   
    %harm_z_score_first = (flat_first_dFF_neuron - harm_mean_first_dFF)/harmonic_std_first_dFF;
    
    %harmonic_zs_pre = [harmonic_zs_pre, mean(harm_z_score_first)];
    %zs_pre = [zs_pre, mean(z_score_first)];
    
    %figure
    %plot(z_score_first, 'LineWidth', 3);
    %hold on 
    %plot(harm_z_score_first, 'LineWidth', 3);
    %title('Z scores vs harmonic z scores');
    %legend({'Z Score', 'Harmonic Z Score'});
    
    % now do harmonic z scores
    close all
end
%%

save('ghrelin_stats.mat', 'means_pre', 'stds_pre', 'zs_pbs', 'zs_ghrelin', ...
    'over_mean_pre', 'std_over_means_pre', 'int_above_mean_thresh_ghrelin', ...
    'int_above_mean_thresh_pbs', 'p_values_first', 'p_values_time', 'excitations', ...
    'excitations_time', 'sig_excited', 'sig_inhibited', 'sig_inhibited_both', ...
    'sig_excited_both', 'sig_excited_time', 'sig_inhibited_time');

