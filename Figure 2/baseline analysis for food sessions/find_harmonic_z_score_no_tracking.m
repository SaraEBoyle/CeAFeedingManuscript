function find_harmonic_z_score()
%% takes a baseline cell n neurons by 1, each entry is an array of o trials by
% p frames. Calculates harmonic z score by using the harmonic mean in all calculations
% instead of arithmetic mean. 

%% TODO
% 

% so, harmonic standard deviation = sqrt(sum(xi - harmonic mean)/N)
% z score = (xi - harmonic mean)/harmonic standard deviation

%% here is an example mouse
close all
clear

exclude_food = 1;

if exclude_food == 1
    trial_start = 1;
    trial_end = 39;
else
    trial_start = 1;
    trial_end = 250;
end

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
zs_pre = [];
zs_post = [];

% % time above dF/F mean thresh done implementing
harm_over_mean_pre = [];
over_mean_pre = [];

% % time above 1 std from dF/F mean thresh done implementing
harmonic_std_over_means_pre = [];
std_over_means_pre = [];

% integral above 1 std + mean done
int_above_mean_thresh_pre = [];
h_int_above_mean_thresh_pre =  [];

for n = 1:size(dFFs, 1)
    first_dFF_neuron = first_dFF{n, :};
    first_dFF_neuron = first_dFF_neuron(:, trial_start:trial_end);
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
    harm_mean_first_dFF = harmmean(flat_first_dFF_neuron);
    mean_first_dFF = mean(flat_first_dFF_neuron);
    
    % save means
    harmonic_means_pre = [harmonic_means_pre, harm_mean_first_dFF];
    means_pre = [means_pre, mean_first_dFF];
    
    % save time dF/F spends above mean
   
    harm_above_mean_thresh_pre = size(find(flat_first_dFF_neuron >= harm_mean_first_dFF), 2)/smallest;
    above_mean_thresh_pre = size(find(flat_first_dFF_neuron >= mean_first_dFF), 2)/smallest;
    
    harm_over_mean_pre = [harm_over_mean_pre, harm_above_mean_thresh_pre];
    over_mean_pre = [over_mean_pre, above_mean_thresh_pre];
    
    % plot means
    %plot([0, smallest], [harm_mean_first_dFF, harm_mean_first_dFF]);
    %plot([0, smallest], [mean_first_dFF, mean_first_dFF]);
    %legend({'Pre dF/F', 'h mean pre', 'mean pre'});
    
    disp(['first arithmetic mean: ', num2str(mean_first_dFF)]);
    disp(['first harmonic mean: ', num2str(harm_mean_first_dFF)]);
    
    % calculate standard deviations
    std_first_dFF = sqrt(sum((flat_first_dFF_neuron-mean_first_dFF).^2)/(smallest-1));
  
    harmonic_std_first_dFF = sqrt(sum((flat_first_dFF_neuron-harm_mean_first_dFF).^2)/(smallest-1));
  
    % save stds
    harmonic_stds_pre = [harmonic_stds_pre, harmonic_std_first_dFF];
    stds_pre = [stds_pre, std_first_dFF];
    
    % calculate above 1 std above mean thresh
    h_std_thresh_pre = (harm_mean_first_dFF + harmonic_std_first_dFF);
    harmonic_std_over_mean_pre = size(find(flat_first_dFF_neuron >= h_std_thresh_pre), 2)/smallest;
  
    std_thresh_pre = (mean_first_dFF + std_first_dFF);
    std_over_mean_pre = size(find(flat_first_dFF_neuron >= std_thresh_pre), 2)/smallest;
    
    % save % time over mean + std threshold
    harmonic_std_over_means_pre = [harmonic_std_over_means_pre, 100* harmonic_std_over_mean_pre];
    std_over_means_pre = [std_over_means_pre, 100* std_over_mean_pre];
    
    disp(['first arithmetic std: ', num2str(std_first_dFF)]);
    disp(['first harmonic std: ', num2str(harmonic_std_first_dFF)]);
    
    %% save integral above threshold
    avg_harm_above_mean_thresh_pre = flat_first_dFF_neuron >= h_std_thresh_pre;
    avg_harm_above_mean_thresh_pre = (sum(flat_first_dFF_neuron(avg_harm_above_mean_thresh_pre)) - (h_std_thresh_pre * sum(avg_harm_above_mean_thresh_pre)))/smallest;
    
    avg_above_mean_thresh_pre = flat_first_dFF_neuron >= std_thresh_pre;
    avg_above_mean_thresh_pre = (sum(flat_first_dFF_neuron(avg_above_mean_thresh_pre)) - (std_thresh_pre * sum(avg_above_mean_thresh_pre)))/smallest;
    
    
    int_above_mean_thresh_pre = [int_above_mean_thresh_pre, avg_above_mean_thresh_pre];
    h_int_above_mean_thresh_pre = [h_int_above_mean_thresh_pre, avg_harm_above_mean_thresh_pre];
    
    % now do arithmetic z score
    z_score_first = (flat_first_dFF_neuron - mean_first_dFF)/std_first_dFF;
    harm_z_score_first = (flat_first_dFF_neuron - harm_mean_first_dFF)/harmonic_std_first_dFF;
    
    harmonic_zs_pre = [harmonic_zs_pre, mean(harm_z_score_first)];
    zs_pre = [zs_pre, mean(z_score_first)];
    
    figure
    plot(z_score_first, 'LineWidth', 3);
    hold on 
    plot(harm_z_score_first, 'LineWidth', 3);
    title('Z scores vs harmonic z scores');
    legend({'Z Score', 'Harmonic Z Score'});
    
    % now do harmonic z scores
    close all
end
%%

save('baseline_only_stats_harmonic.mat', 'harmonic_means_pre', 'means_pre', ...
    'harmonic_stds_pre', 'stds_pre', 'harmonic_zs_pre', 'zs_pre', ...
    'harm_over_mean_pre', 'over_mean_pre', 'harmonic_std_over_means_pre', ...
    'std_over_means_pre', 'int_above_mean_thresh_pre', 'h_int_above_mean_thresh_pre');

