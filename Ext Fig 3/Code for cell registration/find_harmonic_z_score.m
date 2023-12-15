function find_harmonic_z_score(baseline_cell)
%% takes a baseline cell n neurons by 1, each entry is an array of o trials by
% p frames. Calculates harmonic z score by using the harmonic mean in all calculations
% instead of arithmetic mean. 

% so, harmonic standard deviation = sqrt(sum(xi - harmonic mean)/N)
% z score = (xi - harmonic mean)/harmonic standard deviation

%% here is an example mouse
close all
clear
load('responses_in_tracked_cells.mat');
results_folder = uigetdir(pwd);

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
end


% plot whole session from tracked neurons
first_dFF = dFFs{1};
second_dFF = dFFs{2};

% means done implementing
harmonic_means_pre = [];
harmonic_means_post = [];
means_pre = [];
means_post = [];

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

% time above dF/F mean thresh done implementing
harm_over_mean_pre = [];
harm_over_mean_post = [];
over_mean_pre = [];
over_mean_post = [];

% time above 1 std from dF/F mean thresh done implementing
harmonic_std_over_means_pre = [];
harmonic_std_over_means_post = [];
std_over_means_pre = [];
std_over_means_post = [];

% average dF/F above mean + std thresh (so subract thresh)
harmonic_avg_std_over_means_pre = [];
harmonic_avg_std_over_means_post = [];
avg_std_over_means_pre = [];
avg_std_over_means_post = [];




for n = 1:size(indeces_in_both_sessions, 1)
    first_ind = indeces_in_both_sessions(n, 1);
    first_dFF_neuron = first_dFF{first_ind, 1};
    flat_first_dFF_neuron = reshape(first_dFF_neuron', 1, size(first_dFF_neuron, 1) * size(first_dFF_neuron, 2));
    
    second_ind = indeces_in_both_sessions(n, 2);
    second_dFF_neuron = second_dFF{second_ind, 1};
    flat_second_dFF_neuron = reshape(second_dFF_neuron', 1, size(second_dFF_neuron, 1) * size(second_dFF_neuron, 2));
    
    % make sure same length
    smallest = min(size(flat_second_dFF_neuron, 2), size(flat_first_dFF_neuron, 2));
    flat_first_dFF_neuron = flat_first_dFF_neuron(1:smallest);
    flat_second_dFF_neuron = flat_second_dFF_neuron(1:smallest);
    
    %% make sure all numbers are positive
    pooled = [flat_first_dFF_neuron, flat_second_dFF_neuron];
    min_both = min(pooled) - 1;
    flat_first_dFF_neuron = flat_first_dFF_neuron - min_both;
    flat_second_dFF_neuron = flat_second_dFF_neuron - min_both;
    
    % first do dF/F with no change
   
    
    
    
    % calculate means
    harm_mean_first_dFF = harmmean(flat_first_dFF_neuron);
    harm_mean_second_dFF = harmmean(flat_second_dFF_neuron);
    mean_first_dFF = mean(flat_first_dFF_neuron);
    mean_second_dFF = mean(flat_second_dFF_neuron);
    
    % save means
    harmonic_means_pre = [harmonic_means_pre, harm_mean_first_dFF];
    harmonic_means_post = [harmonic_means_post, harm_mean_second_dFF];
    means_pre = [means_pre, mean_first_dFF];
    means_post = [means_post, mean_second_dFF];
    
    % plot pre
    dFF_fig = figure;
    plot(flat_first_dFF_neuron, 'LineWidth', 3, 'Color', 'b');
    hold on 
    plot([0, smallest], [harm_mean_first_dFF, harm_mean_first_dFF], 'Color', 'g', 'LineWidth', 2);
    plot([0, smallest], [mean_first_dFF, mean_first_dFF], 'Color', 'b', 'LineWidth', 2);
    title('pre HFD');
    legend({'Pre HFD', 'Harmonic Mean', 'Arithmetic Mean'});
    
    figure
    plot(flat_second_dFF_neuron, 'LineWidth', 3, 'Color', 'b');
    hold on
    plot([0, smallest], [harm_mean_second_dFF, harm_mean_second_dFF], 'Color', 'g', 'LineWidth', 2);
    plot([0, smallest], [mean_second_dFF, mean_second_dFF], 'Color', 'b', 'LineWidth', 2);
    legend({'Post HFD', 'Harmonic Mean', 'Arithmetic Mean'});
    title('Post HFD');
    % save time dF/F spends above mean
   
    harm_above_mean_thresh_pre = size(find(flat_first_dFF_neuron >= harm_mean_first_dFF), 2)/smallest;
    harm_above_mean_thresh_post = size(find(flat_second_dFF_neuron >= harm_mean_second_dFF), 2)/smallest;
    above_mean_thresh_pre = size(find(flat_first_dFF_neuron >= mean_first_dFF), 2)/smallest;
    above_mean_thresh_post = size(find(flat_second_dFF_neuron >= mean_second_dFF), 2)/smallest;
    
    harm_over_mean_pre = [harm_over_mean_pre, harm_above_mean_thresh_pre];
    harm_over_mean_post = [harm_over_mean_post, harm_above_mean_thresh_post];
    over_mean_pre = [over_mean_pre, above_mean_thresh_pre];
    over_mean_post = [over_mean_post, above_mean_thresh_post];
    
   
    
    % plot means
    
    
    legend({'Pre dF/F', 'Post dF/F', 'h mean pre', ' h mean post', 'mean pre', 'mean post'});
    
    disp(['first arithmetic mean: ', num2str(mean_first_dFF)]);
    disp(['second arithmetic mean: ', num2str(mean_second_dFF)]);
    disp(['first harmonic mean: ', num2str(harm_mean_first_dFF)]);
    disp(['second harmonic mean: ', num2str(harm_mean_second_dFF)]);
    
    % calculate standard deviations
    std_first_dFF = sqrt(sum((flat_first_dFF_neuron-mean_first_dFF).^2)/(smallest-1));
    std_second_dFF = sqrt(sum((flat_second_dFF_neuron-mean_second_dFF).^2)/(smallest-1));
    
    harmonic_std_first_dFF = sqrt(sum((flat_first_dFF_neuron-harm_mean_first_dFF).^2)/(smallest-1));
    harmonic_std_second_dFF = sqrt(sum((flat_second_dFF_neuron-harm_mean_second_dFF).^2)/(smallest-1));
    
    % save stds
    harmonic_stds_pre = [harmonic_stds_pre, harmonic_std_first_dFF];
    harmonic_stds_post = [harmonic_stds_post, harmonic_std_second_dFF];
    stds_pre = [stds_pre, std_first_dFF];
    stds_post = [stds_post, std_second_dFF];
    
    % calculate above 1 std above mean thresh
    h_std_thresh_pre = (harm_mean_first_dFF + harmonic_std_first_dFF);
    harmonic_std_over_mean_pre = size(find(flat_first_dFF_neuron >= h_std_thresh_pre), 2)/smallest;
    h_std_thresh_post = (harm_mean_second_dFF + harmonic_std_second_dFF);
    harmonic_std_over_mean_post = size(find(flat_second_dFF_neuron >= h_std_thresh_post), 2)/smallest;
    
    std_thresh_pre = (mean_first_dFF + std_first_dFF);
    
    std_over_mean_pre = size(find(flat_first_dFF_neuron >= std_thresh_pre), 2)/smallest;
    
    std_thresh_post = (mean_second_dFF + std_second_dFF);
    std_over_mean_post = size(find(flat_second_dFF_neuron >= std_thresh_post), 2)/smallest;
    
    % save std thresholds
    harmonic_std_over_means_pre = [harmonic_std_over_means_pre, harmonic_std_over_mean_pre];
    harmonic_std_over_means_post = [harmonic_std_over_means_post, harmonic_std_over_mean_post];
    std_over_means_pre = [std_over_means_pre, std_over_mean_pre];
    std_over_means_post = [std_over_means_post, std_over_mean_post];
    
     %% find integral over threshold
    avg_harm_above_mean_thresh_pre = flat_first_dFF_neuron >= h_std_thresh_pre;
    avg_harm_above_mean_thresh_pre = .1*(sum(flat_first_dFF_neuron(avg_harm_above_mean_thresh_pre)) - (h_std_thresh_pre * sum(avg_harm_above_mean_thresh_pre)));
    
    avg_harm_above_mean_thresh_post = flat_second_dFF_neuron >= h_std_thresh_post;
    avg_harm_above_mean_thresh_post = .1*(sum(flat_second_dFF_neuron(avg_harm_above_mean_thresh_post)) - (h_std_thresh_post * sum(avg_harm_above_mean_thresh_post)));
    
    
    avg_above_mean_thresh_pre = flat_first_dFF_neuron >= std_thresh_pre;
    avg_above_mean_thresh_pre = .1*(sum(flat_first_dFF_neuron(avg_above_mean_thresh_pre)) - (std_thresh_pre * sum(avg_above_mean_thresh_pre)));
    
    avg_above_mean_thresh_post = flat_second_dFF_neuron >= std_thresh_post;
    avg_above_mean_thresh_post = .1*(sum(flat_second_dFF_neuron(avg_above_mean_thresh_post)) - (std_thresh_post * sum(avg_above_mean_thresh_post)));
    
    figure
    plot(flat_first_dFF_neuron + min_both, 'LineWidth', 3, 'Color', 'b');
    hold on 
    plot([1, smallest], [h_std_thresh_pre + min_both, h_std_thresh_pre + min_both], 'LineWidth', 3, 'Color', 'g');
    plot([1, smallest], [std_thresh_pre + min_both, std_thresh_pre + min_both], 'LineWidth', 3, 'Color', 'b');
    title('Pre HFD dF/F');
    legend({'Pre HFD Signal', 'Harmonic Mean + STD', 'Arithmetic Mean + STD'});
    
    figure
    plot(flat_second_dFF_neuron + min_both, 'LineWidth', 3, 'Color', 'b');
    hold on 
    plot([1, smallest], [h_std_thresh_post + min_both, h_std_thresh_post + min_both], 'LineWidth', 3, 'Color', 'g');
    plot([1, smallest], [std_thresh_post + min_both, std_thresh_post + min_both], 'LineWidth', 3, 'Color', 'b');
    title('Post HFD dF/F');
    legend({'Post HFD Signal', 'Harmonic Mean + STD', 'Arithmetic Mean + STD'});
    
    
    avg_harm_above_mean_thresh_pre = flat_first_dFF_neuron >= h_std_thresh_pre;
    avg_harm_above_mean_thresh_pre = (sum(flat_first_dFF_neuron(avg_harm_above_mean_thresh_pre)) - (h_std_thresh_pre * sum(avg_harm_above_mean_thresh_pre)))/smallest;
    
    avg_above_mean_thresh_pre = flat_first_dFF_neuron >= std_thresh_pre;
    avg_above_mean_thresh_pre = (sum(flat_first_dFF_neuron(avg_above_mean_thresh_pre)) - (std_thresh_pre * sum(avg_above_mean_thresh_pre)))/smallest;
    
    
    int_above_mean_thresh_pre = [int_above_mean_thresh_pre, avg_above_mean_thresh_pre];
    h_int_above_mean_thresh_pre = [h_int_above_mean_thresh_pre, avg_harm_above_mean_thresh_pre];
    
    harmonic_avg_std_over_means_pre = [harmonic_avg_std_over_means_pre, ];
    harmonic_avg_std_over_means_post = [];
    avg_std_over_means_pre = [];
    avg_std_over_means_post = [];
    
    
    
    
    disp(['first arithmetic std: ', num2str(std_first_dFF)]);
    disp(['second arithmetic std: ', num2str(std_second_dFF)]);
    disp(['first harmonic std: ', num2str(harmonic_std_first_dFF)]);
    disp(['second harmonic std: ', num2str(harmonic_std_second_dFF)]);
    
    % now do arithmetic z score
    z_score_first = (flat_first_dFF_neuron - mean_first_dFF)/std_first_dFF;
    z_score_second = (flat_second_dFF_neuron - mean_second_dFF)/std_second_dFF;
    harm_z_score_first = (flat_first_dFF_neuron - harm_mean_first_dFF)/harmonic_std_first_dFF;
    harm_z_score_second = (flat_second_dFF_neuron - harm_mean_second_dFF)/harmonic_std_second_dFF;
    
    
    harmonic_zs_pre = [harmonic_zs_pre, mean(harm_z_score_first)];
    harmonic_zs_post = [harmonic_zs_post, mean(harm_z_score_second)];
    zs_pre = [zs_pre, mean(z_score_first)];
    zs_post = [zs_post, mean(z_score_second)];
    
    figure
    plot(z_score_first, 'LineWidth', 3);
    hold on 
    plot(z_score_second, 'LineWidth', 3);
    title('pre vs post HFD z scores');
    legend({'Pre Z Score', 'Post Z Score'});
    
    figure
    plot(harm_z_score_first, 'LineWidth', 3);
    hold on 
    plot(harm_z_score_second, 'LineWidth', 3);
    title('pre vs post HFD harmonic z scores');
    legend({'Pre Z Score', 'Post Z Score'});
    
    % now do harmonic z scores
    disp(['z score dif: ', num2str(mean(z_score_second) - mean(z_score_first))]);
    
    disp(['harmonic z score dif: ', num2str(mean(harm_z_score_second) - mean(harm_z_score_first))]);
    close all
end
%%

save('baseline_only_tracked_stats.mat', 'avg_harm_above_mean_thresh_pre', 'avg_harm_above_mean_thresh_post', 'avg_above_mean_thresh_pre', 'avg_above_mean_thresh_post', 'harmonic_means_pre', 'harmonic_means_post', ...
    'means_pre', 'means_post', 'harmonic_stds_pre', 'harmonic_stds_post', 'stds_pre', ...
    'stds_post', 'harmonic_zs_pre', 'harmonic_zs_post', 'zs_pre', 'zs_post', ...
    'harm_over_mean_pre', 'harm_over_mean_post', 'over_mean_pre', 'over_mean_post', ...
    'harmonic_std_over_means_pre', 'harmonic_std_over_means_post', ...
    'std_over_means_pre', 'std_over_means_post');

