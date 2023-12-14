function process_extracted_data()
% Extract the licking data
graph_behavior_Ale_V7();

% Segment trials and group by liquid delivered, excluding trials with no
% licks
find_liquid_responses_V2();

% Find neurons that respond to each stimulus
permutation_test_Sara();

% Plot avg neuron responses
plot_average_responsive_neurons();

% Run this right after running pca in previous function
plot_pca_trajectories();
end