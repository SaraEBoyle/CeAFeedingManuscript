%% Pipeline for food analysis
%1. Score food type/timing in boris

%read in data from camera and BORIS and find z scores and get time of food
%contact
find_food_responses();

%Get the average trial by trial responses 
get_trial_by_trial_food_times();

%run permutation test to flag responsive neurons
permutation_test_food()

% Using food times, run PCA on 15 seconds after food contact. Plot
% trajectories
plot_PCA_trajectories_food()


