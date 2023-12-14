%% Makes a response matrix where each row is a different neuron, and 
%% each column is an average response to each stimulus.

liquid_response_window = [41:100]; %4.1 seconds to 14
air_response_window = [41:70]; %4.1 seconds to 14

P = pwd;
S = dir(fullfile(P,'*.mat'));
N = {S.name};
X = ~cellfun('isempty',strfind(N,'z_scores_session'));
z_name = N(X);

S = dir(fullfile(P,'*.mat'));
N = {S.name};
X = ~cellfun('isempty',strfind(N,'responsive_neurons_session'));
r_name = N(X);


all_neurons_to_delete = [];
for z = 1:size(z_name, 2)
    load(z_name{z});
    all_neurons_to_delete = [all_neurons_to_delete, neurons_to_delete];
end
all_neurons_to_delete = unique(all_neurons_to_delete);


%% Fill a matrix with excitatory overlap between all stimuli
% top will be % top liquid excited neurons also excited by side liquid


response_matrix = {};
matrix_liquids = 0;
for n = 1:size(z_name, 2)
    z_title = z_name{n};
    r_title = r_name{n};
    
    first_liquid = [];
    load(z_title);
    load(r_title);
    if ~isempty(first_liquid)
       %add first liquid column
       matrix_liquids = matrix_liquids + 1;
       response_matrix{1, matrix_liquids} = first_liquid;
       
       %add second liquid column
       matrix_liquids = matrix_liquids + 1;
       response_matrix{1, matrix_liquids} = second_liquid;
       
       %add third liquid column
       matrix_liquids = matrix_liquids + 1;
       response_matrix{1, matrix_liquids} = third_liquid;
       response_window = liquid_response_window;
       
       if exist('fourth_liquid', 'var')
           %add fourth liquid column
           matrix_liquids = matrix_liquids + 1;
           response_matrix{1, matrix_liquids} = fourth_liquid;
       end
       
       if exist('fifth_liquid', 'var')
           %add fourth liquid column
           matrix_liquids = matrix_liquids + 1;
           response_matrix{1, matrix_liquids} = fifth_liquid;
       end
       
       if exist('sixth_liquid', 'var')
           %add fourth liquid column
           matrix_liquids = matrix_liquids + 1;
           response_matrix{1, matrix_liquids} = sixth_liquid;
       end
       
    else
       matrix_liquids = matrix_liquids + 1;
       %adjust top title
       response_matrix{1, matrix_liquids} = 'Air';
       response_window = air_response_window;
    end
    % now populate the matrix
    NeuronZScores_mean(all_neurons_to_delete, :) = [];
    z_means = NeuronZScores_mean;
    if isempty(z_means{1,1})
        z_means(:,1) = [];
    end
    for row = 1:size(z_means, 1)
        for col = 1:size(z_means, 2)
            current_res = z_means{row,col};
            avg = mean(current_res(response_window));
            z_means{row,col} = avg;
        end
    end
    if ~isempty(first_liquid)
        % original
        %response_matrix(2:(1 + size(z_means, 1)), matrix_liquids - 2) = z_means(:, 1);
        %response_matrix(2:(1 + size(z_means, 1)), matrix_liquids - 1) = z_means(:, 2);
        %response_matrix(2:(1 + size(z_means, 1)), matrix_liquids) = z_means(:, 3);
        %% Changed 5/12/22 to handle 6 liquids instead of 3
        for liq = 1:matrix_liquids
            response_matrix(2:(1 + size(z_means, 1)), liq) = z_means(:, liq);
        end
    else
        response_matrix(2:(1 + size(z_means, 1)), matrix_liquids) = z_means(:, 1);
    end
end

%% Sort the matrix based on fat response
fats = cell2mat(response_matrix(2:end, 1));
[~, inds] = sort(fats);
response_matrix_no_titles = response_matrix(2:end, :);
sorted_response_matrix_no_titles = response_matrix_no_titles(inds, :);
response_matrix(2:(1 + size(z_means, 1)), 1:matrix_liquids) = sorted_response_matrix_no_titles;

save('single_neuron_avg_response_matrices.mat', 'response_matrix', 'liquid_response_window', 'air_response_window')
halp = [];





