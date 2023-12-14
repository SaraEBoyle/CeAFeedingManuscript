%% combine all sets of single neuron response matrices from different mice into a super set
clear variables
P = pwd;
S = dir(fullfile(P,'*.mat'));
N = {S.name};
X = ~cellfun('isempty',strfind(N,'matrices'));
files_to_combine = N(X);

full_avg_avg_responses = {};
full_avg_responses = {};
full_responses = {};

for file_ind = 1:size(files_to_combine, 2)
    load(fullfile(P,files_to_combine{file_ind}));
    if file_ind == 1
        full_avg_avg_responses = avg_avg_response_matrix;
        full_avg_responses = avg_response_matrix;
        full_responses = full_response_matrix;
    else
        % add size of new array minus 1 taking into account headers
        full_rows = size(full_avg_avg_responses, 1);
        new_rows = size(avg_avg_response_matrix, 1) - 1;
        full_avg_avg_responses((full_rows + 1):(full_rows + new_rows), :) = avg_avg_response_matrix(2:end, :);
        full_avg_responses((full_rows + 1):(full_rows + new_rows), :) = avg_response_matrix(2:end, :);
        full_responses((full_rows + 1):(full_rows + new_rows), :) = full_response_matrix(2:end, :);
    end
end



save('all_mouse_liquid_air_responses.mat', 'full_avg_avg_responses', 'full_avg_responses', ...
    'full_responses', 'air_response_window', 'liquid_response_window');