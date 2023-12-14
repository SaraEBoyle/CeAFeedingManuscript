%% Finds the percentage of X responsive neurons that respond to all liquids 
P = pwd;
S = dir(fullfile(P,'*.mat'));
N = {S.name};
X = ~cellfun('isempty',strfind(N,'responsive_neurons_session'));
%load(N{1});

%% Fill a matrix with excitatory overlap between all stimuli
% top will be % top liquid excited neurons also excited by side liquid
name = N(X);

excitatory_response_matrix = {};
matrix_liquids = 1;
for n = 1:size(name, 2)
    title = name{n};
    first_liquid = [];
    load(title);
    if ~isempty(first_liquid)
       matrix_liquids = matrix_liquids + 1;
       %adjust top title
       excitatory_response_matrix{1, matrix_liquids} = first_liquid;
       %adjust side title
       excitatory_response_matrix{matrix_liquids, 1} = first_liquid;
       % add neurons excited by liquid
       excitatory_response_matrix{matrix_liquids, matrix_liquids} = first_excited_neurons;
       %add row/column
       matrix_liquids = matrix_liquids + 1;
       %add second liquid
       excitatory_response_matrix{1, matrix_liquids} = second_liquid;
       excitatory_response_matrix{matrix_liquids, 1} = second_liquid;
       excitatory_response_matrix{matrix_liquids, matrix_liquids} = second_excited_neurons;
       
       % add third liquid
       matrix_liquids = matrix_liquids + 1;
       excitatory_response_matrix{1, matrix_liquids} = third_liquid;
       excitatory_response_matrix{matrix_liquids, 1} = third_liquid;
       excitatory_response_matrix{matrix_liquids, matrix_liquids} = third_excited_neurons;
    else
       matrix_liquids = matrix_liquids + 1;
       %adjust top title
       excitatory_response_matrix{1, matrix_liquids} = 'Air';
       %adjust side title
       excitatory_response_matrix{matrix_liquids, 1} = 'Air';
       % add neurons excited by liquid
       excitatory_response_matrix{matrix_liquids, matrix_liquids} = air_excited_neurons;
       
    end
end

%% NOW populate the matrix
for point = 2:size(excitatory_response_matrix, 1)
    focus_liquid = excitatory_response_matrix{point, point};
    %populate the matrix
    for column = (point + 1):size(excitatory_response_matrix, 1)
        second_liquid = excitatory_response_matrix{column, column};
        %Above diagonal (top neurons|side)/side neurons
        top_given_side = intersect(focus_liquid, second_liquid);
        top_entry = size(top_given_side, 2)/size(focus_liquid, 2);
        excitatory_response_matrix{point, column} = top_entry;
        %Below diagonal (side neurons|top)/side neurons
        bottom_entry = size(top_given_side, 2)/size(focus_liquid, 2);
        excitatory_response_matrix{column, point} = bottom_entry;
        
    end
    excitatory_response_matrix{point, point} = 1;
end

%% Repeat using inhibition
inhibitory_response_matrix = {};
matrix_liquids = 1;
for n = 1:size(name, 2)
    title = name{n};
    first_liquid = [];
    load(title);
    if ~isempty(first_liquid)
       matrix_liquids = matrix_liquids + 1;
       %adjust top title
       inhibitory_response_matrix{1, matrix_liquids} = first_liquid;
       %adjust side title
       inhibitory_response_matrix{matrix_liquids, 1} = first_liquid;
       % add neurons excited by liquid
       inhibitory_response_matrix{matrix_liquids, matrix_liquids} = first_inhibited_neurons;
       %add row/column
       matrix_liquids = matrix_liquids + 1;
       %add second liquid
       inhibitory_response_matrix{1, matrix_liquids} = second_liquid;
       inhibitory_response_matrix{matrix_liquids, 1} = second_liquid;
       inhibitory_response_matrix{matrix_liquids, matrix_liquids} = second_inhibited_neurons;
       
       % add third liquid
       matrix_liquids = matrix_liquids + 1;
       inhibitory_response_matrix{1, matrix_liquids} = third_liquid;
       inhibitory_response_matrix{matrix_liquids, 1} = third_liquid;
       inhibitory_response_matrix{matrix_liquids, matrix_liquids} = third_inhibited_neurons;
    else
       matrix_liquids = matrix_liquids + 1;
       %adjust top title
       inhibitory_response_matrix{1, matrix_liquids} = 'Air';
       %adjust side title
       inhibitory_response_matrix{matrix_liquids, 1} = 'Air';
       % add neurons excited by liquid
       inhibitory_response_matrix{matrix_liquids, matrix_liquids} = air_inhibited_neurons;
       
    end
end

%% NOW populate the matrix
for point = 2:size(inhibitory_response_matrix, 1)
    focus_liquid = inhibitory_response_matrix{point, point};
    %populate the matrix
    for column = (point + 1):size(inhibitory_response_matrix, 1)
        second_liquid = inhibitory_response_matrix{column, column};
        %Above diagonal (top neurons|side)/side neurons
        top_given_side = intersect(focus_liquid, second_liquid);
        top_entry = size(top_given_side, 2)/size(focus_liquid, 2);
        inhibitory_response_matrix{point, column} = top_entry;
        %Below diagonal (side neurons|top)/side neurons
        bottom_entry = size(top_given_side, 2)/size(focus_liquid, 2);
        inhibitory_response_matrix{column, point} = bottom_entry;
        
    end
    inhibitory_response_matrix{point, point} = 1;
end

%{
%% Repeat to check excitation/inhibition overlap
exinhibitory_response_matrix = {};
matrix_liquids = 1;
for n = 1:size(name, 2)
    title = name{n};
    first_liquid = [];
    load(title);
    if ~isempty(first_liquid)
       matrix_liquids = matrix_liquids + 1;
       %adjust top title
       exinhibitory_response_matrix{1, matrix_liquids} = first_liquid;
       %adjust side title
       exinhibitory_response_matrix{matrix_liquids, 1} = first_liquid;
       % add neurons excited by liquid
       exinhibitory_response_matrix{matrix_liquids, matrix_liquids} = first_inhibited_neurons;
       %add row/column
       matrix_liquids = matrix_liquids + 1;
       %add second liquid
       exinhibitory_response_matrix{1, matrix_liquids} = second_liquid;
       exinhibitory_response_matrix{matrix_liquids, 1} = second_liquid;
       exinhibitory_response_matrix{matrix_liquids, matrix_liquids} = second_inhibited_neurons;
       
       % add third liquid
       matrix_liquids = matrix_liquids + 1;
       exinhibitory_response_matrix{1, matrix_liquids} = third_liquid;
       exinhibitory_response_matrix{matrix_liquids, 1} = third_liquid;
       exinhibitory_response_matrix{matrix_liquids, matrix_liquids} = third_inhibited_neurons;
    else
       matrix_liquids = matrix_liquids + 1;
       %adjust top title
       exinhibitory_response_matrix{1, matrix_liquids} = 'Air';
       %adjust side title
       exinhibitory_response_matrix{matrix_liquids, 1} = 'Air';
       % add neurons excited by liquid
       exinhibitory_response_matrix{matrix_liquids, matrix_liquids} = air_inhibited_neurons;
       
    end
end

%% NOW populate the matrix
for point = 2:size(exinhibitory_response_matrix, 1)
    focus_liquid = exinhibitory_response_matrix{point, point};
    %populate the matrix
    for column = (point + 1):size(exinhibitory_response_matrix, 1)
        second_liquid = exinhibitory_response_matrix{column, column};
        %Above diagonal (top neurons|side)/side neurons
        top_given_side = intersect(focus_liquid, second_liquid);
        top_entry = size(top_given_side, 2)/size(focus_liquid, 2);
        exinhibitory_response_matrix{point, column} = top_entry;
        %Below diagonal (side neurons|top)/side neurons
        bottom_entry = size(top_given_side, 2)/size(focus_liquid, 2);
        exinhibitory_response_matrix{column, point} = bottom_entry;
        
    end
    exinhibitory_response_matrix{point, point} = 1;
end
%}
save('overlap_response_matrices.mat', 'inhibitory_response_matrix', 'excitatory_response_matrix')
halp = [];





