%% TODO
% make it so you can graph the responses based on session number
close all
session_ind = 1; %which session to analyze
response_type = 1; %1 = all, 2 = excited, 3 = inhibited, 4 = unique
data_type = 1; %1 is z score 2 is dF/F



P = pwd;
S = dir(fullfile(P,'*.mat'));
N = {S.name};
X = ~cellfun('isempty',strfind(N,'Deliver'));
names = N(X);
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

%find trial number for each session
P = pwd;
NeurondFF_mean = [];
load(fullfile(P,names{session_ind}));
load(horzcat('z_scores_session_', num2str(session_ind), '.mat'));
load(horzcat('responsive_neurons_session_', num2str(session_ind), '.mat'));

%if an airpuff session
liquids = 1;
trial_n = [1, 2, 3];

first_neuron_avgs = [];
%for n = intersect(first_excited_neurons, first_only_neurons)
%1 = all, 2 = excited, 3 = inhibited, 4 = unique
if response_type == 1
    first_neurons = first_responding_neurons;
elseif response_type == 2
    first_neurons = first_excited_neurons;
elseif response_type == 3
    first_neurons = first_inhibited_neurons;
elseif response_type == 4
    first_neurons = first_only_neurons;
end
if data_type == 1
    data_to_use = NeuronZScores_mean;
elseif data_type == 2
    if ~isempty(NeurondFF_mean)
        data_to_use = NeurondFF_mean;
    else
        NeurondFF_mean = {};
        for row = 1:size(NeurondFF, 1)
            for col = 1:size(NeurondFF, 2)
                if ~isempty(NeurondFF{row, col})
                    NeurondFF_mean{row, col} = mean(NeurondFF{row, col});
                else
                    NeurondFF_mean{row, col} = [];
                end
            end
        end
        data_to_use = NeurondFF_mean;
    end
end
%% TODO align to delivery time
for n = first_neurons
    response = data_to_use{n, 1};
    %plot(response);
    %close all
    first_neuron_avgs = vertcat(first_neuron_avgs, response);
end

figure;
hold on
first_neuron_avg = mean(first_neuron_avgs) - mean(mean(first_neuron_avgs(:, 1:49)));
plot(first_neuron_avg, 'g');
%title(horzcat(first_liquid, ' Average Responsive Neurons'));

second_neuron_avgs = [];
%for n = intersect(second_excited_neurons, second_only_neurons)
if response_type == 1
    second_neurons = second_responding_neurons;
elseif response_type == 2
    second_neurons = second_excited_neurons;
elseif response_type == 3
    second_neurons = second_inhibited_neurons;
elseif response_type == 4
    second_neurons = second_only_neurons;
end

for n = second_neurons
    response = data_to_use{n, 2};
    second_neuron_avgs = vertcat(second_neuron_avgs, response);
end

%figure;
second_neuron_avg = mean(second_neuron_avgs) - mean(mean(second_neuron_avgs(:, 5:49)));
plot(second_neuron_avg, 'm');
%title(horzcat(second_liquid, ' Average Responsive Neurons'));

%figure;
third_neuron_avgs = [];
%for n = intersect(third_excited_neurons, third_only_neurons)
if response_type == 1
    third_neurons = third_responding_neurons;
elseif response_type == 2
    third_neurons = third_excited_neurons;
elseif response_type == 3
    third_neurons = third_inhibited_neurons;
elseif response_type == 4
    third_neurons = third_only_neurons;
end

for n = third_neurons
    response = data_to_use{n, 3};
    %plot(response);
    %close all
    
    third_neuron_avgs = vertcat(third_neuron_avgs, response);
end
if and(~isempty(third_neuron_avgs), size(third_neuron_avgs, 1))
    third_neuron_avg = mean(third_neuron_avgs) - mean(mean(third_neuron_avgs(:, 5:30)));
    plot(third_neuron_avg, 'c');
    %title(horzcat(third_liquid, ' Average Responsive Neurons'));
    legend({first_food, 'Low Fat', third_food}, 'Location','northeast');
else
    legend({first_food, second_food}, 'Location','northeast');
    third_neuron_avg = [];
end

switch response_type
    case 1
        type = 'all';
    case 2
        type = 'excited';
    case 3    
        type = 'inhibited';
end

save(horzcat(type, '_average_responses.mat'), 'third_neuron_avg', 'second_neuron_avg', 'first_neuron_avg', 'third_neuron_avgs', 'second_neuron_avgs', 'first_neuron_avgs');

savefig(horzcat('Session ', num2str(session_ind), ' ', type, '.fig'));