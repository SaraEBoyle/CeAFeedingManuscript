
%% Plot neurons present in both- remove deleted neuron data so indeces match cell reg output

%figure('Position', [100 100 400 350],'Name','Summary','numbertitle','off');
%hold on
close all
time = 0:.1:24.9;
fr_file = 'E:\Pkcd GRIN imaging\GRIN 1m\11.13.21 fr high fat sucrose palm fed3\z_scores_session_1.mat';
s_file = 'E:\Pkcd GRIN imaging\GRIN 1m\11.12.21 sated high fat sucrose fed3\z_scores_session_1.mat';

load(fr_file);
fr_NeuronZScores_aligned = NeuronZScores_aligned;
fr_neurons_to_delete = neurons_to_delete;
fr_NeuronZScores_aligned(fr_neurons_to_delete, :) = [];

load(s_file);
s_NeuronZScores_aligned = NeuronZScores_aligned;
s_neurons_to_delete = neurons_to_delete;
s_NeuronZScores_aligned(s_neurons_to_delete, :) = [];

load('responses_in_tracked_cells.mat');

for neur_ind = 1:size(indeces_in_both_sessions, 1)
    figure
    fr_ind = indeces_in_both_sessions(neur_ind, 1);
    s_ind = indeces_in_both_sessions(neur_ind, 2);
    
    fr_neuron_fat = fr_NeuronZScores_aligned{fr_ind, 1};
    s_neuron_fat = s_NeuronZScores_aligned{s_ind, 1};
    
    fr_smoothed_fat = [];
    for n = 1:size(fr_neuron_fat, 1)
        fr_smoothed_fat(n, :) = smooth(fr_neuron_fat(n, :));
    end
    
    s_smoothed_fat = [];
    for n = 1:size(s_neuron_fat, 1)
        s_smoothed_fat(n, :) = smooth(s_neuron_fat(n, :));
    end

    fr_response = fr_smoothed_fat(:, 1:250);
    fr_avg_response = mean(fr_response,1); %mean all responses
    fr_err = std(fr_response)/sqrt(size(fr_response,1)); %std error
    
    s_response = s_smoothed_fat(:, 1:250);
    s_avg_response = mean(s_response,1); %mean all responses
    s_err = std(s_response)/sqrt(size(s_response,1)); %std error

    cc = [0.2 0.8 0 %green
        0.8 0.2 0 %red
        .6 .8 1 %light blue
        1.0000 0.5000 0.8000 %pink
        0.4900 0.2800 0.5600 %purple
        0.9300 0.6900 0.1300 %orange
        0.8000 0.8000 0.8000 %gray
        0.5 0.5 0.5];
    
    subplot(2,1,1);
    AreaPlot(time,smooth(fr_avg_response,3)',smooth(fr_err)',cc(1,:),0.4,1);
    title(['Neuron # ', num2str(fr_ind)]);
    
    subplot(2,1,2); 
    AreaPlot(time,smooth(s_avg_response,3)',smooth(s_err)',cc(2,:),0.4,1);
    title(['Neuron # ', num2str(s_ind)]);
    
    
    
end






