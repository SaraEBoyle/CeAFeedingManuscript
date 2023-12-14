%% Calculate z score and screen neurons for quality
% show the neuron shape
% Author: Sara Boyle 2022

clear,clc
close all
%addpath(genpath('Z:\Work\WorkCode'))
%% extract file info
new_screen_neurons = 1;
calculate_z_score_from_all_trial_types = 1;
load_old_screen = 0;
num_to_plot = 3; % how many pellet types to plot
baseline_time = 3; % first X s (duration of baseline)
target_frame = 60;
trial_length = 25;

nam = pwd;
[filepath,file_name,ext] = fileparts(nam);
%%
%load('FrameID.mat')
load('neuron.mat')
signal_id = 2; % 1, C; 2, C_raw; 3, S
P = pwd;

%% Parameters

R = 28; % for colorspace
C = 28; % for colorspace
ColorMap = zeros(64,3);
ColorMap(64-R+1:64,1) = linspace(0,1,R); % R
ColorMap(1:C,2) = linspace(1,0,C); % G
ColorMap(1:C,3) = linspace(1,0,C); % B

if ~load_old_screen

psth_folder_name = 'psth_neuron_session';
%%
% C_raw correspond to a scaled version of DF, which is a metric used in most calcium imaging literature. 
% C is its denoised version, and S is the inferred spiking activity.
% signal_id = 2; % 1, C; 2, C_raw; 3, S
C = neuron.C;
C_raw = neuron.C_raw; % scaled dF/F
S = neuron.S; % ubferred spiking activity
A = neuron.A; % shape of neuron


frame_rate = 10;
frame_length = trial_length * frame_rate;
frames = size(C,2);
trial_num_recorded = frames/frame_length;

%% Get some basic info
disp('Analyze baseline only');

%% Get times of motor movement from bpod

actual_frames_recorded = size(neuron.C,2);

disp(['Total imaging frames: ', num2str(actual_frames_recorded)])
disp(['Neuron number (DAQ): ', num2str(size(C,1))])

switch signal_id
    case 1
        SignalD = C;
    case 2
        SignalD = C_raw;
    case 3
        SignalD = S;
end
%%

frames = size(SignalD,2);
trial_IDs = zeros(frames, 1);

for p = 1:trial_num_recorded
    trial_IDs(((p - 1)*frame_length + 1):(p*frame_length)) = p;
end

if size(trial_IDs, 1) ~= size(SignalD, 2)
    disp('Size of frames recorded is inconsistent with the number of trials.')
        if size(trial_IDs, 1) > size(SignalD, 2)
            % SignalID is created with predicted frame size
        end
end

%% so here the order of things is wrong because there are no sucrose trials
%to make the code resilient, take into account all possible trial labels 

NeuronS = []; %Wrong thing, all the same
for neuron_id = 1:size(SignalD,1)
    signal_df = SignalD(neuron_id,:); % each neuron, whole session   
    for tt=1
        NeuronS{neuron_id,tt} = []; %Save avg response from neuron from each trial type
        %if sum(trial_n(tt))>0 %if not an image trial
            for k=1:trial_num_recorded %for each trial
                %Get all frames of the right trial number (sel)
                if tt == 1
                    sel = ((trial_IDs==k));
                    NeuronS{neuron_id,tt}{k,1} = [signal_df(sel');(.1:.1:trial_length)];
                end
            end
        %end
    end
end

%%
NeurondFF = [];
NeuronZScores = [];
NeuronZScores_mean = [];
total_stds = {};
total_means = {};
for neuron_id = 1:size(NeuronS,1)
    %% delete the empty cells to get real trial info
    for ii=1:1
        signal_m = [];
        if isempty(NeuronS{neuron_id,ii})
            continue;
        end
        to_delete = [];
        for k=1:numel(NeuronS{neuron_id,ii})
            temp_d = NeuronS{neuron_id,ii}{k};
            if isempty(temp_d)
                to_delete = horzcat(to_delete, k);
                continue
            end
            signal_m(k,:) = temp_d(1,:);
            time_x = temp_d(2,:);
        end
        signal_m(to_delete, :) = [];
        base_sel = and(((1:length(time_x)) < baseline_time*10), ((1:length(time_x)) > 2)); %Exclude first half second (bleaching)
        base_mean0 = mean(signal_m(:,base_sel),2);
        std_fodder = signal_m(:,base_sel);
        base_mean = mean(base_mean0);
        base_std = mean(std(signal_m(:,base_sel))');
        total_stds{neuron_id, ii} = std(signal_m(:,base_sel)');
        total_raw_baselines{neuron_id, ii} = std_fodder;
        total_means{neuron_id, ii} = base_mean0;
        NeuronZScores{neuron_id,ii} = (signal_m - base_mean)./base_std;
        NeuronZScores_mean{neuron_id,ii} = (mean(signal_m,1) - base_mean)./base_std;
        NeurondFF{neuron_id,ii} = signal_m;
        mean_NeurondFF{neuron_id,ii} = mean(signal_m, 1);
    end
end

%% Calculate z scores based on entire session
if calculate_z_score_from_all_trial_types
    for neuron_id = 1:size(NeuronS,1)
        %calculate std and mean baseline for each neuron
        all_baselines = total_raw_baselines(neuron_id, :);
        first = reshape(all_baselines{1,1}.',1,[]);
        all_baselines = first;
        neuron_std = std(all_baselines);
        neuron_mean = mean(all_baselines);
        for type = 1:1
            signal_m = NeurondFF{neuron_id,type};
            NeuronZScores{neuron_id,ii} = (signal_m - neuron_mean)./neuron_std;
            NeuronZScores_mean{neuron_id,ii} = (mean(signal_m,1) - neuron_mean)./neuron_std;
        end
    end
end


%% TODO calculate z score with NeuronSum and total_stds, total_means
individual_max = [];
individual_min = [];
for neuron_num = 1:size(total_means, 1)
    %Find std of the neuron in all trial types
    full_std = [];
    for p = 1:1
        current_std = total_stds{neuron_num, p};
        full_std = horzcat(full_std, current_std);
    end
    individual_std = mean(full_std);
    
    full_means = [];
    for p = 1:1
        current_mean = total_means{neuron_num, p};
        full_means = horzcat(full_means, current_mean');
    end
    individual_mean = mean(full_means);
    
    for trial_cat = 1:1
        raw_data = NeurondFF{neuron_num,trial_cat};
        if isempty(raw_data)
            continue
        end
        z_score = (raw_data - individual_mean)./individual_std;
        %NeuronZScores{neuron_num,trial_cat} = z_score;
        %NeuronZScores_mean{neuron_num,trial_cat} = (mean(z_score,1) - individual_mean)./individual_std; 
    end
    
    full_maxes = [];
    for p = 1:1
        raw_data = NeuronZScores{neuron_num,p};
        if isempty(raw_data)
            continue
        end
        current_max = max(mean(raw_data, 1));
        full_maxes = horzcat(full_maxes, current_max);
    end
    individual_max(neuron_num) = max(full_maxes);
    
    full_mins = [];
    for p = 1:1
        raw_data = NeuronZScores{neuron_num,p};
        if isempty(raw_data)
            continue
        end
        current_min = min(mean(raw_data, 1));
        full_mins = horzcat(full_mins, current_min);
    end
    individual_min(neuron_num) = min(full_mins);
    
end
%%
% ColorMap = exp_colormap('red-green',64);




%% Plot the individual neurons
if ~exist(psth_folder_name, 'dir')
    mkdir(psth_folder_name);
end


%% for display purposes, align to when pellet reaches mouse
%NeuronZScores
adjusted_axis = time_x;

% get the frame closest to 0
[zero_times, zero_inds] = min(abs(adjusted_axis)');

% Align heatmaps to time 0
NeuronZScores_aligned = {};
NeuronZScores_baseline_normed = {};
% for each trial type
for type_ind = 1:size(NeuronZScores, 2)
    %for each neuron
    for neur = 1:size(NeuronZScores, 1)
        neuron_type_response = NeuronZScores{neur, type_ind};
        % for each trial
        full_neuron = [];
        normed_neuron = [];
        for trial = 1:size(neuron_type_response, 1)
            % pad frames on front/back to get an aligned heat map
            single_trial = neuron_type_response(trial, :);
            single_trial = single_trial - mean(single_trial(3:30));
            normed_neuron = [normed_neuron;single_trial];
            full_neuron = [full_neuron;single_trial];
        end
        NeuronZScores_aligned{neur, type_ind} = full_neuron;
        NeuronZScores_baseline_normed{neur, type_ind} = normed_neuron;
        mean_NeuronZScores_aligned{neur, type_ind} = mean(full_neuron);
    end
end

if new_screen_neurons == 1
    neurons_to_delete = [];
    for neuron_id = 1:size(NeuronZScores,1)   
    
        f = figure('Position', [100 100 1500 800],'Name','Summary','numbertitle','off');
        f.Renderer = 'painters';
        for trial_i=1:1
            %sm = NeuronZScores_aligned{neuron_id,trial_i};
            
            if trial_i > size(NeuronZScores_aligned, 2)
                continue;
            end
            
            sm = NeuronZScores_aligned{neuron_id,trial_i};
            if isempty(sm)
                continue;
            end
        
            Y = mean(sm,1);
            Y_err = std(sm)/sqrt(size(sm,1));        
            h(trial_i) = subplot(4,1,trial_i +1);
            %AreaPlot(time_x-time_x(1),Y,Y_err,'k',0.4,1); changed because I
            %dont have that function
            plot(time_x-time_x(1),Y);
            ylim([individual_min(neuron_id), individual_max(neuron_id)]);
            hold on
            
            xlabel('Time (s)','FontSize', 8)
            ylabel('z-score(dF)','FontSize', 8)
            set(gca,'TickDir', 'out','xlim',[0,time_x(end)-time_x(1)],'xtick',0:2:time_x(end)-time_x(1),'FontSize', 8,'box','off');
            title(['Neuron',num2str(neuron_id)],'FontSize', 6)

            ax1 = subplot(4,num_to_plot,trial_i);
            imagesc(time_x-time_x(1),1:size(sm,1),sm);
            colormap(ax1,ColorMap);
            colorbar;
            % xlabel('Time (s)','FontSize', 8)
            ylabel('Trial #','FontSize', 8)
            set(gca,'TickDir', 'out','xlim',[0,time_x(end)-time_x(1)],'xtick',0:2:time_x(end)-time_x(1),'YDir','reverse','FontSize', 8,'box','off');        
            caxis([-6 6])      
        end
        linkaxes(h,'xy');
    
        subplot(4,3,[7,8,9])
        neuron_trace = neuron.C_raw(neuron_id, :);
        %plot((1:length(neuron_trace))/10, neuron.C(neuron_id, :)*max(neuron.A(:, neuron_id)),'k');
        
        plot((1:length(neuron_trace))/10,neuron_trace);
        % neuron.viewNeurons(neuron_id);
        % neuron_shape = neuron.reshape(neuron.A(:, neuron_id), 2);
        xlabel(['Neuron ID: ',num2str(neuron_id)],'FontSize', 16)
        set(gca,'TickDir', 'out','FontSize', 8,'box','off');

        ax2 = subplot(4,3,[10]);
        Amask = (neuron.A~=0);
        neuron.image(neuron.A(:, neuron_id).*Amask(:, neuron_id)); %
        colormap(ax2,parula)
        axis equal; axis off;

        ax3 = subplot(4,3,[12]);
        ctr = neuron.estCenter();      %neuron's center
        gSiz = neuron.options.gSiz;        % maximum size of a neuron
        neuron.image(neuron.A(:, neuron_id).*Amask(:, neuron_id)); %
        colormap(ax3,parula)
        axis equal; axis off;
        x0 = ctr(neuron_id, 2);
        y0 = ctr(neuron_id, 1);
        if ~isnan(x0)
            xlim(x0+[-gSiz, gSiz]*2);
            ylim(y0+[-gSiz, gSiz]*2);
        end
    
    
        keep = input('Keep this neuron? Enter for yes, n for no, q for quit','s');

        if strcmp(keep, 'n') %if 1 == sum(ismember(old_to_delete, neuron_id))%
            neurons_to_delete = horzcat(neurons_to_delete, neuron_id);
            if ~exist([psth_folder_name, '\deleted'], 'dir')
                mkdir([psth_folder_name, '\deleted']);
            end
            print([pwd,'\',psth_folder_name,'\deleted\NeuronPSTH_',num2str(neuron_id)],'-dpng');
        elseif strcmp(keep, 'q')
            return
        else
            print([pwd,'\',psth_folder_name,'\NeuronPSTH_',num2str(neuron_id)],'-dpng');
            savefig([pwd,'\',psth_folder_name,'\NeuronPSTH_',num2str(neuron_id)]);
        end
        close all
    end
    
    close;
    save('neurons_to_delete.mat','neurons_to_delete')
end

%%
% ColorMap = exp_colormap('red-green',64);
else
    load('z_scores_session_1.mat');
    load('neurons_to_delete.mat');
end
A = {};
for k=1:size(NeuronZScores,2)
    %if ~ismember(neurons_to_delete, k)
        A0 = NeurondFF(:,k);
        A0(neurons_to_delete) = [];
        A{k} = cell2mat(A0);
    %end
end
%Sort by 3 -2, fat


first_sort_ind = 1;
time_x = .1:.1:trial_length;
a_sort0 = A{first_sort_ind}(:,time_x>=5 & time_x<=(220));
a_sort = mean(a_sort0,2);
[~,idx]=sortrows(a_sort);

figure('Position', [100 100 1200 400],'Name','Summary','numbertitle','off');
for k=(1):size(NeuronZScores,2)
    if isempty(A{k})
        continue
    end
    A{k} = A{k}(idx,:);
    
    subplot(2,size(NeuronZScores_aligned,2),k);
    imagesc(time_x-time_x(1),1:size(A{k},1),A{k});
    colormap(ColorMap);
    colorbar;
    xlabel('Time (s)','FontSize', 18)
    ylabel('Neuron #','FontSize', 18)
    caxis([-50 50]);
    set(gca,'TickDir', 'out',...
        'YDir','reverse','FontSize', 16,'box','off');
    
    subplot(2,size(NeuronZScores_aligned,2),k +size(NeuronZScores_aligned,2))
    plot(time_x-time_x(1),mean(A{k},1),'k','LineWidth',2)
    xlabel('Time (s)','FontSize', 18)
    ylabel('Norm activity','FontSize', 18)
    ylim([-1 3])
    xlim([0,time_x(end)-time_x(1)]);
    set(gca,'TickDir', 'out','FontSize', 16,'box','off');
    
end 

if 0
    NeuronZScores(neurons_to_delete, :) = [];
    NeurondFF(neurons_to_delete, :) = [];
    NeuronZScores_mean(neurons_to_delete, :) = [];
end
save(horzcat('z_scores_session_', num2str(1), '.mat'), 'neurons_to_delete', ...
    'NeuronZScores', 'NeurondFF', 'NeuronZScores_mean', 'NeuronZScores_aligned', ...
    'mean_NeuronZScores_aligned', 'NeuronZScores_baseline_normed', 'target_frame');

    print(gcf,'Heatmap_sorted','-dpng');
    saveas(gcf,'Heatmap_sorted.fig');


