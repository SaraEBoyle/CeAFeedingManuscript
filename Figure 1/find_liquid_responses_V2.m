%% Analyzing calcium activities related to behavioral events
% Sync GRIN imaging brain responses to a panel of 6 liquids + airpuff
% stimuli
% Author: Sara Boyle

close all

%% extract file info
screen_neurons = 1; % change to 1 if you haven't screened yet, 0 if you want to load your old screen
nam = pwd;
[filepath,file_name,ext] = fileparts(nam);
%%
%load('FrameID.mat')
load('neuron.mat')
signal_id = 2; % 1, C; 2, C_raw; 3, S
P = pwd;
S = dir(fullfile(P,'*.mat'));
N = {S.name};
X = ~cellfun('isempty',strfind(N,'Deliver'));
names = N(X);
session_ind = 1; %which session to analyze
trials_with_no_licks = [];
baseline_time = 4; % first X s (duration of baseline)

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
total_trial_num = [];
for h = 1:size(names, 2)
    %load each session
    load(names{h});
    if isfield(SessionData, 'TrialTypes')
        %only liquid trials have this
        trials = SessionData.TrialTypes;
        trial_n = unique(SessionData.TrialTypes);

        %Get array of trial numbers
        TrialTypes = SessionData.TrialTypes;
        bad_inds = find(TrialTypes == -1);
        TrialTypes(bad_inds) = [];
        total_trial_num(h) = size(TrialTypes, 2);
    else
        %if air session
        total_trial_num(h) = SessionData.nTrials;
    end
    
end


load(fullfile(P,names{session_ind}));

if and(isfield(SessionData.TrialSettings, 'first_liquid'), session_ind == 1)
    %if 1st session
    if exist('trials_with_no_licks_first_session.mat', 'file')
        load(horzcat('trials_with_no_licks_first_session.mat'));
    else
        load(horzcat('trials_with_no_licks.mat'));
    end
elseif and(isfield(SessionData.TrialSettings, 'first_liquid'), session_ind == 2)
    load(horzcat('trials_with_no_licks_second_session.mat'));
elseif ~isfield(SessionData.TrialSettings, 'first_liquid')
    trials_with_no_licks = [];
end

%% Parameters


if isfield(SessionData.TrialSettings, 'first_liquid')
    first_liquid = SessionData.TrialSettings.first_liquid;
    second_liquid = SessionData.TrialSettings.second_liquid;
    third_liquid = SessionData.TrialSettings.third_liquid;
    fourth_liquid = SessionData.TrialSettings.fourth_liquid;
    fifth_liquid = SessionData.TrialSettings.fifth_liquid;
    sixth_liquid = SessionData.TrialSettings.sixth_liquid;
    title_text = {first_liquid, second_liquid, third_liquid, fourth_liquid, fifth_liquid, sixth_liquid, 'Airpuff Trials'};
else
    title_text = {'Airpuff Trials'};
end
stim_on = baseline_time+0; % timestamp of CS (or other interested time points) (only for sorting neurons)
stim_off = stim_on+4; % timestamp of CS-off (only for sorting neurons)
X_line = [stim_on, stim_off]; % reference line for CS and US



psth_folder_name = horzcat('psth_neuron_session_', num2str(session_ind));
%%
% C_raw correspond to a scaled version of DF, which is a metric used in most calcium imaging literature. 
% C is its denoised version, and S is the inferred spiking activity.
% signal_id = 2; % 1, C; 2, C_raw; 3, S
C = neuron.C;
C_raw = neuron.C_raw; % scaled dF/F
S = neuron.S; % ubferred spiking activity
A = neuron.A; % shape of neuron

trial_length = 14;
frame_rate = 10;
frame_length = trial_length * frame_rate;
frames = size(C,2);
trial_num_recorded = frames/frame_length;

%% Get some basic info
if isfield(SessionData.TrialSettings, 'first_liquid')
    trials = SessionData.TrialTypes;
    trial_n = unique(SessionData.TrialTypes);

    if ismember(-1, trial_n)
        pallet_cleansers = 1; % if 1, you included trials with no imaging
    else
        pallet_cleansers = 0;
    end
    disp(horzcat('1: ', first_liquid));
    disp(horzcat('2: ', second_liquid));
    disp(horzcat('3: ', third_liquid));
    disp(horzcat('4: ', fourth_liquid));
    disp(horzcat('5: ', fifth_liquid));
    disp(horzcat('6: ', sixth_liquid));
    if ismember(9, trial_n)
        disp(horzcat('9: ', 'Airpuff'));
    end

    %Get array of trial numbers

    TrialTypes = SessionData.TrialTypes;
    bad_inds = find(TrialTypes == -1);
    TrialTypes(bad_inds) = [];
    trial_num_bpod = size(TrialTypes, 2);
else
    % for airpuff session
    pallet_cleansers = 0;
    trial_num_bpod = SessionData.nTrials;
    TrialTypes = ones(1, trial_num_bpod);
end
actual_frames = total_trial_num(session_ind)*trial_length*frame_rate;

disp(['Total imaging frames: ', num2str(size(neuron.C,2))])
disp(['Imaging frames to be analyzed: ', ... 
    num2str(actual_frames)])
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
R = 28; % for colorspace
C = 28; % for colorspace
ColorMap = zeros(64,3);
ColorMap(64-R+1:64,1) = linspace(0,1,R); % R
ColorMap(1:C,2) = linspace(1,0,C); % G
ColorMap(1:C,3) = linspace(1,0,C); % B

%%




%Delete the imaging frames not in current session
if trial_num_bpod ~= trial_num_recorded
    disp('Mutiple sessions read.');
    disp(horzcat('Truncating data to ', names{session_ind}));
    if session_ind == 1
        disp('Analyzing first session... ');
        actual_trials = total_trial_num(1);
        SignalD(:, (actual_frames + 1):end) = [];
    elseif session_ind == 2
        disp('Analyzing second session... ');
        actual_trials = total_trial_num(2);
        prev_frames = total_trial_num(1)*trial_length*frame_rate;
        
        SignalD(:, 1:(prev_frames)) = [];
        %if there's a third session, delete that too
        if actual_frames < size(SignalD,2)
            disp('Found third session. Deleting frames...')
            SignalD(:, (actual_frames + 1):end) = [];
        end
        if ~isfield(SessionData.TrialSettings, 'first_liquid')
            %if air session
            trial_n = 1;
        end
    elseif session_ind == 3
        disp('Analyzing third air puff session');
        actual_trials = total_trial_num(3);
        prev_frames = total_trial_num(1)*trial_length*frame_rate + total_trial_num(2)*trial_length*frame_rate;
        SignalD(:, 1:(prev_frames)) = [];
        trial_n = 1;
        TrialTypes = ones(SessionData.nTrials, 1);
        pallet_cleansers = 0;
    end
    %total_trial_num(session_ind);
else
    disp('Analyzing only session... ');
    actual_trials = total_trial_num;
end

frames = size(SignalD,2);
trial_IDs = zeros(frames, 1);

for p = 1:actual_trials
    trial_IDs(((p - 1)*frame_length + 1):(p*frame_length)) = p;
end

NeuronS = []; %Wrong thing, all the same
for neuron_id = 1:size(SignalD,1)
    signal_df = SignalD(neuron_id,:); % each neuron, whole session   
    for tt=1:numel(trial_n) %For each type of trial
        NeuronS{neuron_id,tt} = []; %Save avg response from neuron from each trial type
        if sum(trial_n(tt))>0 %if not an image trial
            for k=1:actual_trials %for each trial
                if isempty(find(trials_with_no_licks == k))
                %Get all frames of the right trial number (sel)
                if trial_n(tt) == TrialTypes(k)
                    sel = ((trial_IDs==k));
                    NeuronS{neuron_id,tt}{k,1} = [signal_df(sel');(.1:.1:14)];
                end
                else
                    
                end
            end
        end
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
    for ii=1:size(NeuronS,2)
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
        base_sel = and(((1:length(time_x)) < baseline_time*10), ((1:length(time_x)) > 5)); %Exclude first half second (bleaching)
        base_mean0 = mean(signal_m(:,base_sel),2);
        base_mean = mean(base_mean0);
        base_std = mean(std(signal_m(:,base_sel))');
        total_stds{neuron_id, ii} = std(signal_m(:,base_sel)');
        total_means{neuron_id, ii} = base_mean0;
        NeuronZScores{neuron_id,ii} = (signal_m - base_mean)./base_std;
        NeuronZScores_mean{neuron_id,ii} = (mean(signal_m,1) - base_mean)./base_std;
        NeurondFF{neuron_id,ii} = signal_m;     
    end
end




%% for display purposes, align to when pellet reaches mouse
%NeuronZScores
%{'Fat 5%','Sucrose 12.5%','Xanthan Gum','Quinine','Water','Mineral Oil'}
trial_inds = TrialTypes;
fat_inds = find(trial_inds == 1);
sucrose_inds = find(trial_inds == 2);
XG_inds = find(trial_inds == 3);
quinine_inds = find(trial_inds == 4);
water_inds = find(trial_inds == 5);
MO_inds = find(trial_inds == 6);

load('lick_data_first_session.mat');


fat_empties = find(cell2mat(cellfun(@(x) isempty(x), first_offset, 'UniformOutput',false)) == 1);
sucrose_empties = find(cell2mat(cellfun(@(x) isempty(x), second_offset, 'UniformOutput',false)) == 1);
XG_empties = find(cell2mat(cellfun(@(x) isempty(x), third_offset, 'UniformOutput',false)) == 1);
quinine_empties = find(cell2mat(cellfun(@(x) isempty(x), fourth_offset, 'UniformOutput',false)) == 1);
water_empties = find(cell2mat(cellfun(@(x) isempty(x), fifth_offset, 'UniformOutput',false)) == 1);
MO_empties = find(cell2mat(cellfun(@(x) isempty(x), sixth_offset, 'UniformOutput',false)) == 1);

fat_frames = cellfun(@(x) round(x*10), first_offset, 'UniformOutput',false);
sucrose_frames = cellfun(@(x) round(x*10), second_offset, 'UniformOutput',false);
XG_frames = cellfun(@(x) round(x*10), third_offset, 'UniformOutput',false);
quinine_frames = cellfun(@(x) round(x*10), fourth_offset, 'UniformOutput',false);
water_frames = cellfun(@(x) round(x*10), fifth_offset, 'UniformOutput',false);
MO_frames = cellfun(@(x) round(x*10), sixth_offset, 'UniformOutput',false);

fat_frames(fat_empties) = []; 
sucrose_frames(sucrose_empties) = [];
XG_frames(XG_empties) = [];
quinine_frames(quinine_empties) = [];
water_frames(water_empties) = [];
MO_frames(MO_empties) = [];

target_time = 4;
target_frame = 40;
adjusted_axis = time_x - target_time;

% get the frame closest to 0
[zero_times, zero_inds] = min(abs(adjusted_axis)');
zeros{1}= fat_frames;
zeros{2}= sucrose_frames;
zeros{3}= XG_frames;
zeros{4}= quinine_frames;
zeros{5}= water_frames;
zeros{6}= MO_frames;
zeros{7}= cell(1, 50);

% Align heatmaps to time 0
NeuronZScores_aligned = {};
NeuronZScores_baseline_normed = {};
% for each trial type
for type_ind = 1:(size(NeuronZScores, 2))
    %for each neuron
    if type_ind == 1
        continue
    end
    if type_ind ~= 8
        food_frames = zeros{type_ind - 1};
    end
    for neur = 1:size(NeuronZScores, 1)
        neuron_type_response = NeuronZScores{neur, type_ind};
        %align_frames = zeros{type_ind};
        % for each trial
        full_neuron = [];
        normed_neuron = [];
        for trial = 1:size(neuron_type_response, 1)
            % pad frames on front/back to get an aligned heat map
            single_trial = neuron_type_response(trial, :);
            single_trial = single_trial - mean(single_trial(3:30));
            normed_neuron = [normed_neuron;single_trial];
            if type_ind ~= 8
                food_frame = food_frames{trial};
            else
                food_frame = target_frame;
            end
            frame_dif = food_frame - target_frame;
            if frame_dif > 0
            % if positive, shift single trial left by removing frames from
            % beginning and tacking them on at end, just duplicating the
            % difference
                single_trial((end + 1):(end + frame_dif)) = single_trial((end - frame_dif + 1):end);
                single_trial(1:frame_dif) = [];
            elseif frame_dif < 0
            % if negative, shift single trial right by adding on to beggining, then removing from end
                to_add = single_trial(1:abs(frame_dif));
                single_trial((end + frame_dif + 1):end) = [];
                single_trial = [single_trial, to_add];
            end
            full_neuron = [full_neuron;single_trial];
        end
        NeuronZScores_aligned{neur, type_ind} = full_neuron;
        NeuronZScores_baseline_normed{neur, type_ind} = normed_neuron;
        mean_NeuronZScores_aligned{neur, type_ind} = mean(full_neuron);
    end
end













%% TODO calculate z score with NeuronSum and total_stds, total_means
individual_max = [];
individual_min = [];
for neuron_num = 1:size(total_means, 1)
    %Find std of the neuron in all trial types
    full_std = [];
    for p = 1:size(total_stds, 2)
        current_std = total_stds{neuron_num, p};
        full_std = horzcat(full_std, current_std);
    end
    individual_std = mean(full_std);
    
    full_means = [];
    for p = 1:size(total_means, 2)
        current_mean = total_means{neuron_num, p};
        full_means = horzcat(full_means, current_mean');
    end
    individual_mean = mean(full_means);
    
    for trial_cat = 1:size(total_stds, 2)
        raw_data = NeurondFF{neuron_num,trial_cat};
        if isempty(raw_data)
            continue
        end
        z_score = (raw_data - individual_mean)./individual_std;
        NeuronZScores{neuron_num,trial_cat} = z_score;
        %figure
        %plot(mean(z_score));
        NeuronZScores_mean{neuron_num,trial_cat} = mean(z_score,1);
    end
    %close all
    full_maxes = [];
    for p = 1:size(total_means, 2)
        raw_data = NeuronZScores{neuron_num,p};
        if isempty(raw_data)
            continue
        end
        current_max = max(mean(raw_data, 1));
        full_maxes = horzcat(full_maxes, current_max);
    end
    individual_max(neuron_num) = max(full_maxes);
    
    full_mins = [];
    for p = 1:size(total_means, 2)
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

if screen_neurons == 1
    neurons_to_delete = [];
end

for neuron_id = 1:size(NeuronZScores,1)   
    
    figure('Position', [100 100 1500 800],'Name','Summary','numbertitle','off');
    for trial_i=1:size(NeuronZScores,2)
        %sm = NeuronZScores{neuron_id,trial_i};
        sm = NeuronZScores_aligned{neuron_id,trial_i};
        if isempty(sm)
            continue;
        end
        
        Y = mean(sm,1);                         %4    7 - 1  , current liquid ind + 6
        Y_err = std(sm)/sqrt(size(sm,1));       %row, column, index 
        h(trial_i - pallet_cleansers) = subplot(4,(size(NeuronZScores,2) - pallet_cleansers),trial_i - pallet_cleansers + (size(NeuronZScores,2) - pallet_cleansers));
        %AreaPlot(time_x-time_x(1),Y,Y_err,'k',0.4,1); changed because I
        %dont have that function
        plot(time_x-time_x(1),Y);
        ylim([individual_min(neuron_id), individual_max(neuron_id)]);
        hold on
        for xk=1:length(X_line)
            plot([X_line(xk),X_line(xk)],[min(Y),max(Y)],'b--');
        end
        xlabel('Time (s)','FontSize', 8)
        ylabel('z-score(dF)','FontSize', 8)
        set(gca,'TickDir', 'out','xlim',[0,time_x(end)-time_x(1)],'xtick',0:2:time_x(end)-time_x(1),'FontSize', 8,'box','off');
        title(['Neuron',num2str(neuron_id),'-',title_text{trial_i - pallet_cleansers}],'FontSize', 6)
        
        ax1 = subplot(4,(size(NeuronZScores,2) - pallet_cleansers),trial_i - pallet_cleansers);
        
        % subtract baseline for betting visualization
        base_subtract = [];
        for line = 1:size(sm,1)
            base_subtract(line, :) = sm(line, :) - mean(sm(line, 3:(baseline_time*10 - 2)));
        end
        
        imagesc(time_x-time_x(1),1:size(sm,1),base_subtract);
        colormap(ax1,ColorMap);
        %colorbar;
        % xlabel('Time (s)','FontSize', 8)
        ylabel('Trial #','FontSize', 8)
        set(gca,'TickDir', 'out','xlim',[0,time_x(end)-time_x(1)],'xtick',0:2:time_x(end)-time_x(1),'YDir','reverse','FontSize', 8,'box','off');        
        caxis([-6 6])      
    end
    linkaxes(h,'xy');
    
    subplot(4,3,[7,8,9])
    neuron_trace = neuron.C(neuron_id, :);
    plot((1:length(neuron_trace))/10, neuron.C(neuron_id, :)*max(neuron.A(:, neuron_id)),'k');
    % plot((1:length(neuron_trace))/10,neuron_trace);
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
    
    if screen_neurons == 1
        keep = input('Keep this neuron? Enter for yes, n for no, q for quit','s');

        if strcmp(keep, 'n')
            neurons_to_delete = horzcat(neurons_to_delete, neuron_id);
            print([pwd,'\',psth_folder_name,'\NeuronPSTH_',num2str(neuron_id)],'-dpng');
        elseif strcmp(keep, '')
            print([pwd,'\',psth_folder_name,'\NeuronPSTH_',num2str(neuron_id)],'-dpng');
        elseif strcmp(keep, 'q')
            return
        end
    end
    
    close;
    
end

%% Here you want to preserve info. Don't delete the neurons! Keep for later
% ColorMap = exp_colormap('red-green',64);

NeuronZScores_copy = NeuronZScores_aligned;
NeuronZScores_mean_copy = mean_NeuronZScores_aligned;

NeuronZScores_copy(neurons_to_delete, :) = [];
NeuronZScores_mean_copy(neurons_to_delete, :) = [];

A = {};
for k=1:size(NeuronZScores_copy,2)
    
        A0 = NeuronZScores_mean_copy(:,k);
        A{k} = cell2mat(A0);
    
end
%Sort by 3 -2, fat

if pallet_cleansers
    first_sort_ind = 2;
else
    first_sort_ind = 1;
end
a_sort0 = A{first_sort_ind}(:,time_x>=stim_on & time_x<=(stim_off));
a_sort = mean(a_sort0,2);
[~,idx]=sortrows(a_sort);

figure('Position', [100 100 1200 400],'Name','Summary','numbertitle','off');
for k=(pallet_cleansers + 1):size(NeuronZScores,2)
    
    A{k} = A{k}(idx,:);
    
    subplot(2,size(NeuronZScores,2),k - pallet_cleansers);
    imagesc(time_x-time_x(1),1:size(A{k},1),A{k});
    colormap(ColorMap);
    colorbar;
    xlabel('Time (s)','FontSize', 10)
    ylabel('Neuron #','FontSize', 10)
    caxis([-3 3])
    set(gca,'TickDir', 'out',...
        'YDir','reverse','FontSize', 10,'box','off');
    title([title_text{k - pallet_cleansers}]);
    
    subplot(2,size(NeuronZScores,2),k - pallet_cleansers +size(NeuronZScores,2))
    plot(time_x-time_x(1),mean(A{k},1),'k','LineWidth',2)
    xlabel('Time (s)','FontSize', 10)
    ylabel('Norm activity','FontSize', 10)
    ylim([-1 3])
    xlim([0,time_x(end)-time_x(1)]);
    set(gca,'TickDir', 'out','FontSize', 16,'box','off');
    
end 

if 0
    NeuronZScores(neurons_to_delete, :) = [];
    NeurondFF(neurons_to_delete, :) = [];
    NeuronZScores_mean(neurons_to_delete, :) = [];
end
save(horzcat('z_scores_session_', num2str(session_ind), '.mat'), 'neurons_to_delete', ...
    'NeuronZScores', 'NeurondFF', 'NeuronZScores_mean', 'NeuronZScores_aligned', ...
    'mean_NeuronZScores_aligned', 'NeuronZScores_baseline_normed');

if isfield(SessionData.TrialSettings, 'first_liquid')
    print(gcf,horzcat('Heatmap_', first_liquid, '_sorted'),'-dpng');
    saveas(gcf,horzcat('Heatmap_', first_liquid, '_sorted.fig'));
else
    print(gcf,'Heatmap_sorted','-dpng');
    saveas(gcf,'Heatmap_sorted.fig');
end

