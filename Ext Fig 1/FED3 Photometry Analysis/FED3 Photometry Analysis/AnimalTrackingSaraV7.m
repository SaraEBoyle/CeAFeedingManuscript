%% animal trajectory tracking from .avi video
% ---> Inputs: video 
% ---> Outputs: mouse_parameter (column1: Xcenter; column1: Ycenter; column3: Distance; column1: Velocity)
% ---> Sara Boyle sboyle@cshl.edu
%% clear workspace
%% TODO Count up number of seconds spent nose poking
%% TODO Sound off tracking is sus
%% TODO- potentially another measure- average point in sound delivery where mouse spends 90% of time on platform
%Idea is that this can be a measure of how well the mouse has learned. If
%less than 18 the mouse has learned
%% TODO- separate platform crosses that are first- ie first platform entry after sound.
%% Find first platform exit after a shock- purple has a long lasting increase in activity to shock. Will she stay on the platform until the signal returns to baseline?
%% First platform exit after sound ends
%% Number of nose pokes
%% Amount of water delivered
%% TODO platform entries not during sound
%% TODO get first water reward

%%TODO classify time spent in the corners
function [] = AnimalTrackingSaraV7(bleaching_correction, event_groups, pellets, codes, pre_times, post_times, baseline_times, control_fiber, signal_fibers, fiber_names, photometry_data, key_words)
%This is a list of all possible events to align to- do not change it
event_labels = pellets;
bpod_key_word = key_words{1};
analog_key_word = key_words{2};
signal_key_word = key_words{3};
background_key_word = key_words{4};
event_names = {'Food Delivered', 'Food Lifted'};

%% photometry variables
if strcmp(bleaching_correction, 'exponential')
    %This corrects by fitting an exponential
    correct_for_bleaching = 1;

    %This corrects by using a sliding window
    moving_avg = 0;
else
    %This corrects by fitting an exponential
    correct_for_bleaching = 0;

    %This corrects by using a sliding window
    moving_avg = 1;
end


frame_rate = 30;
fps_photometry = 10;

pre_ITI_length = 150;
seconds_to_extend = 0; 
sound_on = 0;
sound_only = 0;
trigger1 = 1;
plot_single_trials = 0;
plot_heatmap_and_average = 1;
save_graphs = 1;
volume_to_check = 1000;
stim_duration = .2;
sound_offset = 0;
time_after_sound_onset_to_check = 0;
record_sample = 1;
plot_location = 1;
behavior_data = 1;
flip_arena = 1; % just leave this on 1

file = horzcat(pwd, '\');
tem = dir(horzcat(file, '*', signal_key_word, '*.csv')); 
if isempty(tem)
    photometry_data = 0;
end

%% select EthoVision data from excel
% clear nam
halp = dir(file);
for thing = 1:size(halp, 1)
    if and(and(contains(halp(thing).name, background_key_word), ~contains(halp(thing).name, '.fig')), ~contains(halp(thing).name, '.mat'))
        nam = [file, halp(thing).name];  % full name of the data file
        [~, file_nm, file_type] = fileparts(nam);
        break
    end
end
start_file = (fullfile(file, horzcat('*', background_key_word)));

%if ~exist('nam', 'var') || isempty(nam)
%    [file_nm, dir_nm] = uigetfile(fullfile('*.png'));
%    if dir_nm==0
%        fprintf('no file was selected. Try again.\N');
%        
%        return;
%    end
%    nam = [dir_nm, file_nm];  % full name of the data file
%    [~, file_nm, file_type] = fileparts(nam);
%end

[photometry_start, start_times] = get_start_times(pwd, analog_key_word);

%% Load in data
if photometry_data == 1
    [start_times, photo_data, SessionData] = load_data(bpod_key_word, signal_key_word, analog_key_word, file, photometry_data);
end

%% Read in video file
file = pwd;
file = horzcat(file, '\');
%video = readframe(vidObj, 1); %reads only the frames specified by index .

close all;
 
pellet_delivered_times = [];
pellet_delivered_ind = [];
food_lifted_times = [];
food_lifted_ind = [];
pick_up_times = [];
pick_up_ind = [];

tem = dir(horzcat(file, '*', bpod_key_word, '*.mat')); 
if isempty(tem)
    tem = dir(horzcat(file, '*', bpod_key_word, '*.mat'));
    if isempty(tem)
        tem = dir(horzcat(file, '*behavior*.mat'));
    end
end
if isempty(tem)
    error('set signal_key_word to a phrase that appears in the bpod file')
end

name = tem(1).name;
folder = tem(1).folder;
full_path = horzcat(folder, '\', name);
disp(name);
load(full_path); 

%Find time of events
for q = 1:min(SessionData.nTrials, size(start_times, 2))
    trial_start = start_times(q);
    %Record nose poke times
       
    %Record food delivery times
    if isfield(SessionData.RawEvents.Trial{1, q}.States, 'FoodDelivered')
        if ~isnan(eval(horzcat('SessionData.RawEvents.Trial{1, q}.States.', 'FoodDelivered', '(1)')))
            reward_time = eval(horzcat('SessionData.RawEvents.Trial{1, q}.States.', 'FoodDelivered'));
            reward_time = reward_time(:,1);
            if ~isnan(reward_time)        
                pellet_delivered_times = vertcat(pellet_delivered_times, reward_time + trial_start);
                pellet_delivered_ind = horzcat(pellet_delivered_ind, q); 
            end
        end
    end

    %Record pellet lifted times
    if isfield(SessionData.RawEvents.Trial{1, q}.States, 'FoodLifted')
        if ~isnan(eval(horzcat('SessionData.RawEvents.Trial{1, q}.States.', 'FoodLifted', '(1)')))
            lifted_time = eval(horzcat('SessionData.RawEvents.Trial{1, q}.States.', 'FoodLifted'));
            food_lifted_time = lifted_time(:,1);

            if ~isnan(lifted_time)        
                food_lifted_times = vertcat(food_lifted_times, food_lifted_time + trial_start - 5);
                food_lifted_ind = horzcat(food_lifted_ind, q);               
            end
        end
    end
end

%photometry data isn't pulled in until here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
event_times = {};
pellet_order = event_groups{2, 1};
for x = 1:size(pellets, 2)
    current_food_ind = [];
    for d = 1:size(pellet_order, 1)
        if strcmp(pellet_order{d}, codes{1, x})
            current_food_ind = horzcat(current_food_ind, d);
        end
    end
    event_times{x} = current_food_ind;
end

save(horzcat(file, 'event_times.mat'), 'event_times', 'pellets');

raw_signal_dFFs = {};
raw_control_dFFs = {};
event_group_ind = 1;


for fiber = 1:size(signal_fibers, 2)
    %for each fiber

    if photometry_data == 1
        %%So it's the event time minus the first time value in data to get the
        %%timestamp in the video

        for signal = 1:2
            %For signal and control
            [~, photo_data, ~] = load_data(bpod_key_word, signal_key_word, analog_key_word, file, photometry_data);
            if signal == 1
                fiber_name = horzcat(fiber_names{signal_fibers(1, fiber)}, ' signal');
                %to get signal data
                %Trim your data if needed
                [photo_data, auto_flo] = trim_data(photo_data, signal_fibers(1, fiber), 1);

                %De-interleave and subtract autofluorescence
                [gcamp, iso] = de_interleave(1, 0, file, fiber_name, photo_data, signal_fibers(1, fiber), auto_flo, 1, signal);
%% TODO MAKE CORRECT BLEACHING SAVE EACH WHOLE SESSION IN RIGHT PLACE 
                %correct bleaching
                [gcamp, iso] = correct_bleaching(horzcat(file, fiber_name, '\'), correct_for_bleaching, gcamp, iso, moving_avg, frame_rate, 1, 0, pellet_delivered_times);

            elseif signal == 2
                %to get control channel
                fiber_name = horzcat(fiber_names{signal_fibers(1, fiber)}, ' control');
                current_control_fiber = control_fiber(1, fiber);
                if signal_fibers(1, fiber) == current_control_fiber%(fiber)
                    %if it's the same channel as signal
                    gcamp = iso; %just swap the channels
                else
                    %otherwise you have to re analyze
                    %fiber = current_control_fiber;
                    %Trim your data if needed

                    %signal_file = dir(horzcat(dir_nm, '*', signal_key_word, '*.csv'));
                    %photo_data = dlmread(horzcat(dir_nm, signal_file.name));
                    [IAmLazy, photo_data, alsotrash, trash] = load_data(bpod_key_word, signal_key_word, analog_key_word, file, photometry_data);

                    [photo_data, auto_flo] = trim_data(photo_data, current_control_fiber, 1);

                    %De-interleave and subtract autofluorescence
                    [gcamp, iso] = de_interleave(1, 0, file, fiber_name, photo_data, current_control_fiber, auto_flo, 1, signal);

                    %correct bleaching
                    [gcamp, iso] = correct_bleaching(horzcat(file, fiber_name, '\'), correct_for_bleaching, gcamp, iso, moving_avg, frame_rate, 1, 0);

                end

            end

            for ev_ind = 1:size(event_times, 2)
                %for each event
                events_to_analyze = event_times{1, ev_ind};
                pre_time = pre_times{ev_ind};
                post_time = post_times{ev_ind};
                baseline_time = baseline_times{ev_ind};

                seconds_per_trial = pre_time + post_time;
                %fiber_name = fiber_names{signal_fibers(1, fiber)};

                if and(event_group_ind == 1, signal == 1)
                    raw_signal_dFFs{fiber} = gcamp;
                elseif and(event_group_ind == 1, signal == 2)
                    raw_control_dFFs{fiber} = gcamp;
                end

                for event_type_ind = [1, 2]
                    if isempty(event_times{ev_ind})
                        continue
                    else
                        alignment = event_labels{ev_ind};     
                    end
                    if ~exist(horzcat(file, fiber_name), 'dir')
                        mkdir(horzcat(file, fiber_name));
                        mkdir(horzcat(file, fiber_name, '\', alignment));
                        mkdir(horzcat(file, fiber_name, '\', alignment, '\', event_names{1, event_type_ind}));
                        
                    elseif ~exist(horzcat(file, fiber_name, '\', alignment), 'dir')
                        mkdir(horzcat(file, fiber_name, '\', alignment));
                        mkdir(horzcat(file, fiber_name, '\', alignment, '\', event_names{1, event_type_ind}));
                    elseif ~exist(horzcat(file, fiber_name, '\', alignment, '\', event_names{1, event_type_ind}), 'dir')
                        mkdir(horzcat(file, fiber_name, '\', alignment, '\', event_names{1, event_type_ind}));
                    end
                    filesub = horzcat(file, fiber_name, '\', alignment, '\', event_names{1, event_type_ind});
                    cd(filesub);


                    trials = {};
                    iso_trials = {};
                    OG_trials = {};
                    pre_ITIs = {};
                    post_ITIs = {};
                    US_frame = {}; 
                    no_US_delivered = [];
                    US_delivered = [];
                    ind = 1;
                    US_ind = 1;
                    US_frames = [];
                    CS_frames = [];
                    checked_trial_types = [];
                    other_side = [];
                    %Just for US receiving
                    event_trials = {};
                    event_i = 1;
                    event_pre_ITIs = {};
                    event_post_ITIs = {};
                    event_frame = {};          %Record frame air puff happens                              
                    event_delivered = [];

                    %Just for no US trials
                    no_event_trials = {};
                    no_event_i = 1;
                    no_event_pre_ITIs = {};
                    no_event_post_ITIs = {};

                    if event_type_ind == 1
                        %First do delivery times
                        if events_to_analyze(end) > size(pellet_delivered_times, 1)
                            events_to_analyze(end) = [];
                        end
                        times_delivered = pellet_delivered_times(events_to_analyze);
                        event_times_temp = pellet_delivered_times(events_to_analyze);
                        if events_to_analyze(end) <= length(food_lifted_ind)
                            event_trial_ind = food_lifted_ind(events_to_analyze);
                        else
                            events_to_analyze(end) = [];
                            event_trial_ind = food_lifted_ind(events_to_analyze);
                        end
                    elseif event_type_ind == 2
                        %Then do food lifted times
                        times_delivered = food_lifted_times(events_to_analyze);
                        event_times_temp = food_lifted_times(events_to_analyze);
                        event_trial_ind = food_lifted_ind(events_to_analyze);
                    end

                    %start_times = alt_start_times; 
                    start_times = times_delivered - pre_time;
                    end_times = times_delivered + post_time;
                    evs = 1;
                    %TODO INVESTIGATE
                    gcamp_clone = gcamp;
                    if ~isempty(end_times)
                        if end_times(1) > gcamp(end, 1)
                            start_times = [];
                            event_trial_ind = [];
                            event_times_temp = [];
                        end
                    end
                    for y = 1:size(start_times, 1)
                        %Segment the trials based on start and end times
                        %pre_ITI_length = 5;
                        post_ITI_length = 150;   
                        %if evs > length(event_times)

                        starty = start_times(y);
                        endy = end_times(y);
                        indices1 = find(gcamp(:,1) >= starty);
                        if isempty(indices1)
                            disp('This event happened after the end of the recording...')
                            break
                        end
                        if 1
                            indices2 = find(iso(:,1) >= starty);
                            trial2= iso(indices2,:);
                        end
                        indices3 = find(gcamp_clone(:,1) >= starty);
                        trial1=gcamp(indices1,:);
                        trial3 = gcamp_clone(indices3, :);

                        indices1 = find(trial1(:,1) <= endy);
                        if 1
                            indices2 = find(trial2(:,1) <= endy);
                            trial2=trial2(indices2,:);
                        end
                        indices3 = find(trial3(:,1) <= endy);
                        trial1=trial1(indices1,:);
                        if size(trial1, 1) < fps_photometry*seconds_per_trial - 2
                            event_times_temp(1) = [];
                            disp('Event too soon. Had to skip first event')
                            continue
                        end
                        trial3=trial3(indices3,:);
                        pre_ITI=gcamp(find(gcamp(:,1) >= (start_times(y) - pre_ITI_length)),:); %Bigger than start - 15
                        pre_ITI=pre_ITI(find(pre_ITI(:,1) <= (start_times(y))), :); %smaller than start of trial?
                        post_ITI= gcamp(find(gcamp(:,1) >= end_times(y)), :); %Greater than end
                        post_ITI_clone= post_ITI(find(post_ITI(:,1) <= (end_times(y) + post_ITI_length)),:);  %Less than end + ITI
                        post_ITI= post_ITI(find(post_ITI(:,1) <= (end_times(y) + post_ITI_length)),:);  %Less than end + ITI

                        US_ind = event_i;
                        if ~isempty(event_times_temp)            %If not too many USs
                            if y ~= 1
                                if  and(event_times_temp(1) < starty, event_times_temp(1) > end_times(y - 1) + seconds_to_extend)
                                    disp(horzcat('The mouse licked too late into trial ', num2str(y - 1), '. Skipping that event'));
                                    event_times_temp = event_times_temp(2:end);
                                    to_take_out = find(event_trial_ind == y - 1);
                                    event_trial_ind(to_take_out:end) = event_trial_ind(to_take_out:end) - 1;
                                    event_trial_ind(to_take_out) = [];
                                    ind = ind - 1;
                                end
                            else
                                if size(trial1, 1) < ((pre_time + post_time) * 9.5)
                                    event_times_temp(1) = [];
                                    disp('recording is off. Had to skip first event');
                                    continue
                                end
                            end
                            US_frames = find(trial1(:, 1) >= event_times_temp(1)); %Find US frames
                             if isempty(US_frames)
                                US_frames = find(post_ITI_clone(:, 1) >= event_times_temp(1)); %Find US frames
                                US_frames = US_frames + size(trial1, 1);
                                if ~isempty(US_frames)
                                    if US_frames(1) > size(trial1,1) + seconds_to_extend*10 -20
                                        disp(horzcat('The mouse licked too late into trial ', num2str(event_trial_ind(y)), '. Skipping that event'));
                                        event_times_temp = event_times_temp(2:end);
                                        to_take_out = find(event_trial_ind == event_trial_ind(y));
                                        event_trial_ind(to_take_out:end) = event_trial_ind(to_take_out:end) - 1;
                                        event_trial_ind(to_take_out) = [];
                                        ind = ind - 1;
                                        US_frames = find(trial1(:, 1) >= event_times_temp(1)); %Find US frames
                                        if isempty(US_frames)
                                            disp('This mouse licked too late into previous trial then too late into current')
                                        end
                                    end   
                                end
                            end
                        end


                        if or(~isempty(US_frames), sound_only)                               %If there's a US
                            if sound_only
                                US_frame{ind} = sound_frame{ind};
                            else
                                if US_frames(1, 1) ~= 1
                                    US_frame{ind} = US_frames(1, 1);                 %Find the US onset frame
                                else
                                    US_frame{ind} = 0;
                                end
                            end                
                            US_delivered = horzcat(US_delivered, y);         %record trial # with US
                            trials{ind} = trial1(:, [1 2]);                  %Record trial data in total list
                            if trigger1
                                iso_trials{ind} = trial2(:, [1 2]);
                            end
                            OG_trials{ind} = trial3(:, [1 2]);
                            pre_ITIs{ind} = pre_ITI;                         %Record pre ITI data
                            post_ITIs{ind} = post_ITI;                       %Record post ITI data

                            if ismember(y, US_delivered)
                                event_trials{event_i} = trial1(:, [1 2]);        %Record trial as air trial
                                event_pre_ITIs{event_i} = pre_ITI;               %Record air pre ITI
                                event_post_ITIs{event_i} = post_ITI;             %Record are post ITI
                                event_frame{event_i} = US_frame{ind};          %Record frame air puff happens                              
                                event_delivered = horzcat(event_delivered, event_i);
                                event_i = event_i + 1;
                                event_times_temp = event_times_temp(2:end);

                            end
                            US_frames = [];                                  %Clear US frames
                            ind = ind + 1;
                        else                                                 %If there's no US
                            if ind == 0
                                ind = 1;
                            end
                            US_frame{ind} = 0;
                            other_side = horzcat(other_side, event_trial_ind(y));
                            no_US_delivered = horzcat(no_US_delivered, ind);

                            if ismember(event_trial_ind(y), no_event_trial_ind)
                                no_event_trials{no_event_i} = trial1(:, [1 2]);        %Record trial as air trial
                                no_event_pre_ITIs{event_i} = pre_ITI;               %Record air pre ITI
                                no_event_post_ITIs{event_i} = post_ITI;             %Record are post ITI
                                no_event_i = no_event_i + 1;
                            end
                            trials{ind} = trial1(:, [1 2]);
                            if trigger1
                                iso_trials{ind} = trial2(:, [1 2]);
                            end
                            OG_trials{ind} = trial3(:, [1 2]);
                    %        checked_trial_types = horzcat(checked_trial_types, trial_types(event_trial_ind(y)));
                            pre_ITIs{ind} = pre_ITI;
                            post_ITIs{ind} = post_ITI;
                            ind = ind + 1;
                        end
                    end

                    %% Plot single trials 
                    if plot_single_trials
                        if ~isempty(event_trial_ind)
                            plot_single_trials_V2(1, event_trials, alignment, event_labels, stimulus_times_temp, stim_duration);
                        end
                    end
                    if plot_single_trials
                        if ~isempty(no_event_trial_ind)
                            plot_single_trials_V2(0, no_event_trials, alignment, event_labels, no_stimulus_times_temp, stim_duration);
                        end
                    end

                    %% Sort and align trials.
                    close all
                    if plot_heatmap_and_average
                        if ~exist('Alignment Data\', 'dir')
                            mkdir(horzcat(file, 'Alignment Data\'));
                        end
                        align_to_this = alignment;
                        if and(~isempty(event_trial_ind), ~isempty(event_trials))
                            good_pre_ITIs = pre_ITIs;
                            good_post_ITIs = post_ITIs;

                            [US_trials, US_target_frame, true_seconds_per_trial, sec_offsets, true_sec_offsets] = sort_and_align_trials_V7(event_trials, good_pre_ITIs, good_post_ITIs, ...
                                seconds_per_trial, US_frame, 10, seconds_to_extend);

                            frame_offsets = sec_offsets * 10;
                            trial_start_times = start_times;
                            trial_order = 1:size(start_times, 1);
                            save(horzcat(file, 'Alignment Data\', alignment, ' alignment data.mat'), 'sec_offsets', 'frame_offsets', 'event_trial_ind', 'true_sec_offsets', 'event_times', 'trial_start_times', 'trial_order');
                        end
                        no_US_trials = [];
                        if ~isempty(event_trial_ind)
                            US_trials_percent = US_trials * 100;
                            baselines = [];
                            for r = 1:size(US_trials_percent, 1)
                                cur = US_trials_percent(r, :);
                                baseline = cur(1:((baseline_time - 1)*fps_photometry - 5));
                                baselines(r, :) = baseline;
                            end
                            plot_average('010', fiber_name, US_trials_percent, true_seconds_per_trial, US_target_frame, 0, 0, alignment, .2);
                            ylabel('dF/F (%)', 'FontSize', 20);
                            if strcmp(alignment, 'left')
                                title(horzcat(fiber_name, ': % dF/F in ', left_spout, ' trials'));
                            elseif strcmp(alignment, 'right')
                                title(horzcat(fiber_name, ': ', right_spout, ' trials'));
                            else
                                title(horzcat(fiber_name, ': % dF/F in US delivered ', alignment, ' trials'));
                            end

                            if save_graphs
                                if strcmp(alignment , 'left')
                                    savefig('left_average_dFF.fig')
                                elseif strcmp(alignment , 'right')
                                    savefig('right_average_dFF.fig')
                                else
                                    savefig(horzcat(alignment, '_average.fig'));
                                end
                            end
                        else
                            disp(horzcat('Found no ', alignment, ' trials. Exited.'))
                            cd(file);

                        end

                        %% Plot trial by trial responses


                        baseline_no_US_trials = [];
                        las_baseline_no_US_trials = [];
                        no_las_baseline_no_US_trials = [];
                        z_no_US_trials = [];
                        %end
                        if ~isempty(event_trial_ind)            
                            plot_heatmap(US_trials_percent, fiber_name, true_seconds_per_trial, US_target_frame, 0, alignment, pre_time);
                            title(horzcat(fiber_name, ': ', alignment, ' trials'));

                            if strcmp(alignment, 'left')
                                savefig('left_delivery_heatmap.fig');
                            elseif strcmp(alignment, 'right')
                                savefig('right_delivery_heatmap.fig');
                            else
                                savefig(horzcat(alignment, '_heatmap.fig'));
                            end              
                        end

                        increments = size(US_trials_percent, 2);
                        if increments == 0
                            increments = size(no_US_trials, 2);
                        end
                        x_axis = ((1:increments)/increments)*true_seconds_per_trial;
                        if volume_to_check == 1000
                            laser_time = x_axis(US_target_frame);
                            x_axis = x_axis - laser_time;
                            if ~isempty(event_trial_ind)
                                if strcmp(alignment, 'left')
                                    left_trials = US_trials_percent;
                                    save('dFF_data.mat', 'left_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines', 'pre_time', 'post_time');
                                elseif strcmp(alignment, 'right')
                                    right_trials = US_trials_percent;
                                    save('dFF_data.mat', 'right_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines', 'pre_time', 'post_time');
                                else
                                    save('dFF_data.mat', 'US_trials_percent', 'US_target_frame', 'x_axis', 'baselines', 'auto_flo', 'pre_time', 'post_time');
                                end
                            end

                            %US
                            if isempty(event_trial_ind)
                                baseline_US_trials = [];
                            else
                                for x = 1:size(US_trials_percent, 2)
                                    if isnan(US_trials_percent(1,x))
                                        US_trials_percent(1,x) = 0;
                                    end
                                end
                                baseline_US_trials = US_trials_percent(:,20:(US_target_frame - 20));
                            end
                            no_event_trial_ind = [];
                            if and(~isempty(no_event_trial_ind), ~isempty(no_US_trials))
                                %baseline_trials = vertcat(baseline_US_trials, baseline_no_US_trials);
                                baseline_trials = baselines;
                                std_avg = mean(std(baseline_trials));
                                avg = mean(mean(baseline_trials));
                                diffs = (no_US_trials - avg);
                                z_no_US_trials = (diffs ./ std_avg);
                                diffs = (US_trials_percent - avg);
                                z_US_trials = (diffs ./ std_avg); 
                                plot_average(010, fiber_name, z_US_trials, true_seconds_per_trial, US_target_frame, CS_target_frame, sound_only, alignment, stim_duration);
                                title(horzcat(fiber_name, ': ', alignment, ' trials'));
                                ylabel('z-score(dF/F)', 'FontSize', 20);
                                if strcmp(alignment, 'left')
                                    savefig('left_z_score_average.fig');
                                elseif strcmp(alignment, 'right')
                                    savefig('right_z_score_average.fig');
                                else
                                    savefig('US_z_score_average.fig');
                                end
                                if and(~isempty(no_event_trial_ind), ~isempty(no_US_trials))
                                    plot_average(010, fiber_name, z_no_US_trials, true_seconds_per_trial, CS_target_frame, CS_target_frame, sound_only, alignment, stim_duration);
                                    title(horzcat(fiber_name, ': ', alignment, ' trials'));
                                    ylabel('z-score(dF/F)', 'FontSize', 20);
                                    if strcmp(alignment, 'left')
                                        z_right_trials = z_no_US_trials;
                                        z_left_trials = z_US_trials;
                                        save('z_scores.mat', 'z_left_trials', 'x_axis', 'US_target_frame');                                 
                                    elseif strcmp(alignment, 'right')
                                        z_right_trials = z_US_trials;
                                        z_left_trials = z_no_US_trials;
                                        save('z_scores.mat', 'z_right_trials', 'x_axis', 'US_target_frame');
                                    else
                                        save('z_scores.mat', 'z_US_trials', 'z_no_US_trials', 'x_axis', 'US_target_frame');
                                    end
                                end
                            elseif(~isempty(event_trial_ind))
                                %%baseline_trials = vertcat(baseline_US_trials);





                                baseline_trials = baselines;
                                std_avg = mean(std(baseline_trials));
                                avg = mean(mean(baseline_trials));
                                diffs = (US_trials_percent - avg);
                                z_US_trials = (diffs ./ std_avg); 
                                %post_laser = round(mean(end_CS_frame)) + 1;
                                %laser_effect = mean(z_US_trials(:,post_laser:post_laser + 20));
                    %                    save('z_scores.mat', 'z_US_trials', 'z_no_US_trials', 'x_axis');

                                plot_average('010', fiber_name, z_US_trials, true_seconds_per_trial, US_target_frame, 0, 0, alignment, stim_duration);
                                title(horzcat(fiber_name, ': ', alignment, ' trials'));
                                ylabel('z-score(dF/F)', 'FontSize', 20);
                                if strcmp(alignment, 'left')
                                    savefig('left_z_score_average.fig');
                                    save('dFF_data.mat', 'left_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines', 'pre_time', 'post_time');
                                    z_left_trials = z_US_trials;
                                    save('z_scores.mat', 'z_left_trials', 'x_axis', 'US_target_frame');
                                elseif strcmp(alignment, 'right')
                                    savefig('right_z_score_average.fig');
                                    save('dFF_data.mat', 'right_trials', 'US_target_frame', 'x_axis', 'baselines', 'pre_time', 'post_time');
                                    z_right_trials = z_US_trials;
                                    save('z_scores.mat', 'z_right_trials', 'x_axis', 'US_target_frame');
                                else
                                    %save('dFF_data.mat', 'right_trials', 'US_target_frame', 'x_axis', 'baselines');
                                    save(horzcat(filesub, '\z_scores.mat'), 'z_US_trials', 'x_axis', 'US_target_frame');
                                    z_right_trials = z_US_trials;
                                    savefig(horzcat(filesub, '\US_z_score_average.fig'));
                                end
                                if and(~isempty(no_event_trial_ind), ~isempty(no_US_trials))
                                    plot_average(010, fiber_name, z_no_US_trials, true_seconds_per_trial, CS_target_frame, CS_target_frame, sound_only, alignment, stim_duration);
                                    title(horzcat(fiber_name, ': ', alignment, ' trials'));
                                    ylabel('z-score(dF/F)', 'FontSize', 20);
                                    if strcmp(alignment, 'left') 
                                        save('dFF_data.mat', 'left_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines', 'pre_time', 'post_time');
                                        z_right_trials = z_no_US_trials;
                                        z_left_trials = z_US_trials;
                                        save('z_scores.mat', 'z_left_trials', 'x_axis', 'US_target_frame');
                                    elseif strcmp(alignment, 'right')
                                        save('dFF_data.mat', 'right_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines', 'pre_time', 'post_time');
                                        z_right_trials = z_US_trials;
                                        z_left_trials = z_no_US_trials;
                                        save('z_scores.mat', 'z_right_trials', 'x_axis', 'US_target_frame');
                                    else
                                        save('z_scores.mat', 'z_US_trials', 'z_no_US_trials', 'x_axis', 'US_target_frame');
                                    end
                                end
                            end    
                        else
                            if align_to_US
                                if ~isnan(US_target_frame)
                                    laser_time = x_axis(US_target_frame);
                                else
                                    laser_time = x_axis(72);
                                    US_target_frame = 72;
                                end
                            else
                                laser_time = x_axis(CS_target_frame);
                            end
                            if strcmp(alignment, 'right')
                                right_trials = US_trials_percent;
                            elseif strcmp(alignment, 'left')
                                left_trials = US_trials_percent;
                            end
                            %x_axis = x_axis - laser_time;
                            laser_time = x_axis(US_target_frame);
                            x_axis = x_axis - laser_time;
                            if strcmp(alignment, 'right')
                                save('dFF_data.mat', 'right_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines', 'pre_time', 'post_time');
                            elseif strcmp(alignment, 'left')
                                save('dFF_data.mat', 'left_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines', 'pre_time', 'post_time');
                            else
                                save('dFF_data.mat', 'US_trials_percent', 'no_US_trials', 'US_target_frame', 'CS_target_frame', 'x_axis', 'baselines', 'auto_flo', 'pre_time', 'post_time');
                            end

                            %US
                            if or(strcmp(alignment, 'right'), strcmp(alignment, 'left'))
                                %baseline_no_US_trials = no_US_trials(:,1:(US_target_frame - 20));
                                baseline_US_trials = US_trials_percent(:,1:(US_target_frame - 20));
                            end
                            %baseline_trials = vertcat(baseline_US_trials, baseline_no_US_trials);
                            baseline_trials = baselines;
                            std_avg = mean(std(baseline_trials));
                            avg = mean(mean(baseline_trials));
                            %diffs = (no_US_trials - avg);
                            %z_no_US_trials = (diffs ./ std_avg);
                            diffs = (US_trials_percent - avg);
                            z_US_trials = (diffs ./ std_avg); 
                            if ~isempty(event_trial_ind)
                                plot_average(010, fiber_name, z_US_trials, true_seconds_per_trial, US_target_frame, CS_target_frame, sound_only, alignment, stim_duration);
                                if strcmp(alignment, 'left')
                                    correct_title = 'left';
                                elseif strcmp(alignment, 'right')
                                    correct_title = 'right';
                                else
                                    correct_title = alignment;
                                end

                                title(horzcat(fiber_name, ': ', correct_title, ' trials'));
                                ylabel('z-score(dF/F)', 'FontSize', 20);
                                if save_graphs
                                    if strcmp(alignment, 'left')
                                            savefig('left_z_score_average.fig');
                                            z_left_trials = z_US_trials;
                                    elseif strcmp(alignment, 'right')
                                        savefig('right_z_score_average.fig');
                                        z_right_trials = z_US_trials;
                                    else
                                        savefig('US_z_score_average.fig');
                                    end
                                end
                            end

                            if strcmp(alignment, 'left')
                                save('z_scores.mat', 'z_left_trials', 'x_axis', 'US_target_frame');
                            elseif strcmp(alignment, 'right')
                                save('z_scores.mat', 'z_right_trials', 'x_axis', 'US_target_frame');
                            else
                                save('z_scores.mat', 'z_US_trials', 'z_no_US_trials', 'x_axis', 'US_target_frame');
                            end
                        end
                    end
                            % It takes 0.0001 s for arduino to read analog input. .00025 s for bpod
                            % state changes. .0001 to do output action
                            cd(file);
                end
            end
        %random_message(2);
        end
    end
    event_group_ind = event_group_ind + 1;
end

save(horzcat(file, 'raw_data.mat'), 'raw_signal_dFFs', 'raw_control_dFFs');

end
function [photometry_start_time, start_times] = get_start_times(file, analog_key_word)
    file = horzcat(file, '\'); 
    err = 1;
    start_file = dir(horzcat(file, '*', analog_key_word, '*.csv'));
    if ~isempty(start_file)
        photometry_start = dlmread(horzcat(file, start_file.name));
        %photometry_start(1:1177104) = []; %% DELETE THIS LATER
    else
        photometry_start = [0];
    end
    photometry_start_time = photometry_start(err);
    if 1
        %Convert analog in file to seconds
        photometry_start_clone = photometry_start;
        photometry_start_clone = (photometry_start_clone(err:end) - photometry_start(err,1))/1000;
        if photometry_start_clone == 0
            error('Change the analog_key_word to something in the name of your file.');
        end
        if photometry_start_clone(2) - photometry_start_clone(1) > 1
            photometry_start_clone(1) = [];
            photometry_start(1,:) = [];
            photometry_start_clone = photometry_start_clone - photometry_start_clone(1);
        end
        photometry_start_clone = photometry_start_clone';

        prev = 1;
        cur = 1;
        alt_start_times = [photometry_start_clone(1)];
    else
        %Convert centroid and analog file to seconds
        photometry_start_clone = photometry_start;
        photometry_start_clone = (photometry_start_clone - photometry_start(1,1))/1000;
        if photometry_start_clone(2) - photometry_start_clone(1) > 1
            photometry_start_clone(1) = [];
            photometry_start(1,:) = [];
            photometry_start_clone = photometry_start_clone - photometry_start_clone(1);
        end
        photometry_start_clone = photometry_start_clone';

        prev = 1;
        cur = 1;
        alt_start_times = [photometry_start_clone(1)];
    end

    for s = 2:size(photometry_start_clone,2)
        dif = photometry_start_clone(s) - photometry_start_clone(prev);
        if dif > 5
            alt_start_times = horzcat(alt_start_times, photometry_start_clone(s));
            disp(s);
        end
        prev = prev + 1;
    end

    %% Use this to keep track of matlab events
    %set 0 point for imaging beginning. This coincides with trial 1.
    startRef = photometry_start(1);

    %convert time to relative time (s)
    start_times = alt_start_times;
end

function [platform_matrix, sound_data_time_adjusted] = find_sound_platform_position(tracking_data, background_image)
%input- tracking data N x 3 (x,y,time from onset of sound) output- matrix N
%x 2 (boolean, timestamp)

    sound_data_time_adjusted = tracking_data;
    prev = 0;
    index = 1;
    for y = 1:size(tracking_data, 1)
        cur = tracking_data(y, 3);
        if prev == 0
            prev = cur;
            continue
        end
        if cur > prev + 20
            sound_data_time_adjusted(index:(y-1), 3) = sound_data_time_adjusted(index:(y-1), 3) - sound_data_time_adjusted(index, 3);
            index = y;
            prev = cur;
        end
    end

    platform_matrix = zeros(size(tracking_data, 1), 2);
    for i = 1:size(tracking_data, 1)
        vert_num = round(tracking_data(i, 2));
        horz_num = round(tracking_data(i, 1));
        time = tracking_data(i, 3);
        
        
        %I2 = im2double(background_image);
        %I2((abs(vert_num - 3)):(vert_num + 3), (abs(horz_num - 3)):(horz_num + 3), 1) = 0;
        %imshow(I2);
        
        % Find if the mouse is on the platform
        if vert_num > size(background_image, 1)/2 + 10
            %if in bottom half
            if horz_num < (size(background_image, 2)/2 - 10)
                %if in left half
                %disp('In platform zone');
                platform_matrix(i, 1:2) = [1, time];
            else
                %disp('Not in platform zone');
                platform_matrix(i, 1:2) = [0, time];
            end
        else
            %disp('Not in platform zone');
            platform_matrix(i, 1:2) = [0, time];
        end
    end
end
