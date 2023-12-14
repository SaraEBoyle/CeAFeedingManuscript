function [hf_eating_lengths, s_eating_lengths] = get_eating_length()
    %% reads in behavior .csv file generated in BORIS and finds feeding length
    
    %% Load in all the data
    P = pwd;
    
    % load behavior and remove garbage from the beginning
    behavior = readtable('behavior.csv');
    times = table2cell(behavior(:,1));
    time_ind = strfind(times, 'Time');
    start_ind = find(~cellfun('isempty', time_ind));
    behavior(1:start_ind(2), :) = [];
    
    % Load the timestamps for the video data
    
    S = dir(fullfile(P,'*.csv'));
    N = {S.name};
    X = ~cellfun('isempty',strfind(N,'timestamp'));
    timestamps = N(X);
    timestamps = timestamps{1};
    
    vid_times = readtable(timestamps);
    vid_times = table2cell(vid_times(:, 17));
    
    translated_times = []; %translate bonsai's timestamps to seconds
    vid_starts = [];       %Keep track of video start times relative to beginning
    for ind = 1:size(vid_times, 1)
        bad_time = vid_times{ind};
        bad_time(1:11) = []; 
        hour = str2double(bad_time(1:2));
        min = str2double(bad_time(4:5));
        sec = str2double(bad_time(7:(end - 6)));
        total_sec = hour*60*60 + min*60 + sec;
        if ~isempty(translated_times)
            last_secs = translated_times(end);
            if total_sec - last_secs > 5
                vid_starts = [vid_starts;total_sec];
            end
        else
            vid_starts = total_sec;
        end
        translated_times = [translated_times;total_sec];
    end
    first_time = translated_times(1);
    translated_times = translated_times - first_time;
    vid_starts = vid_starts - first_time;
    
    %% find all the event time stamps
    events = table2array(behavior(:, 6));
    event_names = unique(events);
    high_fat_times = [];
    low_fat_times = [];
    missed_times = [];
    start_times = [];
    end_times = [];
    double_high_fat = [];
    double_different = [];
    double_sucrose = [];
    blue_dropped = [];
    pink_dropped = [];
    rejected = [];
    motor_turns = [];
    finishes = [];
    all_foods = [];
    all_food_identities = [];
    
    labels = {};
    count_fin = 0;
    for event_ind = 1:size(events, 1)
        event = events{event_ind, 1};
        switch event
            case 'High Fat'
                high_fat_times = [high_fat_times;str2double(cell2mat(behavior{event_ind, 1}))];
                labels{event_ind} = 'High Fat';
                all_foods = [all_foods; str2double(cell2mat(behavior{event_ind, 1}))];
                all_food_identities = [all_food_identities; 1];
                count_fin = 1;
            case 'Low Fat'
                low_fat_times = [low_fat_times;str2double(cell2mat(behavior{event_ind, 1}))];
                labels{event_ind} = 'Low Fat';
                all_foods = [all_foods; str2double(cell2mat(behavior{event_ind, 1}))];
                all_food_identities = [all_food_identities; 2];
                count_fin = 1;
            case 'Missed'
                missed_times = [missed_times;str2double(cell2mat(behavior{event_ind, 1}))];
                labels{event_ind} = 'Missed';
            case 'start trial'
                start_times = [start_times;str2double(cell2mat(behavior{event_ind, 1}))];
            case 'end trial'
                end_times = [end_times;str2double(cell2mat(behavior{event_ind, 1}))];
            case '2 pink pellets'
                double_high_fat = [double_high_fat;str2double(cell2mat(behavior{event_ind, 1}))];
                labels{event_ind} = 'Double High Fat';
            case '2 different pellets'
                double_different = [double_different;str2double(cell2mat(behavior{event_ind, 1}))];
                labels{event_ind} = 'Double HF and S';
            case '2 white pellets'
                double_sucrose = [double_sucrose;str2double(cell2mat(behavior{event_ind, 1}))];
                labels{event_ind} = 'Double S';
            case 'blue dropped'
                blue_dropped = [blue_dropped;str2double(cell2mat(behavior{event_ind, 1}))];
                labels{event_ind} = 'Blue dropped';
            case 'pink dropped'
                pink_dropped = [pink_dropped;str2double(cell2mat(behavior{event_ind, 1}))];
                labels{event_ind} = 'Pink dropped';
            case 'reject pellet'
                rejected = [rejected;str2double(cell2mat(behavior{event_ind, 1}))];
                labels{event_ind} = 'rejected';
            case 'motor turn'
                motor_turns = [motor_turns;str2double(cell2mat(behavior{event_ind, 1}))];
                labels{event_ind} = 'motor turn';
            case 'pellet finished'
                if count_fin
                    finishes = [finishes;str2double(cell2mat(behavior{event_ind, 1}))];
                    labels{event_ind} = 'pellet finished';
                    count_fin = 0;
                end
        end
    end
    
    %out = unique(cellfun(@num2str,labels,'uni',0));
    %halp = not(cellfun(@isempty, out));
    %unique_labels = out(halp);
    %% find all the trial lengths
    logitec_trial_lengths = end_times - start_times;
    
    %% Find the bpod timing data
    % timer 4 is LED, that is first
    %then timer 2 which is wire 1, the logitech cams through bonsai
    %then timer 1 which is bnc for 25 secs
    %the first run I ended the trial with timer 2, which controls logitec
    %% TODO
    %1. find the difference in timers 2 and 1, you'll have to add that to
    %the video timing data
    
    %2. Find the video frame of each event
    
    %3. load the timestamp data. Lookup actual times based on frames. You
    %will have to do this relative to the beginning time since this was
    %done using 2 different computers and so times won't match exactly
    
    %4. Collect food times = vid_food time - vid_start_time - timer dif
    %translated_times = translated_times - first_time;
    %vid_starts = vid_starts - first_time;
    
    
    %% Convert BORIS times to frame # by multiplying time by 30
    
    %translated_times is the bonsai timestamp recording
    %# of entries should match video frames
    high_fat_frames = round(high_fat_times * 30);
    translated_high_fat_times = translated_times(high_fat_frames);
    
    low_fat_frames = round(low_fat_times * 30);
    translated_low_fat_times = translated_times(low_fat_frames);
    
    missed_frames = round(missed_times * 30);
    translated_missed_times = translated_times(missed_frames);
    
    start_frames = round(start_times * 30);
    translated_start_times = translated_times(start_frames);
    
    end_frames = round(end_times * 30);
    translated_end_times = translated_times(end_frames);
    
    % get double times
    double_high_fat_frames = round(double_high_fat * 30);
    translated_double_high_fat = translated_times(double_high_fat_frames);
    
    double_different_frames = round(double_different * 30);
    translated_double_different = translated_times(double_different_frames);
    
    double_sucrose_frames = round(double_sucrose * 30);
    translated_double_sucrose = translated_times(double_sucrose_frames);
    
    blue_dropped_frames = round(blue_dropped * 30);
    translated_blue_dropped = translated_times(blue_dropped_frames);
    
    pink_dropped_frames = round(pink_dropped * 30);
    translated_pink_dropped = translated_times(pink_dropped_frames);
    
    rejected_frames = round(rejected * 30);
    translated_rejected = translated_times(rejected_frames);
    
    finished_frames = round(finishes * 30);
    translated_finishes = translated_times(finished_frames);
    
    all_foods_frames = round(all_foods * 30);
    translated_all_foods = translated_times(all_foods_frames);
    
    translated_eating_lengths = translated_finishes - translated_all_foods(1:(size(translated_finishes, 1)));
    
    f_trials = find(all_food_identities(1:size(translated_finishes, 1)) == 1);
    s_trials = find(all_food_identities(1:size(translated_finishes, 1)) == 2);
    hf_eating_lengths = translated_eating_lengths(f_trials);
    s_eating_lengths = translated_eating_lengths(s_trials);
    
    % sometimes I mess up marking finished times. delete trials with
    % mistakes
    hf_del_inds = find(hf_eating_lengths > 50);
    s_del_inds = find(s_eating_lengths > 50);
    hf_eating_lengths(hf_del_inds) = [];
    s_eating_lengths(s_del_inds) = [];
    
    
    %these I keep separate
    motor_frames = round(motor_turns * 30);
    translated_motor = translated_times(motor_frames);
    
    if ~isempty(translated_motor)
        relative_motor_time_boris = translated_motor - translated_start_times;
    else
        relative_motor_time_boris = 4*ones(size(translated_start_times, 1), size(translated_start_times, 2));
    end
    % start and end times as translated/scored through logitech
    bonsai_video_lengths = translated_end_times - translated_start_times;
    
    %Get trial numbers for each type of delivery
    all_food = [translated_high_fat_times; translated_low_fat_times; translated_missed_times; translated_double_high_fat; translated_double_different; translated_double_sucrose;translated_blue_dropped;translated_pink_dropped;translated_rejected];
    [sorted_food, sorted_inds] = sort(all_food);
    all_frames = [high_fat_frames; low_fat_frames; missed_frames; double_high_fat_frames; double_different_frames; double_sucrose_frames; blue_dropped_frames; pink_dropped_frames; rejected_frames];
    [sorted_frames, sorted_f_inds] = sort(all_frames);
    if sum(sorted_f_inds ~= sorted_inds) > 0
        disp('Uh oh. There is an issue translating the times of pellet delivery');
        times_sorted_by_frames = translated_times(sorted_frames);
        sorted_food = times_sorted_by_frames;
        sorted_inds = sorted_f_inds;
        disp('Saved pellet times according to frames, not translated times')
    end
    
    
    
    %Make a list of trial types, 1 = high fat, 2 = low fat, 3 = missed
    high_fat_block = ones(size(translated_high_fat_times, 1), 1);
    low_fat_block = 2*ones(size(translated_low_fat_times, 1), 1);
    missed_block = 3*ones(size(translated_missed_times, 1), 1);
    double_hf_block = 4*ones(size(translated_double_high_fat, 1), 1);
    double_hf_s_block = 5*ones(size(translated_double_different, 1), 1);
    double_s_block = 6*ones(size(translated_double_sucrose, 1), 1);
    blue_dropped_block = 7*ones(size(translated_blue_dropped, 1), 1);
    pink_dropped_block = 8*ones(size(translated_pink_dropped, 1), 1);
    translated_rejected_block = 9*ones(size(translated_rejected, 1), 1);
    %translated_motor_block = 10*ones(size(translated_motor, 1), 1);
    
    labels = {'High Fat', 'Sucrose', 'Missed', 'Double High Fat', 'Double HF and S', 'Double S', 'Blue Dropped', 'Pink Dropped', 'Rejected'};
    
    trial_ind_block = vertcat(high_fat_block, low_fat_block, missed_block, ...
        double_hf_block, double_hf_s_block, double_s_block, blue_dropped_block, ...
        pink_dropped_block, translated_rejected_block);
    trial_inds = trial_ind_block(sorted_inds);
    
    %Find time of food reception relative to first frame saved through
    %bonsai
    relative_sorted_food_times = sorted_food - translated_start_times;
    %around 8.
    %Find bpod distance between timer 2 start and end. The difference
    %between the bonsai video lengths and that is probably the offset to
    %use, or maybe not. I could usually see the light turn on
    avg_hf_length = mean(hf_eating_lengths);
    avg_s_length = mean(s_eating_lengths);
    save('eating lengths.mat', 'hf_eating_lengths', 's_eating_lengths', 'avg_hf_length', 'avg_s_length');
    
    
    
end