function [] = align_to_event_frames()
%% Written by Sara Boyle (sboyle@cshl.edu)
%% Aligns to seconds of event because I didn't do this bonsai thing right
    %% Aligns photometry data to poke events, given the video frames
    
    %Setup: set keyword(s) of frame files and of video time file
    bleaching_correction = 1; %1 is exponential, 0 is sliding average, I recommend 1
    pre_times = 5; %Time before event to plot
    post_times = 10; % Time after event to plot
    baseline_times = 3; %This is how much to use as baseline, can't be bigger than pre_times
    control_fiber = [1]; % Use these settings for 2 fibers
    signal_fibers = [1];
    fiber_names = {'CeL Left', 'CeL Right'}; % Change these to label the folders
    % Read csv from Boris
    %T = readtable('boris.csv'); %change this to match the Boris file name
    [num,txt,raw] = xlsread('boris.csv'); %change this to match the Boris file name
    X = raw;
    %% Don't edit below here
    file = pwd;
    file = horzcat(file, '\');
    
    pre_times = ones(1, 10) .* pre_times;
    post_times = ones(1, 10) .* post_times;
    baseline_times = ones(1, 10) .* baseline_times;
    
    %% load data from Boris
    
    
    %Convert weird table to cell
    %X = table2cell(T);
    
    % Search for the row that starts with "Time"
    for w = 1:size(X, 1)
        current = X{w, 1};
        if strcmp(current, 'Time')
            start_row = w;
            break
        end
    end
    %Remove junk from beginning
    X(1:(start_row - 1), :) = [];
    
    %Find unique Behavior events
    for w = 1:size(X, 2)
        current = X{1, w};
        if strcmp(current, 'Behavior')
            behavior_ind = w;
            break
        end
    end
    
    %Find whether the event is a start or stop or point
    for w = 1:size(X, 2)
        current = X{1, w};
        if strcmp(current, 'Status')
            status_ind = w;
            break
        end
    end
    
    behavior_column = X(:, behavior_ind);
    behavior_names = unique(behavior_column);
    
    %Remove "Behavior"
    behavior_names(1) = [];
    behavior_names = behavior_names';
    behavior_column(1) = [];
    
    %Find status
    status_column = X(:, status_ind);
    status_column(1) = [];
    
    %Combine behavior and status to get names
    event_names = {};
    for q = 1:size(status_column, 1)
        current_behavior = behavior_column{q};
        current_status = status_column{q};
        event_names{1, q} = horzcat(current_behavior, ' ', current_status);
    end
    
    behavior_names = unique(event_names);
    
    behavior_times = X(2:end,1);
    loaded_events = behavior_names;
    
    %for each event name, find times, create loaded_events cell
    for name_ind = 1:size(behavior_names, 2)
        name = behavior_names{name_ind};
        cur_times = [];
        for cur_name_ind = 1:size(event_names, 2)
            cur_name = event_names{cur_name_ind};
            if strcmp(cur_name, name)
                cur_time = behavior_times{cur_name_ind};
                cur_times = horzcat(cur_times, cur_time);
            end
        end
        loaded_events{2, name_ind} = cur_times;
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %% Use the keywords to extract data
    %{
    for y = 1:size(keywords, 2)
        keyword = keywords{y};
        tem = dir(horzcat(file, '*', keyword, '*.xlsx')); 
        if isempty(tem)
            disp(horzcat('Keyword ', keyword, ' cannot be found'));
        end
        
        if size(tem, 1) > 1
            for x = 1:size(tem, 1)
                current = tem(x, 1);
                
                if strcmp(current.name(1, [1,2]), 'no')
                    disp('delete it');
                    tem(x) = [];
                end
            end
        end
        
        %Save keyword and imported data
        imported_data{1, y} = keyword;
        imported_data{2, y} = xlsread(tem.name);
        
    end
    %}
    %% Now import time file and convert the time to seconds
    %tem = dir(horzcat(file, '*', timing_file_keyword, '*.csv')); 
    %data = dlmread(horzcat(file, tem.name));
    
    %timing_data = (data - data(1,1))\1000;
    timing_data = 1/30;
    halp = 0:timing_data:780;
    save('loaded_events.mat', 'loaded_events');
    AnimalTrackingSaraV7(bleaching_correction, halp, loaded_events, pre_times, post_times, baseline_times, control_fiber, signal_fibers, fiber_names, 1)
end