function [] = align_to_event_frames()
%% Written by Sara Boyle (sboyle@cshl.edu)
%% Aligns to seconds of event because I didn't do this bonsai thing right
    %% Aligns photometry data to poke events, given the video frames
    
    %Setup: set keyword(s) of frame files and of video time file
    keywords = {'white chocolate', 'wc drops'};
    bleaching_correction = 1; %1 is exponential, 0 is sliding average
    pre_times = [20, 20, 20, 20, 20, 20, 20, 20]; %Should have same number as excel files
    post_times = [40, 40, 40, 40, 40, 40, 40, 40]; % This is seconds before event to plot
    baseline_times = [10, 10, 10, 10, 10, 10, 10, 10]; %This is how much to use as baseline
    control_fiber = [1, 2]; % If using gCaMP use these settings
    signal_fibers = [1, 2];
    fiber_names = {'CeL Left', 'CeL Right'};
    
    %% Don't edit below here
    file = pwd;
    file = horzcat(file, '\');
    
    %% load data from Boris
    % Read csv from Boris
    %T = readtable('behavior.csv'); 
    T = xlsread('white chocolate.xlsx');
    %Convert weird table to cell
    %X = table2cell(T);
    % Search for the row that starts with "Time"
    
    %for each event name, find times, create loaded_events cell
    loaded_events{1, 1} = 'White chocolate'; %names, top row
    loaded_events{2, 1} = T;    
    
    Q = xlsread('wc drops.xlsx');
    loaded_events{1, 2} = 'WC Drops'; %names, top row
    loaded_events{2, 2} = Q;  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
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
    oldAnimalTrackingSaraV7(bleaching_correction, halp, loaded_events, pre_times, post_times, baseline_times, control_fiber, signal_fibers, fiber_names, 1)
end