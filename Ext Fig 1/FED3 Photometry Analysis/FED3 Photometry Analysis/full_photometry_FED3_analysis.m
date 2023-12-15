function [] = full_photometry_FED3_analysis()
    %% KEY UPDATES- can now select multiple folders to analyze at once. You
    %% don't have to select the background image anymore, select folders instead.
    %% As always, make sure the key words are set, set the events you want to 
    %% plot, and make sure you set your fibers and fiber names.

    clearvars -except BpodSystem
    
    %% SET THESE VARIABLES
    bpod_key_word = 'FED3';
    analog_key_word = 'analog';
    signal_key_word = 'signal';
    %pellets = {'Low Fat', 'High Fat'}; % These need to match BORIS exactly
    %codes = {'Low Fat', 'High Fat'}; %
    pellets = {'HFHS', 'chow'}; % These need to match BORIS exactly
    codes = {'HFHS', 'chow'}; %
    file = horzcat(pwd, '\');
    tem = dir(horzcat(file, '*boris.csv')); % change "ior" to include name of BORIS file
    
    %Time before stimuli to check
    pre_times = {20, 20, 20};

    %Time after stimuli to check
    post_times = {40, 40, 40};

    %which timepoints to consider baseline
    baseline_times = {10, 10, 10}; 

    %Names of fibers you're checking. Order matters here
    fiber_names = {'CeL Left', 'CeL Right'};
    
    %List fibers you want to check in the analysis
    signal_fibers = [1, 2];
    
    %Your control cha1nnel- if using isosbestic channel it will be the same
    %as your signal fiber
    control_fiber = [1, 2];
    
    %% DON'T CHANGE BELOW HERE
    loaded_events = {}; 
    
    for event_ind = 1%:size(tem, 1)
        event_data = tem(event_ind);
        name = event_data.name;
        %[~, event_time] = xlsread(name);
        %halp = textscan(fid, '%f', 'Delimiter','\,')
        [event_time_pre] = readtable(name);
        event_time = table2cell(event_time_pre(17:end, 6));
        name = name(1:(end - 5));
        loaded_events{1, event_ind} = name;
        loaded_events{2, event_ind} = event_time;
    end    
    
    background_key_word = 'back';
    og_file = pwd;
    folders_to_analyze = uigetdir2(pwd);
    if isempty(folders_to_analyze)
        folders_to_analyze = {og_file};
    end
    for fold = 1:size(folders_to_analyze, 2)
        
        cd(folders_to_analyze{fold});
        file = pwd;
        %If 0, just do behavior data, no photometry
        start_file = dir(horzcat(file, '\', '*', signal_key_word, '*.csv'));
        if isempty(start_file)
            analyze_photometry_data = 0; 
        else
            analyze_photometry_data = 1; 
        end
        %exponential
        key_words = {bpod_key_word, analog_key_word, signal_key_word, background_key_word};
        AnimalTrackingSaraV7('exponential', loaded_events, pellets, codes, pre_times, post_times, baseline_times, control_fiber, signal_fibers, fiber_names, analyze_photometry_data, key_words)
        %AnimalTrackingSaraV5({16}, {5}, {5}, {5}, 2, 1, fiber_names, analyze_photometry_data, key_words)
    end
    cd(og_file);
    calculate_z_score_w_consistent_baselinesV2();
    remove_high_artifact_trials_V3()
    calculate_z_score_w_consistent_baselinesV2();
end