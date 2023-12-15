% Uses cell_registered_struct to call a cell index map, then prints the 
%cells only in saline trials, the cells only in amphetamine trials, cells 
%from both. Prints cells that respond to CS or US in both and either. You
%can also plot average baseline response from all neurons over the 3
%seconds before the CS, and individual average baseline amplitude in saline
%condition and amph condition. Can plot cells that respond to US in both.
%Finally, you can choose to graph a slopegraph of individual neuron
%responses.

%Returns ordered_events_sal, ordered_events_amph, average_baselines_sal,
%average_baselines_amph

%% Load files
%Choose the file from cell reg
load('/Users/sboyle/Documents/Li Lab/Image analysis/SST6/cell_registration/cellRegistered.mat');

%Choose files with your dF/F data. Important to do this in order
files = {'/Users/sboyle/Documents/Li Lab/Image analysis/SST6/onlysaline/', ...
    '/Users/sboyle/Documents/Li Lab/Image analysis/SST6/amph/', ...
    '/Users/sboyle/Documents/Li Lab/Image analysis/SST5/onlysaline/', ...
    '/Users/sboyle/Documents/Li Lab/Image analysis/SST5/amph/'
    };

%Load response pools from Response_Pools_SB.mat
response_pools_1 = load('/Users/sboyle/Documents/Li Lab/Image analysis/SST6/onlysaline/response_pools.mat');
response_pools_2 = load('/Users/sboyle/Documents/Li Lab/Image analysis/SST6/amph/response_pools.mat');
%response_pools_3 = load('/Users/sboyle/Documents/Li Lab/Image analysis/SST5and6onlywatertrials/Cell_reg/Testing_response_pool/response_pools_3.mat');
%response_pools_4 = load('/Users/sboyle/Documents/Li Lab/Image analysis/SST5and6onlywatertrials/Cell_reg/Testing_response_pool/response_pools_4.mat');
pools = [response_pools_1, response_pools_2];

%% Define variables
%Variables you need to define each time
n_water_Trials = 20; %Assuming same # in saline and amphetamine trials
n_shock_Trials = 15;
window = 69:99;%69:99; %The area to test dF. 
framerate = 10;
trial_duration = 20; %in seconds

%Initializing variables
cell_to_index_map = cell_registered_struct.cell_to_index_map;
cells_only_in_sal = [];
cells_only_in_amph = [];
cells_in_both = []; 
cs_responders = [];
us_responders = [];
total_water_cs_responders = [];
water_cs_responders_sal = [];
water_cs_responders_amph = [];
total_water_us_responders = [];
total_shock_cs_responders = [];
total_shock_us_responders = [];
shock_us_responders_amph = [];
shock_cs_responders_amph = [];
shock_us_responders_sal = [];
shock_cs_responders_sal = [];
water_us_responders_sal = [];
water_us_responders_amph = [];
water_cs_responders_both = [];
water_us_responders_both = [];
 
% Variables for event counting
f1 = (files(1));
f1 = load(fullfile(cell2mat(f1), 'Sorted.mat'));
f2 = (files(2));
f2 = load(fullfile(cell2mat(f2), 'Sorted.mat'));
responsive = 0; % 0 for all neurons, 1 if checking responsive, 2 for unresponsive
DF_amph = f2.DF;
DF_sal = f1.DF;
num_window_frames = size(window, 2);
FrameLength_water = size(DF_sal,2)/n_water_Trials;
%FrameLength_shock = size(DF_sal,2)/n_water_Trials;
total_frames_water = n_water_Trials*FrameLength_water;
total_frames_amph = n_shock_Trials*framerate*trial_duration;
frames_per_trial = trial_duration*framerate;

%% Group neurons by response to US and CS and by treatment
% Make array of responsive neurons for each dataset. odd = sal, even = amph
count = 1;
for y=pools 
    if mod(count,2) == 1
        water_cs_responders{count} = y.water_cs_response_pool; 
        water_us_responders{count} = y.water_us_response_pool;
    end
    if mod(count,2) == 0
        water_cs_responders{count} = y.water_cs_response_pool; 
        water_us_responders{count} = y.water_us_response_pool;
        shock_cs_responders{count} = y.shock_cs_response_pool; 
        shock_us_responders{count} = y.shock_us_response_pool;
    end
    count = count + 1;
end
% Records which neurons respond to CS or US using their cellreg index
for x = 1:size(cell_to_index_map, 1)
    sal_index = cell_to_index_map(x, 1);
    amph_index = cell_to_index_map(x, 2);
    indeces = num2str([double(sal_index>0), double(amph_index>0)]*1);
    switch indeces
        case '1  0'
            % In sal, not amph
            water_cs_sal = sum(ismember(cell2mat(water_cs_responders(1)), sal_index))>0; %If the neuron responds to cs in sal
            water_us_sal = sum(ismember(cell2mat(water_us_responders(1)), sal_index))>0; %If the neuron responds to us in sal
            %shock_cs_sal = sum(ismember(cell2mat(shock_cs_responders_sal(1)), sal_index))>0; %If the neuron responds to shock cs in amph
            %shock_us_sal = sum(ismember(cell2mat(shock_us_responders(1)), sal_index))>0; %If the neuron responds to shock us in amph
            cells_only_in_sal = horzcat(cells_only_in_sal, x);
            if water_cs_sal
                total_water_cs_responders = horzcat(total_water_cs_responders, x);
                water_cs_responders_sal = horzcat(water_cs_responders_sal, x);
            end
            if water_us_sal
                total_water_us_responders = horzcat(total_water_us_responders, x);
                water_us_responders_sal = horzcat(water_us_responders_sal, x);            
            end
            %shock
            %if shock_cs_sal
            %    total_shock_cs_responders = horzcat(total_shock_cs_responders, x);
            %    water_cs_responders_sal = horzcat(water_cs_responders_sal, x);
            %end
            %if shock_us_sal
            %    total_shock_us_responders = horzcat(total_shock_us_responders, x);
            %    shock_us_responders_sal = horzcat(shock_us_responders_sal, x);            
            %end
        case '0  1'
            % In amph not sal
            water_cs_amph = sum(ismember(cell2mat(water_cs_responders(2)), amph_index))>0; %If responds to cs in amph
            water_us_amph = sum(ismember(cell2mat(water_us_responders(2)), amph_index))>0; %If responds to us in amph        
            shock_cs_amph = sum(ismember(cell2mat(shock_cs_responders(2)), amph_index))>0; %If responds to cs in amph
            shock_us_amph = sum(ismember(cell2mat(shock_us_responders(2)), amph_index))>0; %If responds to us in amph        
            cells_only_in_amph = horzcat(cells_only_in_amph, x);
            if water_cs_amph
                total_water_cs_responders = horzcat(total_water_cs_responders, x);
                water_cs_responders_amph = horzcat(water_cs_responders_amph, x);
            end
            if water_us_amph
                total_water_us_responders = horzcat(total_water_us_responders, x);
                water_us_responders_amph = horzcat(water_us_responders_amph, x);            
            end
            if shock_cs_amph
                total_shock_cs_responders = horzcat(total_shock_cs_responders, x);
                shock_cs_responders_amph = horzcat(shock_cs_responders_amph, x);
            end
            if shock_us_amph
                total_shock_us_responders = horzcat(total_shock_us_responders, x);
                shock_us_responders_amph = horzcat(shock_us_responders_amph, x);            
            end
        case '1  1'
            % In both
            water_cs_sal = sum(ismember(cell2mat(water_cs_responders(1)), sal_index))>0; %If the neuron responds to cs in sal
            water_cs_amph = sum(ismember(cell2mat(water_cs_responders(2)), amph_index))>0; %If responds to cs in amph
            water_us_sal = sum(ismember(cell2mat(water_us_responders(1)), sal_index))>0; %If the neuron responds to us in sal
            water_us_amph = sum(ismember(cell2mat(water_us_responders(2)), amph_index))>0; %If responds to us in amph        
            shock_cs_amph = sum(ismember(cell2mat(shock_cs_responders(2)), amph_index))>0; %If responds to cs in amph
            shock_us_amph = sum(ismember(cell2mat(shock_us_responders(2)), amph_index))>0;
            cells_in_both = horzcat(cells_in_both, x);
            if water_cs_amph
                total_water_cs_responders = horzcat(total_water_cs_responders, x);
                water_cs_responders_amph = horzcat(water_cs_responders_amph, x);
            end
            if water_us_amph
                total_water_us_responders = horzcat(total_water_us_responders, x);
                water_us_responders_amph = horzcat(water_us_responders_amph, x);            
            end
            if water_cs_sal
                total_water_cs_responders = horzcat(total_water_cs_responders, x);
                water_cs_responders_sal = horzcat(water_cs_responders_sal, x);
            end%
            if water_us_sal
                total_water_us_responders = horzcat(total_water_us_responders, x);%
                water_us_responders_sal = horzcat(water_us_responders_sal, x);% 
            end
            if and(water_cs_sal, water_cs_amph)
                water_cs_responders_both = horzcat(water_cs_responders_both, x);
                water_cs_responders_sal = horzcat(water_cs_responders_sal, x);
                water_cs_responders_amph = horzcat(water_cs_responders_amph, x);
            end
            if and(water_us_sal, water_us_amph)
                water_us_responders_both = horzcat(water_us_responders_both, x);
            end
            %%%%%% Shock
            if shock_cs_amph
                total_shock_cs_responders = horzcat(total_shock_cs_responders, x);
                shock_cs_responders_amph = horzcat(shock_cs_responders_amph, x);
            end
            if shock_us_amph
                total_shock_us_responders = horzcat(total_shock_us_responders, x);
                shock_us_responders_amph = horzcat(shock_us_responders_amph, x);
            end
        otherwise
            disp('Something is wrong');
    end
end

%% Report results
disp(horzcat(['Cells only in saline trials: ' newline 'Total: ', num2str(size(cells_only_in_sal, 2)) newline], num2str(cells_only_in_sal)));
disp(horzcat(['Cells only in amphetamine trials: ' newline 'Total: ', num2str(size(cells_only_in_amph, 2)) newline], num2str(cells_only_in_amph)));
disp(horzcat('cells that respond to water CS in sal: ', num2str(water_cs_responders_sal)));
disp(horzcat('cells that respond to water CS in amph: ', num2str(water_cs_responders_amph)));
disp(horzcat('cells that respond to water US in sal: ', num2str(water_us_responders_sal)));
disp(horzcat('cells that respond to water US in amph: ', num2str(water_us_responders_amph)));
disp(horzcat('cells that respond to shock CS in amph: ', num2str(shock_cs_responders_amph)));
disp(horzcat('cells that respond to shock US in amph: ', num2str(shock_us_responders_amph)));
disp(horzcat(['Cells in both saline and amphetamine trials: ' newline 'Total: ', num2str(size(cells_in_both, 2)) newline], num2str(cells_in_both)));
disp(horzcat(['Cells that respond to water CS in both: ' newline 'Total: ', num2str(size(water_cs_responders_both, 2)), ' out of ', num2str(size(cells_in_both, 2)), ' total cells in both' newline], num2str(water_cs_responders_both)));
disp(horzcat(['Cells that respond to water US in both: ' newline 'Total: ', num2str(size(water_us_responders_both, 2))]));

%Print neurons that respond to US in both conditions
for n=water_us_responders_both
    disp(horzcat('Saline and amphetamine indeces of US responders: ', '[', num2str(cell_to_index_map(n, [1,2])), ']'));    
end

%% Graph all neurons that respond to US or CS in both conditions
if ~isempty(water_us_responders_both)   
    %If cells respond to US in both sal and amph, plot them 
    for cells = water_us_responders_both
        sal_ind = cell_to_index_map(n, 1);
        amph_ind = cell_to_index_map(n, 2);
        %Get average response in sal condition
        neuron_s = f1.DF(sal_ind,:);
        b = 1;
        e = frames_per_trial;
        av = zeros(1,frames_per_trial);
        for i = (1:n_water_Trials)
            av = av + neuron_s(b:e);
            b = b + frames_per_trial;
            e = e + frames_per_trial;
        end
        %Get average response in amph condition
        neuron_a = f2.DF(amph_ind,:);
        b = 1;
        e = frames_per_trial;
        av2 = zeros(1,frames_per_trial);
        for i = (1:n_water_Trials)
            av2 = av2 + neuron_a(b:e);
            b = b + frames_per_trial;
            e = e + frames_per_trial;
        end
        %Plot average US responses in sal and amph
        f = figure;
        p = uipanel('Parent',f,'BorderType','none'); 
        p.Title = horzcat('Neuron ', num2str(sal_ind), ' US Response'); 
        p.TitlePosition = 'centertop'; 
        p.FontSize = 12;
        p.FontWeight = 'bold';
        subplot(1,1,1, 'Parent',p)       
        plot(vertcat(av, av2).')      
        legend({'Saline', 'Amphetamine'});
        title('dF/F');
        line([100 100], [0 1.8]);
        line([125 125], [0 1.8]);
    end
end
if ~isempty(water_cs_responders_both)       
    %If cells respond to CS in both sal and amph, plot them
    for cells = water_cs_responders_both
        sal_ind = cell_to_index_map(n, 1);
        amph_ind = cell_to_index_map(n, 2);
        %Get average response in sal condition
        neuron_s = f1.DF(sal_ind,:);
        b = 1;
        e = frames_per_trial;
        av = zeros(1,frames_per_trial);
        for i = (1:n_water_Trials)
            av = av + neuron_s(b:e);
            b = b + frames_per_trial;
            e = e + frames_per_trial;
        end
        %Get average response in amph condition
        neuron_a = f2.DF(amph_ind,:);
        b = 1;
        e = frames_per_trial;
        av2 = zeros(1,frames_per_trial);
        for i = (1:n_water_Trials)
            av2 = av2 + neuron_a(b:e);
            b = b + frames_per_trial;
            e = e + frames_per_trial;
        end
        %Plot both average responses
        f = figure;
        p = uipanel('Parent',f,'BorderType','none'); 
        p.Title = horzcat('Neuron ', num2str(sal_ind), ' CS Response'); 
        p.TitlePosition = 'centertop'; 
        p.FontSize = 12;
        p.FontWeight = 'bold';
        subplot(1,1,1, 'Parent',p)       
        plot(vertcat(av, av2).')      
        legend({'Saline', 'Amphetamine'});
        title('dF/F');
        line([100 100], [0 1.8]);
        line([125 125], [0 1.8]);
    end
end
%% Calculate baselines
average_baselines_sal = 0;
average_baselines_amph = 0;
average_baselines_sal_n = 0;
average_baselines_amph_n = 0;
average_neuron_baseline_sal = [];
average_neuron_baseline_amph_w = [];
average_neuron_baseline_amph_s = [];
average_baselines_amph_w = 0;
average_baselines_amph_s = 0;
average_neuron_baseline_only_amph_w = [];
average_neuron_baseline_only_amph_s = [];
average_neuron_baseline_only_sal_w = [];
average_baselines_only_amph_w = 0;
average_baselines_only_amph_s = 0;
average_baselines_only_sal = 0;
% Calculates average dF/F over 3 second baseline period in both conditions
%Record the average dF/F per neuron per trial
sal_per_trial = [];
amph_per_trial_w = [];
amph_per_trial_s = [];
co = 1;
%%{

for y=cells_only_in_sal
sal_index = cell_to_index_map(y, 1);
    sal_neuron = f1.DF(sal_index,:); %saline
    b = 1;
    e = frames_per_trial;
    av = zeros(1,frames_per_trial); %holds sum of trials, frame by frame
    for i = (1:n_water_Trials)
        %for every trial,
        trial = sal_neuron(b:e);
        %Save avg baseline response per neuron per trial
        only_sal_per_trial(co, i) = mean(trial(window)); 
        av = av + trial;
        b = b + frames_per_trial;
        e = e + frames_per_trial;
    end
    %Add's each neuron's average baseline in saline
    average_baselines_only_sal = average_baselines_only_sal + mean(av(window)/n_water_Trials); 
    %Makes a list of all neuron's average baseline in saline
    average_neuron_baseline_only_sal_w = horzcat(average_neuron_baseline_only_sal_w, mean(av(window)/n_water_Trials));
    co = co + 1;
end
co = 1;

for y=cells_only_in_amph
    %Load amph data
    amph_neuron_w = f2.DF(cell_to_index_map(y, 2),1:total_frames_water);
    if n_shock_Trials ~= 0
        amph_neuron_s = f2.DF(cell_to_index_map(y, 2),total_frames_water + 1:total_frames_water + total_frames_amph);
    end
    b = 1;
    e = frames_per_trial;
    av2 = zeros(1,frames_per_trial);
    for i = (1:n_water_Trials)
        %for every trial,
        trial = amph_neuron_w(b:e);
        %Save avg baseline response per neuron per trial, amph, water
        only_amph_per_trial_w(co, i) = mean(trial(window));
        av2 = av2 + trial;
        b = b + frames_per_trial;
        e = e + frames_per_trial;
    end
    
    b = 1;
    e = frames_per_trial;
    av2_s = zeros(1,frames_per_trial);
    for i = (1:n_shock_Trials)
        trial = amph_neuron_s(b:e);
        only_amph_per_trial_s(co, i) = mean(trial(window));
        av2_s = av2_s + trial;
        b = b + frames_per_trial;
        e = e + frames_per_trial;
    end
    %Add's each neuron's average baseline in amph water US
    average_baselines_only_amph_w = average_baselines_only_amph_w + mean(av2(window)/n_water_Trials);
    %Add's each neuron's average baseline in amph shock US
    if n_shock_Trials
        average_baselines_only_amph_s = average_baselines_only_amph_s + mean(av2_s(window)/n_shock_Trials);
    end
    %Makes a list of all neuron's average baseline in amph
    average_neuron_baseline_only_amph_w = horzcat(average_neuron_baseline_only_amph_w, mean(av2(window)/n_water_Trials));
    if n_shock_Trials ~= 0
        average_neuron_baseline_only_amph_s = horzcat(average_neuron_baseline_only_amph_s, mean(av2_s(window)/n_shock_Trials));
    end
    co = co + 1;
end
co = 1;
for y=cells_in_both %adds up dF/F of neurons in both cases
    sal_index = cell_to_index_map(y, 1);
    sal_neuron = f1.DF(sal_index,:); %saline
    b = 1;
    e = frames_per_trial;
    av = zeros(1,frames_per_trial); %holds sum of trials, frame by frame
    for i = (1:n_water_Trials)
        %for every trial,
        trial = sal_neuron(b:e);
        %Save avg baseline response per neuron per trial
        sal_per_trial(co, i) = mean(trial(window)); 
        av = av + trial;
        b = b + frames_per_trial;
        e = e + frames_per_trial;
    end
    %Add's each neuron's average baseline in saline
    average_baselines_sal = average_baselines_sal + mean(av(window)/n_water_Trials); 
    %Makes a list of all neuron's average baseline in saline
    average_neuron_baseline_sal = horzcat(average_neuron_baseline_sal, mean(av(window)/n_water_Trials));
    %Load amph data
    amph_neuron_w = f2.DF(cell_to_index_map(y, 2),1:total_frames_water);
    if n_shock_Trials ~= 0
        amph_neuron_s = f2.DF(cell_to_index_map(y, 2),total_frames_water + 1:total_frames_water + total_frames_amph);
    end
    b = 1;
    e = frames_per_trial;
    av2 = zeros(1,frames_per_trial);
    for i = (1:n_water_Trials)
        %for every trial,
        trial = amph_neuron_w(b:e);
        %Save avg baseline response per neuron per trial, amph, water
        amph_per_trial_w(co, i) = mean(trial(window));
        av2 = av2 + trial;
        b = b + frames_per_trial;
        e = e + frames_per_trial;
    end
    
    b = 1;
    e = frames_per_trial;
    av2_s = zeros(1,frames_per_trial);
    for i = (1:n_shock_Trials)
        trial = amph_neuron_s(b:e);
        amph_per_trial_s(co, i) = mean(trial(window));
        av2_s = av2_s + trial;
        b = b + frames_per_trial;
        e = e + frames_per_trial;
    end
    %Add's each neuron's average baseline in amph water US
    average_baselines_amph_w = average_baselines_amph_w + mean(av2(window)/n_water_Trials);
    %Add's each neuron's average baseline in amph shock US
    if n_shock_Trials
        average_baselines_amph_s = average_baselines_amph_s + mean(av2_s(window)/n_shock_Trials);
    end
    %Makes a list of all neuron's average baseline in amph
    average_neuron_baseline_amph_w = horzcat(average_neuron_baseline_amph_w, mean(av2(window)/n_water_Trials));
    if n_shock_Trials ~= 0
        average_neuron_baseline_amph_s = horzcat(average_neuron_baseline_amph_s, mean(av2_s(window)/n_shock_Trials));
    end
    co = co + 1;
end 
%%}

%To find average baseline dF/F of all neurons, not just those in both
%{
for y=horzcat(cells_in_both,cells_only_in_sal)
    sal_neuron = f1.DF(y,:); %saline
    b = 1;
    e = frames_per_trial;
    av = zeros(1,frames_per_trial); %holds sum of trials, frame by frame
    for i = (1:n_water_trials)
        %for every trial,
        trial = sal_neuron(b:e);
        %Save avg baseline response per neuron per trial
        sal_per_trial(co, i) = mean(trial(69:99)); 
        av = av + trial;
        b = b + frames_per_trial;
        e = e + frames_per_trial;
    end
    %Add's each neuron's average baseline in saline
    average_baselines_sal = average_baselines_sal + mean(av(69:99)/n_water_Trials); 
    %Makes a list of all neuron's average baseline in saline
    average_neuron_baseline_sal = horzcat(average_neuron_baseline_sal, mean(av(69:99)/n_water_Trials));
    %Load amph data
    co = co + 1;
end
co = 1;
for y=horzcat(cells_in_both,cells_only_in_amph)
    amph_neuron_w = f2.DF(cell_to_index_map(y, 2),1:4000);
    amph_neuron_s = f2.DF(cell_to_index_map(y, 2),4001:7000);
    b = 1;
    e = frames_per_trial;
    av2 = zeros(1,frames_per_trial);
    for i = (1:n_water_Trials)
        %for every trial,
        trial = amph_neuron_w(b:e);
        %Save avg baseline response per neuron per trial, amph, water
        amph_per_trial_w(co, i) = mean(trial(69:99));
        av2 = av2 + trial;
        b = b + frames_per_trial;
        e = e + frames_per_trial;
    end
    b = 1;
    e = frames_per_trial;
    av2_s = zeros(1,frames_per_trial);
    for i = (1:n_shock_Trials)
        trial = amph_neuron_s(b:e);
        amph_per_trial_s(co, i) = mean(trial(69:99));
        av2_s = av2_s + trial;
        b = b + frames_per_trial;
        e = e + frames_per_trial;
    end
    %Add's each neuron's average baseline in amph water US
    average_baselines_amph_w = average_baselines_amph_w + mean(av2(69:99)/n_water_Trials);
    %Add's each neuron's average baseline in amph shock US
    average_baselines_amph_s = average_baselines_amph_s + mean(av2_s(69:99)/n_shock_Trials);
    %Makes a list of all neuron's average baseline in amph
    average_neuron_baseline_amph_w = horzcat(average_neuron_baseline_amph_w, mean(av2(69:99)/n_water_Trials));
    average_neuron_baseline_amph_s = horzcat(average_neuron_baseline_amph_s, mean(av2_s(69:99)/n_shock_Trials));
    co = co + 1;
end
%}

%divide by number of neurons
average_baselines_sal = average_baselines_sal/size(average_neuron_baseline_sal, 2);
average_baselines_amph_w = average_baselines_amph_w/size(average_neuron_baseline_amph_w, 2);
if n_shock_Trials ~= 0
    average_baselines_amph_s = average_baselines_amph_s/size(average_neuron_baseline_amph_s, 2);
end
%Plot first few shock trials to see change
disp(horzcat('Average all neuron baseline dF/F in saline water trials: ', num2str(mean(average_baselines_sal))));
disp(horzcat('Average baseline dF/F in amph water trials: ', num2str(mean(average_baselines_amph_w))));
disp(horzcat('Average baseline dF/F in amph shock trials: ', num2str(mean(average_baselines_amph_s))));

%% Counts events in each neuron. Compares between saline and amph conditions
% Counts number of "events" per neuron during baseline firing, determined 
% by dF/F exceeding sd of each neuron's dF/F. To run, change responsive to 
% either 0 (all neurons), 1 (responsive), or 2 (unresponsive).

%% If you're doing responsive or unresponsive neurons
responsive_amph_neurons = response_pools_2.water_us_response_pool;
unresponsive_amph_neurons = 1:size(DF_amph',2);
unresponsive_amph_neurons(:, responsive_amph_neurons) = []; %take out responsive neurons
responsive_sal_neurons = response_pools_1.water_us_response_pool;
unresponsive_sal_neurons = 1:size(DF_sal',2);
unresponsive_sal_neurons(:, responsive_sal_neurons) = []; %take out responsive neurons
%% Setting up variables
if responsive == 0 %if you want to run this for all neurons
    sal_neurons = 1:size(DF_sal,1);
    neuron_num_sal = size(DF_sal,1);
    neuron_num_amph = size(DF_amph,1);
    amph_neurons = 1:size(DF_amph,1);
    neuron_type = ' all';
elseif responsive == 1 %if you want responsive neurons
    sal_neurons = responsive_sal_neurons;
    neuron_num_sal = length(responsive_sal_neurons);
    amph_neurons = responsive_amph_neurons;
    neuron_num_amph = length(responsive_amph_neurons);
    neuron_type = ' responsive';
elseif responsive == 2 %unresponsive
    sal_neurons = unresponsive_sal_neurons;
    neuron_num_sal = length(unresponsive_sal_neurons);
    amph_neurons = unresponsive_amph_neurons;
    neuron_num_amph = length(unresponsive_amph_neurons);
    neuron_type = ' unresponsive';
end

%% Find std of dF/F for each neuron. Set event threshold to 2*std
co = 1;
sd_sal = zeros(size(cells_in_both, 2), 2);
sd_amph_w = zeros(size(cells_in_both, 2), 2);
sd_amph_s = zeros(size(cells_in_both, 2), 2);
for y=cells_in_both
    sal_index = cell_to_index_map(y, 1);
    sal_neuron = f1.DF(sal_index,:); %saline %pulls DF/F for each frame for 1 neuron
    %sal_neuron = sal_neuron';  % trial X frames
    sd_baselines_sal = std(sal_neuron);
    sd_sal(co, 1) = y; sd_sal(co, 2) = sd_baselines_sal; %Here you calculate average baseline per neuron

    amph_neuron_w = f2.DF(cell_to_index_map(y, 2),1:total_frames_water);
    if n_shock_Trials ~= 0
        amph_neuron_s = f2.DF(cell_to_index_map(y, 2),total_frames_water + 1:total_frames_water + total_frames_amph);
        sd_baselines_amph_s = std(amph_neuron_s);
        sd_amph_s(co, 1) = y; sd_amph_s(co, 2) = sd_baselines_amph_s;
    end
    %amph_neuron_w = amph_neuron_w';  % trial X frames
    %amph_neuron_s = amph_neuron_s';
    sd_baselines_amph_w = std(amph_neuron_w);
    
    sd_amph_w(co, 1) = y; sd_amph_w(co, 2) = sd_baselines_amph_w;

    co = co + 1;
end
sd_sal(:,2) = sd_sal(:,2)*2 * 2;
thresholds_sal = sd_sal;
sd_amph_w(:,2) = sd_amph_w(:,2)*2 * 2;
thresholds_amph_w = sd_amph_w;
sd_amph_s(:,2) = sd_amph_s(:,2)*2 * 2;
thresholds_amph_s = sd_amph_s;
%{
for y=cells_in_both
    b = 1;
    e = frames_per_trial;
    av = zeros(1,frames_per_trial); %holds sum of trials, frame by frame
    for i = (1:n_water_Trials)
        %for every trial,
        trial = sal_neuron(b:e);
        %Save avg baseline response per neuron per trial
        sal_per_trial(co, i) = mean(trial(window)); 
        av = av + trial;
        b = b + frames_per_trial;
        e = e + frames_per_trial;
    end
    %Add's each neuron's average baseline in saline
    average_baselines_sal = average_baselines_sal + mean(av(window)/n_water_Trials); 
    %Makes a list of all neuron's average baseline in saline
    average_neuron_baseline_sal = horzcat(average_neuron_baseline_sal, mean(av(window)/n_water_Trials));
    %Load amph data
    amph_neuron_w = f2.DF(cell_to_index_map(y, 2),1:total_frames_water);
    if n_shock_Trials ~= 0
        amph_neuron_s = f2.DF(cell_to_index_map(y, 2),total_frames_water + 1:total_frames_water + total_frames_amph);
    end
    b = 1;
    e = frames_per_trial;
    av2 = zeros(1,frames_per_trial);
    for i = (1:n_water_Trials)
        %for every trial,
        trial = amph_neuron_w(b:e);
        %Save avg baseline response per neuron per trial, amph, water
        amph_per_trial_w(co, i) = mean(trial(window));
        av2 = av2 + trial;
        b = b + frames_per_trial;
        e = e + frames_per_trial;
    end
    
    b = 1;
    e = frames_per_trial;
    av2_s = zeros(1,frames_per_trial);
    if n_shock_Trials ~= 0
        for i = (1:n_shock_Trials)
            trial = amph_neuron_s(b:e);
            amph_per_trial_s(co, i) = mean(trial(window));
            av2_s = av2_s + trial;
            b = b + frames_per_trial;
            e = e + frames_per_trial;
        end
    end
    %Add's each neuron's average baseline in amph water US
    average_baselines_amph_w = average_baselines_amph_w + mean(av2(window)/n_water_Trials);
    %Add's each neuron's average baseline in amph shock US
    if n_shock_Trials ~= 0  
        average_baselines_amph_s = average_baselines_amph_s + mean(av2_s(window)/n_shock_Trials);
        average_neuron_baseline_amph_s = horzcat(average_neuron_baseline_amph_s, mean(av2_s(window)/n_shock_Trials));
    end
    %Makes a list of all neuron's average baseline in amph
    average_neuron_baseline_amph_w = horzcat(average_neuron_baseline_amph_w, mean(av2(window)/n_water_Trials));
    
    co = co + 1;
end 
%}
%% Use the thresholds to count "events" for each neuron during baseline
events_sal = ones(size(cells_in_both, 2),1);
events_amph_w = ones(size(cells_in_both, 2),1);
events_amph_s = ones(size(cells_in_both, 2),1);
c = 1;
for o=cells_in_both   % neuron number    
    %% Saline water trials
    sal_index = cell_to_index_map(o, 1);
    sal_neuron = f1.DF(sal_index,:); %saline %pulls DF/F for each frame for 1 neuron
    p_temp2 = reshape(sal_neuron,FrameLength_water,n_water_Trials); %splits frames into trials
    p_temp3 = p_temp2';  % trial X frames
    p_temp4 = p_temp3(:,window); %looking only at baseline
    p_temp7 = p_temp4 > thresholds_sal(c, 2); %count number of times above threshold
    p_temp8 = p_temp7';
    long_events = reshape(p_temp8,1,num_window_frames*n_water_Trials); %conglomerate into 1D array 620 before
    events = 0;
      for x = 1:num_window_frames*n_water_Trials %count the number of times you pass threshold, 
          %not including times where you passed it the frame before.
          if x == 1
              if long_events(x) == 1
                  events = events + 1;      
              end
              continue
          end
          if and((long_events(x) == 1), (long_events(x - 1) == 0))
              events = events + 1;
          end
      end
      
    events_sal(c) = events;
    sum_events_sal = sum(sum(p_temp7));
    
    %% Amph water trials
    amph_neuron_w = f2.DF(cell_to_index_map(o, 2),1:total_frames_water);
    p_temp2 = reshape(amph_neuron_w,FrameLength_water,n_water_Trials); %splits frames into trials
    p_temp3 = p_temp2';  % trial X frames
    p_temp4 = p_temp3(:,window); %looking only at baseline
    p_temp7 = p_temp4 > thresholds_amph_w(c, 2); %count number of times above threshold
    p_temp8 = p_temp7';
    long_events = reshape(p_temp8,1,num_window_frames*n_water_Trials); %conglomerate into 1D array 620 before
    events = 0;
      for x = 1:num_window_frames*n_water_Trials %count the number of times you pass threshold, 
          %not including times where you passed it the frame before.
          if x == 1
              if long_events(x) == 1
                  events = events + 1;      
              end
              continue
          end
          if and((long_events(x) == 1), (long_events(x - 1) == 0))
              events = events + 1;
          end
      end
      
    events_amph_w(c) = events;
    sum_events_amph_w = sum(sum(p_temp7));
    
    %% Amph shock trial
    if n_shock_Trials ~= 0
        amph_neuron_s = f2.DF(cell_to_index_map(o, 2),4001:7000);
        p_temp2 = reshape(amph_neuron_s,FrameLength_water,n_shock_Trials); %splits frames into trials
        p_temp3 = p_temp2';  % trial X frames
        p_temp4 = p_temp3(:,window); %looking only at baseline
        p_temp7 = p_temp4 > thresholds_amph_s(c, 2); %count number of times above threshold
        p_temp8 = p_temp7';
        long_events = reshape(p_temp8,1,num_window_frames*n_shock_Trials); %conglomerate into 1D array 620 before
        events = 0;
          for x = 1:num_window_frames*n_shock_Trials %count the number of times you pass threshold, 
              %not including times where you passed it the frame before.
              if x == 1
                  if long_events(x) == 1
                      events = events + 1;      
                  end
                  continue
              end
              if and((long_events(x) == 1), (long_events(x - 1) == 0))
                  events = events + 1;
              end
          end

        events_amph_s(c) = events;
        c = c + 1;
        sum_events_amph_s = sum(sum(p_temp7));    
    end
end

events_sal = events_sal/n_water_Trials;
events_amph_w = events_amph_w/n_water_Trials;
if n_shock_Trials ~= 0
    events_amph_s = events_amph_s/n_shock_Trials;
end
% Now you need to use cell_to_index to match these cells
%{
cell_to_index_only_both = cell_to_index_map;
cell_to_index_only_both(horzcat(cells_only_in_sal, cells_only_in_amph), :) = []; % delete rows of cells only in one condition
ordered_events_amph_w = [];
for t = cell_to_index_only_both(:, 2)
    ordered_events_amph_w = horzcat(ordered_events_amph_w, events_amph_w(t)/n_water_Trials);
end

ordered_events_amph_s = [];
for t = cell_to_index_only_both(:, 2)
    ordered_events_amph_s = horzcat(ordered_events_amph_s, events_amph_s(t)/n_water_Trials);
end

ordered_events_sal = cell_to_index_only_both(:, 1);
ordered_events_sal = events_sal(ordered_events_sal)/n_water_Trials;
%}

% events_amph and events_sal are the numbers of events per neuron
disp('Mean events counted per trial')
disp(strcat('Neuron type: ', neuron_type));
disp(strcat('Mean sal baseline events:', num2str(mean(events_sal)/n_water_Trials)));
disp(strcat('Mean amph baseline events in water trials:', num2str(mean(events_amph_w)/n_water_Trials)));
if n_shock_Trials ~= 0
    disp(strcat('Mean amph baseline events in shock trials:', num2str(mean(events_amph_s)/n_shock_Trials)));
end

sall = zeros(size(events_sal, 1), 1);
ampp = zeros(size(events_amph_w, 1), 1) + 1;

events_sal_v = events_sal';
events_amph_v = events_amph_w';
size(events_amph_w, 1)

figure;
scatter(vertcat(sall, ampp), vertcat(events_sal, events_amph_w));
title('Event counts in Sal and Amph');

for u = 1:size(events_amph_w, 1)
    if events_sal(u) >= events_amph_w(u)
        disp(horzcat("Larger sal baseline in neuron # ", u));
    end
end 

figure;
hold on;
plot(mean(sal_per_trial),'LineWidth',8);
plot(mean(amph_per_trial_w),'LineWidth',8);
if n_shock_Trials ~= 0
    plot(mean(amph_per_trial_s),'LineWidth',8);
end
title('Average dF/F during baseline per trial');
xlabel('dF/F per Trial');
ylabel('Trial Number');
legend({'Saline Water','Amph Water', 'Amph Shock'},'Location','southwest')
hold off
%% Plots average baselines across all neurons over 3 seconds
y = figure;
p = uipanel('Parent',y,'BorderType','none'); 
p.Title ='Averaged baseline single neuron responses'; 
p.TitlePosition = 'centertop'; 
p.FontSize = 12;
p.FontWeight = 'bold';

subplot(1,1,1, 'Parent',p)       % add first plot in 2 x 1 grid
combo = vertcat(average_neuron_baseline_sal, average_neuron_baseline_amph_w);
one = plot(combo.');
legend({'Saline', 'Amphetamine'});
title('dF/F');
xlabel('Frame');
ylabel('dF/F magnitude');
set(one,'LineWidth',2);

all_baselines = horzcat(average_neuron_baseline_amph_s, average_neuron_baseline_amph_w, average_neuron_baseline_sal);
mean_all_baselines = mean(all_baselines);

average_neuron_baseline_amph_s_z = (average_neuron_baseline_amph_s - mean_all_baselines)/std(all_baselines);
average_neuron_baseline_amph_w_z = (average_neuron_baseline_amph_w - mean_all_baselines)/std(all_baselines);
average_neuron_baseline_sal_z = (average_neuron_baseline_sal - mean_all_baselines)/std(all_baselines);

%% Graph results
%Note: sal= neurons that appear both in saline and in amphetamine trials,
%during saline water trials. amph w = neurons present in both conditions
%during amph water trials. amph s = neurons present in both conditions
%during amph shock trials. Only saline- Neurons only present in saline
%condition, during saline water trials. only amph w- neurons only present
%in amph condition, during water trials. Amph s- neurons only present in
%amph condition, during shock trials. All neuron water- all neurons in amph
%condition during water. All neuron shock- all neuron in amph condition
%during shock.
if n_shock_Trials ~= 0
    figure
    hold on
    amplitudes = [mean(average_neuron_baseline_sal), ...
        mean(average_neuron_baseline_amph_w), ...
        mean(average_neuron_baseline_amph_s), ...
        mean(average_neuron_baseline_only_sal_w), ...
        mean(average_neuron_baseline_only_amph_w), ...
        mean(average_neuron_baseline_only_amph_s), ...
        mean(horzcat(average_neuron_baseline_only_amph_w, average_neuron_baseline_amph_w)), ...
        mean(horzcat(average_neuron_baseline_only_amph_s, average_neuron_baseline_amph_s))];
    bar(1:8, amplitudes);
    set(gca,'xticklabel',{''; 'sal';'amph w'; 'amph shock'; 'Only saline'; 'Only amph water'; 'Only amph shock'; 'All neuron water'; 'All neuron shock'})
    SEM = [std(average_neuron_baseline_sal)/sqrt(length(average_neuron_baseline_sal)), ...
        std(average_neuron_baseline_amph_w)/sqrt(length(average_neuron_baseline_amph_w)), ...
        std(average_neuron_baseline_amph_s)/sqrt(length(average_neuron_baseline_amph_s)), ...
        std(average_neuron_baseline_only_amph_w)/sqrt(length(average_neuron_baseline_only_amph_w)), ...
        std(average_neuron_baseline_only_sal_w)/sqrt(length(average_neuron_baseline_only_sal_w)), ...
        std(average_neuron_baseline_only_amph_s)/sqrt(length(average_neuron_baseline_only_amph_s)), ...
        std(horzcat(average_neuron_baseline_only_amph_w, average_neuron_baseline_amph_w))/sqrt(length(horzcat(average_neuron_baseline_only_amph_w, average_neuron_baseline_amph_w))), ...
        std(horzcat(average_neuron_baseline_only_amph_s, average_neuron_baseline_amph_s))/sqrt(length(horzcat(average_neuron_baseline_only_amph_s, average_neuron_baseline_amph_s)))];
    errorbar(1:8, amplitudes, SEM,'.')
    title(sprintf('Average dF/F of neurons appearing in both conditions during frames: %d-%d', window(1), window(length(window):length(window))));
else   
    figure
    hold on
    amplitudes = [mean(average_neuron_baseline_sal), mean(average_neuron_baseline_amph_w), mean(average_neuron_baseline_amph_s)];
    bar(1:2, amplitudes);
    set(gca,'xticklabel',{''; '';'sal';''; ''; ''; ''; 'amph w'; '';})
    SEM = [std(average_neuron_baseline_sal)/sqrt(length(average_neuron_baseline_sal)), std(average_neuron_baseline_amph_w)/sqrt(length(average_neuron_baseline_amph_w))];
    errorbar(1:2, amplitudes, SEM,'.')
    title(sprintf('Average dF/F during frames: %d-%d', window(1), window(length(window):length(window))));
end

%%%%%Plots baselines of cells in both
%{
baselines = vertcat(average_neuron_baseline_sal, average_neuron_baseline_amph_w);
g = zeros(2, 26);
g(1, :) = 1; %Make salines 1
g = g(:, 1)';

baselines = baselines';
sally = baselines(:, 1);
ampy = baselines(:, 2);
halp = vertcat(zeros(1, size(sally, 1))', (zeros(1, size(ampy, 1)) + 1)');
baselines = vertcat(baselines(:, 1), baselines(:, 2));
figure;
scatter(halp, baselines);
title('dF/F of saline and amphetamine')

for u = 1:size(average_neuron_baseline_amph_w, 2)
    if average_neuron_baseline_sal(u) > average_neuron_baseline_amph_w(u)
        disp(horzcat("Larger sal baseline in neuron # ", u));
    end
end 

end_amp_results = vertcat(average_neuron_baseline_sal, average_neuron_baseline_amph_w);
tags = [1, 2]';
end_amp_results = horzcat(tags, end_amp_results);
%}
