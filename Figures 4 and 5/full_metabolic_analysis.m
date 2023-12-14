function [] = full_metabolic_analysis()
%% Go through all mice and highlight obvious RER artifacts. Replace these values with
%average data from nearby to preserve time info. Display data and record
%variables to easily plot in prism. Has the option to smooth RER data and
%downsample to make statistical analysis easier (so you don't have to calculate an 
%ANOVA on 1,500 time points). Saves full_graphing_variables.mat for line
%graphing of RER, down_sampled_graphing_variables.mat for downsampled RER (recommended),
%full_trace_data.mat for all raw heat, O2, CO2, food, and drink traces, and 
%full_metabolic_matrix.mat for avg values segmented by light and dark
%periods. Also saves RERParameters.mat to keep track of the settings. x
%axis variables are called "time_align_first_dark" and contain timing data
%with the first dark cycle set as 0. If you see 'ds' within a variable
%name, that means there was downsampling. 'smooth' means there's been
%smoothing on that variable. Something that says "smooth_ds_smooth" means
%there was smoothing, then downsampling, then mild smoothing. All smoothed
%variables are plotted so you can adjust level of smoothing.

% inputs: a .csv file, saved directly from CLAX software, with data from all cages
%   click data, check all, and wait. It should be a large csv. Then click
%   export and save it as all.csv.
%   You also need to specify the exp (experimental group) and con (control
%   group) variables to indicate which cages are which. There is a
%   to_exclude varible which you can specify to prevent plotting of empty
%   cages
%>> Written by Sara Boyle, saraboyle2017@u.northwestern.edu, 2022, Matlab 2017a    
    close all
    
    %% Variables for line plotting
    % if you want no smoothing or downsampling, you can set these all to 1. I'd
    % recommend a small amount of smoothing- I think it makes average
    % calculations less noisy. You can adjust these parameters and look at
    % how your data changes, then choose your favorite
    
    rf = 5; % reduction factor for how much to downsample. Makes statistical analysis easier
    smooth_span = 15; % how big a window to smooth RER, by moving average. 15 would be ~1 hour, since each point is about 4.3 minutes apart
    ds_smooth_span = 3;% how big a window to smooth after down sampling, I recommend 3
    lowest_y = .65; % specify y limit for RER plotting
    highest_y = 1;
    heat_lowest_y = .2;
    heat_highest_y =.75;
    
    %% Variables for average calculations
    % Use this to specify which data to include in averages. Everything is
    % aligned to the first dark cycle, since the first light cycle will be
    % a variable length. To include time before that, use negative numbers
    
    time_from = 0; %start of time span in hours, 0 = first dark cycle, -5 = 5 hours before that
    time_to = 60; % end of time span for avg calculations
    
    %% specify which cages are your experimental or control group
    
    %% B2 NaCh
    %exp = [6, 4, 10, 11]; 
    %con = [5, 3, 9, 12];
    %to_exclude = []; % if you have empty or terribly wonky cages, exclude for easier artifact screening (only for viewing)
    %% B1 NaCh
    %exp = [2, 3, 5, 6, 7]; %B1
    %con = [9, 10, 11, 12];
    %to_exclude = [];
    %% B2 from Ale
    %exp = [4, 5, 8, 9, 10]; %experimental group
    %con = [6, 7, 11]; %control group
    %to_exclude = [1, 2, 3, 12]; % cages to exclude (empty or malfunctioning cages)
    %% B3 TeLC from Ale
    exp = [1, 3, 4, 6, 7, 10]; %experimental group
    con = [2, 5, 8, 9, 11]; %control group
    to_exclude = [12]; % cages to exclude (empty or malfunctioning cages)
  
    %% DON'T CHANGE BELOW HERE
    time_to_plot = 1000; %minutes to plot at one time while skimming data for artifacts
    save('RERParameters.mat', 'rf', 'smooth_span', 'ds_smooth_span', 'time_from', ...
        'time_to', 'exp', 'con', 'to_exclude');
    
    %% Import ALL data
    %% TODO, just get the all csv working because it's waaaaaaayyyyyy easier
    [doubles, b, cells] = xlsread('all.csv');
    %segment into individual cages
    cage_inds = find(doubles(:, 1) == 1);
    cage_inds = [cage_inds;size(doubles,1) + 1];
    %segment into different data types
    
    raw_vo2 = doubles(:, 4);
    raw_do2 = doubles(:, 7);
    raw_vco2 = doubles(:, 9);
    raw_dco2 = doubles(:, 12);
    raw_RERs = doubles(:, 14);
    raw_heat = doubles(:, 15);
    raw_foods = doubles(:, 19);
    raw_drinks = doubles(:, 21);
    raw_x_amb = doubles(:, 23); % OH NO THIS IS 23
    raw_x_tot = doubles(:, 22); % THIS IS 22
    raw_y_amb = doubles(:, 25); %25
    raw_y_tot = doubles(:, 24); %24
    raw_z_tot = doubles(:, 26); 
    
    vo2 = [];
    do2 = [];
    vco2 = [];
    dco2 = [];
    heat = [];
    foods = [];
    drinks = [];
    x_amb = [];
    x_tot = [];
    y_amb = [];
    y_tot = [];
    z_tot = [];
    RERs = [];
    bigness = cage_inds(2) - cage_inds(1);
    for ind = 1:(size(cage_inds, 1) - 1)
        RERtest = raw_RERs(cage_inds(ind):(cage_inds(ind + 1) - 1));
        vo2test = raw_vo2(cage_inds(ind):(cage_inds(ind + 1) - 1));
        do2test = raw_do2(cage_inds(ind):(cage_inds(ind + 1) - 1));
        vco2test = raw_vco2(cage_inds(ind):(cage_inds(ind + 1) - 1));
        dco2test = raw_dco2(cage_inds(ind):(cage_inds(ind + 1) - 1));
        heattest = raw_heat(cage_inds(ind):(cage_inds(ind + 1) - 1));
        foodstest = raw_foods(cage_inds(ind):(cage_inds(ind + 1) - 1));
        drinkstest = raw_drinks(cage_inds(ind):(cage_inds(ind + 1) - 1));
        x_ambtest = raw_x_amb(cage_inds(ind):(cage_inds(ind + 1) - 1));
        x_tottest = raw_x_tot(cage_inds(ind):(cage_inds(ind + 1) - 1));
        y_ambtest = raw_y_amb(cage_inds(ind):(cage_inds(ind + 1) - 1));
        y_tottest = raw_y_tot(cage_inds(ind):(cage_inds(ind + 1) - 1));
        z_tottest = raw_z_tot(cage_inds(ind):(cage_inds(ind + 1) - 1));
        while size(RERtest, 1) < bigness
            RERtest = [RERtest(1);RERtest];
            vo2test = [vo2test(1);vo2test];
            do2test = [do2test(1);do2test];
            vco2test = [vco2test(1);vco2test];
            dco2test = [dco2test(1);dco2test];
            heattest = [heattest(1);heattest];
            foodstest = [foodstest(1);foodstest];
            drinkstest = [drinkstest(1);drinkstest];
            x_ambtest = [x_ambtest(1);x_ambtest];
            x_tottest = [x_tottest(1);x_tottest];
            y_ambtest = [y_ambtest(1);y_ambtest];
            y_tottest = [y_tottest(1);y_tottest];
            z_tottest = [z_tottest(1);z_tottest];
        end
        RERs(:, ind) = RERtest;
        vo2(:, ind) = vo2test;
        do2(:, ind) = do2test;
        vco2(:, ind) = vco2test;
        dco2(:, ind) = dco2test;
        heat(:, ind) = heattest;
        foods(:, ind) = foodstest;
        drinks(:, ind) = drinkstest;
        x_amb(:, ind) = x_ambtest;
        x_tot(:, ind) = x_tottest;
        y_amb(:, ind) = y_ambtest;
        y_tot(:, ind) = y_tottest;
        z_tot(:, ind) = z_tottest;
    end
    
    %19
    %21
    %% Read in timing data
    %[a, b, c] = xlsread('cage1.csv');
    dates = cells(2:(bigness + 1), 3);  
    %% Read in food data
    %[foods, b, trash] = xlsread('food.csv');
    %foods = foods(2:end, 2:end);  
    
    % drink
    %[drinks, b, trash] = xlsread('drink.csv');
    %drinks = drinks(2:end, 2:end);
    
    to_plot = setdiff([exp, con], to_exclude);
    % find timing data relative to the beginning value
    first_point = datevec(dates{1});
    time_vector_min = [0];
    prev_point = 0;
    for num = 2:size(dates, 1)
        cur_point=datevec(dates{num});
        time_interval_min = etime(cur_point,first_point)/60; % find time interval in minutes
        time_vector_min = [time_vector_min, time_interval_min + prev_point];
        prev_point = time_interval_min + prev_point;
        first_point = cur_point;
    end
    disp('Data loaded');
    
    %% Manually flag artifacts and record cleaned and smoothed data
    % read in RER data, sort out timing
    %[RERs, bb, cc] = xlsread('RERs.csv');
    %RERs(:,1) = []; % delete the stupid
    times_isnan = sum(isnan(RERs), 2);
    isnan_inds = find(times_isnan > 0);
    to_delete_x = [time_vector_min(isnan_inds), time_vector_min(isnan_inds)];
    %Plot each 500 segment chunk of all mice, select periods with obvious
    %artifacts, delete or fill
    min_length = time_vector_min(end);
    screen = get(0,'ScreenSize');
   
    if exist('to_delete_RER.mat', 'file')
        choice = input('WARNING: You already selected artifacts to delete. Do you want to redo this? input y to redo or press enter to use old settings: ', 's');
        if ~strcmp(choice, 'y')
            disp('Analyzing with saved artifact information...')
            load('to_delete_RER.mat');       
        else
            disp('Manually removing artifacts...')
            color_array = [1 0 0; 0 1 0; 0 1 1; 0 0 1; 1 1 0; 0 0 0; 0 0.4470 0.7410; ...
                0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; ...
                0.4660 0.6740 0.1880; 0.6350 0.0780 0.1840];
            for time = 0:time_to_plot:(min_length - time_to_plot)
                start_ind = find(time_vector_min >= time);
                start_ind = start_ind(1);
                end_ind = find(time_vector_min < time + time_to_plot);
                end_ind = end_ind(end);
                plot_time = time_vector_min(start_ind: end_ind);
                figure
                hold on
                set(gcf,'position',[100,200,2500,450])
                xlabel('minutes');

                for cage_id = to_plot
                    plot(plot_time, RERs(start_ind:end_ind, cage_id), 'Color', color_array(cage_id, :));
                end
                set(gcf,'position',[100,100,screen(3) - 200,screen(4)/3]);
                title('Click on left, then right side of each artifact to delete. When finished with whole graph, press enter');
                [x_in,y_in]=ginput(); % specify number of points if you can
                to_delete_x = [to_delete_x, x_in'];
                close all
            end
            save('to_delete_RER.mat', 'to_delete_x', 'to_plot');
        end
    else
        disp('Manually removing artifacts...')
            color_array = [1 0 0; 0 1 0; 0 1 1; 0 0 1; 1 1 0; 0 0 0; 0 0.4470 0.7410; ...
                0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560; ...
                0.4660 0.6740 0.1880; 0.6350 0.0780 0.1840];
            for time = 0:time_to_plot:(min_length - time_to_plot)
                start_ind = find(time_vector_min >= time);
                start_ind = start_ind(1);
                end_ind = find(time_vector_min < time + time_to_plot);
                end_ind = end_ind(end);
                plot_time = time_vector_min(start_ind: end_ind);
                figure
                hold on
                set(gcf,'position',[100,200,2500,450])
                xlabel('minutes');

                for cage_id = to_plot
                    plot(plot_time, RERs(start_ind:end_ind, cage_id), 'Color', color_array(cage_id, :));
                end
                set(gcf,'position',[100,100,screen(3) - 200,screen(4)/3]);
                title('Click on left, then right side of each artifact to delete. When finished with whole graph, press enter');
                [x_in,y_in]=ginput(); % specify number of points if you can
                to_delete_x = [to_delete_x, x_in'];
                close all
            end
            save('to_delete_RER.mat', 'to_delete_x', 'to_plot');
    end
    
    %% Get indeces of light and dark cycles
    cycle = cells(2:(bigness + 1), 4);
    lights = strfind(cycle, 'Light');
    darks = strfind(cycle, 'Dark');
    yournumber=0;
    darks = cellfun(@(x) x*yournumber,darks,'un',0);
    yournumber=1;
    lights = cellfun(@(x) x*yournumber,lights,'un',0);
    cycles = [];
    for ind = 1:size(darks, 1)
        d = darks{ind};
        L = lights{ind};
        if isempty(d)
            d = 0;
        end
        if isempty(L)
            L = 0;
        end
        cycles(ind) = d + L;
    end    
    
    disp('Now replacing artifacts with previous values...');
    all_del_inds = [];
    for ind = 1:2:(size(to_delete_x, 2) - 1)
        first = to_delete_x(ind);
        last = to_delete_x(ind + 1);
        first_ind = find(time_vector_min >= first);
        last_ind = find(time_vector_min < last);
        del_ind = intersect(first_ind, last_ind);
        if last == time_vector_min(end)
            %if end of file has NaN
            del_ind = [del_ind, first_ind];
        end
        all_del_inds = [all_del_inds, del_ind];
    end
    
    pre = figure;
    plot(time_vector_min, RERs(:, to_plot));
    ylim([lowest_y, highest_y]);
    set(gcf,'position',[100,150 + screen(4)/3,screen(3) - 200,screen(4)/3]);
    title('Pre deleted artifacts')
    %% Instead of deleting, fill with average of previous 5 values. Preserves timing info
    prev = 0;
    prev_avg = 0;
    for del = all_del_inds
        if abs(del - prev) > 1
            if del - 5 < 0
                too_smol = del;
                % find the next ind that's not to be deleted and take avg
                % after that
                while ismember(too_smol, all_del_inds)
                    too_smol = too_smol + 1;
                end
                five_avg = RERs((too_smol):(too_smol + 5), :);
                five_avg = mean(five_avg);
            else
                five_avg = RERs((del - 5):(del - 1), :);
                five_avg = mean(five_avg);
            end
            
            RERs(del, :) = five_avg;
            prev = del;
            prev_avg = five_avg;
        else
            RERs(del, :) = prev_avg;
            prev = del;
        end
    end
    
    light_inds = find(cycles == 1);
    dark_inds = find(cycles == 0);
    post = figure;
    ylim([lowest_y, highest_y]);
    plot(time_vector_min, RERs(:, to_plot));
    set(gcf,'position',[100,100,screen(3) - 200,screen(4)/3]);
    title('Post deleted artifacts'); 
    input('Press enter to confirm')
    
    %% Smooth and/or down sample RERs using chosen parameters
    smoothed_RERs = [];
    for i = 1:size(RERs, 2)
        %Just smooth using moving average, window size = smooth_span
        sub = smooth(time_vector_min, RERs(:, i), smooth_span, 'moving');
        smoothed_RERs(:, i) =  sub;
        
        % for this, downsample the already smoothed RERs
        down_samp_smoothed_RERs(:, i) = downsample(sub, rf);
        
        %Downsample data that has been smoothed using the smooth_span, then
        %smooth that again using a mild smoothing (ds_smooth_span), just to reduce
        %jerkiness of downsampling
        smoothed_down_samp_smoothed_RERs(:, i) = smooth(downsample(sub, rf), ds_smooth_span, 'moving');
        
        %Also downsample light cycle data in the same way
        ds_cycles = downsample(cycles, rf);
    end
   
    %% Align RER data to the first dark dark cycle, convert to hrs for readability
    time_vector_hr = time_vector_min./60;
    down_samp_time_vector_hr = downsample(time_vector_hr, rf);
    first_dark = find(cycles == 0);
    first_dark = first_dark(1);
    time_align_1st_dark_hr = time_vector_hr - time_vector_hr(first_dark);
    
    close all
    
    %% Get start and end times of dark cycles
    ds_dark_cycles = ~ds_cycles;
    dark_cycles = ~cycles;
    %get all inds included in dark cycle
    all_dark_inds = find(dark_cycles == 1);
    % first ind will be start of first cycle
    dark_starts = [all_dark_inds(1)];
    
    dark_ends = [];
    prev = all_dark_inds(1);
    for p = all_dark_inds(2:end)
        cur = p;
        if cur - prev > 10
            dark_starts = [dark_starts, cur];
            dark_ends = [dark_ends, prev];
        end
        prev = cur;
    end
    dark_ends = [dark_ends, all_dark_inds(end)];
    
    ds_time_align_1st_dark_hr = downsample(time_align_1st_dark_hr, rf);
    
    %save('RERData.mat', 'RERs', 'dark_cycles', 'ds_dark_cycles', 'smoothed_RERs','ds_time_align_1st_dark_hr', 'smoothed_down_samp_smoothed_RERs', 'down_samp_smoothed_RERs', 'ds_cycles', 'down_samp_time_vector_hr', 'light_inds', 'dark_inds', 'time_vector_min', 'time_vector_hr', 'exp', 'con', 'cycles');
    
    %% Use only times in specified interval for average calculations. No downsampling
    afters = find(time_align_1st_dark_hr >= time_from); %get indeces after interval
    befores = find(time_align_1st_dark_hr <= time_to); %get indeces before interval
    insides = intersect(befores, afters); %get indeces within time window
    
    inside_lights = intersect(light_inds, insides); %get light indeces in window
    inside_darks = intersect(dark_inds, insides); %get light indeces in window
    
    %RER
    RER_avg_light = smoothed_RERs(inside_lights, :);
    RER_avg_dark = smoothed_RERs(inside_darks, :);
    
    %Heat
    heat_avg_light = heat(inside_lights, :);
    heat_avg_dark = heat(inside_darks, :);
    
    end_food = find(time_align_1st_dark_hr >= time_to);
    if isempty(end_food)
        disp('Your session is not long enough. Reduce time_to variable');
    end
    end_food_ind = end_food(1);
    
    start_food = find(time_align_1st_dark_hr <= time_from);
    start_food_ind = start_food(end);
    food_highest_y = round(max(max(foods))) + 3; % max y value for food graph
    
    total_food = foods(end_food_ind, :) - foods(start_food_ind, :);
    total_drink = drinks(end_food_ind, :) - drinks(start_food_ind, :);
    total_x_amb = sum(x_amb(start_food_ind:end_food_ind, :));
    total_x_amb_light = sum(x_amb(inside_lights, :));
    total_x_amb_dark = sum(x_amb(inside_darks, :));
    
    total_y_amb = sum(y_amb(start_food_ind:end_food_ind, :));
    total_y_amb_light = sum(y_amb(inside_lights, :));
    total_y_amb_dark = sum(y_amb(inside_darks, :));
    
    total_z = sum(z_tot(start_food_ind:end_food_ind, :));
    total_z_light = sum(z_tot(inside_lights, :));
    total_z_dark = sum(z_tot(inside_darks, :));
    
    total_x = sum(x_tot(start_food_ind:end_food_ind, :));
    total_x_light = sum(x_tot(inside_lights, :));
    total_x_dark = sum(x_tot(inside_darks, :));
    
    total_y = sum(y_tot(start_food_ind:end_food_ind, :));
    total_y_light = sum(y_tot(inside_lights, :));
    total_y_dark = sum(y_tot(inside_darks, :));
    
    total_ambulatory_bb = total_x_amb + total_y_amb;
    total_ambulatory_bb_light = total_x_amb_light + total_y_amb_light;
    total_ambulatory_bb_dark = total_x_amb_dark + total_y_amb_dark;
    
    total_non_ambulatory_bb = total_x - total_x_amb + total_y - total_y_amb;
    total_non_ambulatory_bb_light = total_x_light - total_x_amb_light + total_y_light - total_y_amb_light;
    total_non_ambulatory_bb_dark = total_x_dark - total_x_amb_dark +  total_y_dark - total_y_amb_dark;
    
    %% Save animal ids
    full_metabolic_matrix = {'Experimental Group';'Cage ID';'RER';'Light';'Dark';'Total';'Heat'; ...
        'Light';'Dark';'Total';'vO2';'Light';'Dark';'Total';'vCO2';'Light';'Dark'; ...
        'Total';'dO2';'Light'; 'Dark';'Total';'dCO2';'Light';'Dark';'Total'; ...
        'Food';'Total';'Water';'Total';'X beam breaks';'Light'; 'Dark'; 'Total'; ...
        'X ambulatory beam breaks';'Light';'Dark';'Total';'Y beam breaks'; 'Light'; ...
        'Dark';'Total';'Y ambulatory beam breaks';'Light'; ...
        'Dark';'Total';'Z beam breaks';'Light';'Dark';'Total'; ...
        'Total Ambulatory Beam Breaks';'Light';'Dark';'Total'; ...
        'Total Non-Ambulatory Beam Breaks';'Light';'Dark';'Total'};
    to_include = [sort(exp), sort(con)];
    for y = 2:(1 + size(to_include, 2))
        %% TODO save data in a better format add each evg measure in dark and light
        id = to_include(y - 1);
        if ismember(id, exp)
            full_metabolic_matrix{1, y} = 'Experimental';
        else
            full_metabolic_matrix{1, y} = 'Control';
        end
        full_metabolic_matrix{2, y} = id;
        %Input avg RERs
        full_metabolic_matrix{4, y} = mean(RER_avg_light(:, id));
        full_metabolic_matrix{5, y} = mean(RER_avg_dark(:, id));
        full_metabolic_matrix{6, y} = mean(smoothed_RERs([inside_lights,inside_darks], id));
        %Input avg heats
        full_metabolic_matrix{8, y} = mean(heat_avg_light(:, id));
        full_metabolic_matrix{9, y} = mean(heat_avg_dark(:, id));
        full_metabolic_matrix{10, y} = mean(heat([inside_lights,inside_darks], id));
        
        %Input avg vO2
        full_metabolic_matrix{12, y} = mean(vo2(inside_lights, id));
        full_metabolic_matrix{13, y} = mean(vo2(inside_darks, id));
        full_metabolic_matrix{14, y} = mean(vo2([inside_lights,inside_darks], id));
        
        %Input avg vCO2
        full_metabolic_matrix{16, y} = mean(vco2(inside_lights, id));
        full_metabolic_matrix{17, y} = mean(vco2(inside_darks, id));
        full_metabolic_matrix{18, y} = mean(vco2([inside_lights,inside_darks], id));
        
        %Input avg dO2
        full_metabolic_matrix{20, y} = mean(do2(inside_lights, id));
        full_metabolic_matrix{21, y} = mean(do2(inside_darks, id));
        full_metabolic_matrix{22, y} = mean(do2([inside_lights,inside_darks], id));
        
        %Input avg dCO2
        full_metabolic_matrix{24, y} = mean(dco2(inside_lights, id));
        full_metabolic_matrix{25, y} = mean(dco2(inside_darks, id));
        full_metabolic_matrix{26, y} = mean(dco2([inside_lights,inside_darks], id));
        
        %Input food
        full_metabolic_matrix{28, y} = total_food(:, id);
        
        %Input water
        full_metabolic_matrix{30, y} = total_drink(:, id);
        
        %Input avg x beam breaks
        full_metabolic_matrix{32, y} = total_x_light(:, id);
        full_metabolic_matrix{33, y} = total_x_dark(:, id);
        full_metabolic_matrix{34, y} = total_x(:, id);
        
        %input x ambulatory beam breaks
        full_metabolic_matrix{36, y} = total_x_amb_light(:, id);
        full_metabolic_matrix{37, y} = total_x_amb_dark(:, id);
        full_metabolic_matrix{38, y} = total_x_amb(:, id);
        
        %Input avg y beam breaks
        full_metabolic_matrix{40, y} = total_y_light(:, id);
        full_metabolic_matrix{41, y} = total_y_dark(:, id);
        full_metabolic_matrix{42, y} = total_y(:, id);
        
        %input y ambulatory beam breaks
        full_metabolic_matrix{44, y} = total_y_amb_light(:, id);
        full_metabolic_matrix{45, y} = total_y_amb_dark(:, id);
        full_metabolic_matrix{46, y} = total_y_amb(:, id);
        
        %Input avg z beam breaks
        full_metabolic_matrix{48, y} = total_z_light(:, id);
        full_metabolic_matrix{49, y} = total_z_dark(:, id);
        full_metabolic_matrix{50, y} = total_z(:, id);
        
        % Input ambulatory beam breaks
        full_metabolic_matrix{52, y} = total_ambulatory_bb_light(:, id);
        full_metabolic_matrix{53, y} = total_ambulatory_bb_dark(:, id);
        full_metabolic_matrix{54, y} = total_ambulatory_bb(:, id);
        
        % Input non ambulatory beam breaks
        full_metabolic_matrix{56, y} = total_non_ambulatory_bb_light(:, id);
        full_metabolic_matrix{57, y} = total_non_ambulatory_bb_dark(:, id);
        full_metabolic_matrix{58, y} = total_non_ambulatory_bb(:, id);
        
    end
    save('full_metabolic_matrix.mat', 'full_metabolic_matrix');
    % RER Get averages of experimental and controls in light and dark cycles
    exp_RER_avg_light = mean(RER_avg_light(:, exp));
    con_RER_avg_light = mean(RER_avg_light(:, con));
    exp_RER_avg_dark = mean(RER_avg_dark(:, exp));
    con_RER_avg_dark = mean(RER_avg_dark(:, con));
    %heat_matrix(4, (1 + exp)) = exp_RER_avg_light;
    %Heat
    exp_heat_avg_light = mean(heat_avg_light(:, exp));
    con_heat_avg_light = mean(heat_avg_light(:, con));
    exp_heat_avg_dark = mean(heat_avg_dark(:, exp));
    con_heat_avg_dark = mean(heat_avg_dark(:, con));
    
    %% Save convenient variables for graphing
    smooth_ds_smooth_exp_RERs = smoothed_down_samp_smoothed_RERs(:, exp);
    smooth_ds_smooth_con_RERs = smoothed_down_samp_smoothed_RERs(:, con);
    avg_smooth_ds_smooth_exp_RERs = mean(smoothed_down_samp_smoothed_RERs(:, exp), 2);
    avg_smooth_ds_smooth_con_RERs = mean(smoothed_down_samp_smoothed_RERs(:, con), 2);
    
    smooth_exp_RERs = smoothed_RERs(:, exp);
    smooth_con_RERs = smoothed_RERs(:, con);
    avg_smooth_exp_RERs = mean(smoothed_RERs(:, exp), 2);
    avg_smooth_con_RERs = mean(smoothed_RERs(:, con), 2);
    ds_dark_cycles = int8(ds_dark_cycles);
    save('down_sampled_graphing_variables.mat', 'smooth_ds_smooth_exp_RERs', 'smooth_ds_smooth_con_RERs', 'ds_dark_cycles', 'ds_time_align_1st_dark_hr');
    save('full_graphing_variables.mat', 'smooth_exp_RERs', 'smooth_con_RERs', 'dark_cycles', 'time_align_1st_dark_hr');
    ds_dark_cycles = logical(ds_dark_cycles);
    
    %% Create function for plotting individual lines
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot individual smoothed RER lines
    params = {exp, con, lowest_y, highest_y, dark_starts, dark_ends, time_from, ...
        time_to, screen, smooth_span, time_align_1st_dark_hr};
    plot_individual_lines(time_align_1st_dark_hr, smoothed_RERs, params);
    %bottom
    set(gcf,'position',[100,100,screen(3) - 200,screen(4)/3]);
    title(['Smoothed RERs with window size: ', num2str(smooth_span), ', no down sampling']);
    savefig('smoothed_RER_individuals.fig');
    
    %% Plot individual smoothed RER lines with downsampling
    plot_individual_lines(ds_time_align_1st_dark_hr, smoothed_down_samp_smoothed_RERs, params);
    %top
    set(gcf,'position',[100,150 + screen(4)/3,screen(3) - 200,screen(4)/3]);
    title(['Smoothed RERs with window size: ' num2str(smooth_span), ', with down sampling factor: ', num2str(rf)]);
    savefig('ds_smoothed_RER_individuals.fig');
    input('Press enter to continue');
    close all
    %% Plot individual heat lines with no downsampling
    heat_params = {exp, con, heat_lowest_y, heat_highest_y, dark_starts, dark_ends, time_from, ...
        time_to, screen, smooth_span, time_align_1st_dark_hr};
    plot_individual_lines(time_align_1st_dark_hr, heat, heat_params);
    set(gcf,'position',[100,150 + screen(4)/3,screen(3) - 200,screen(4)/3]);
    title('Heat with no down sampling or smoothing');
    savefig('heat_individuals.fig');
    
    food_params = {exp, con, 0, food_highest_y, dark_starts, dark_ends, time_from, ...
        time_to, screen, smooth_span, time_align_1st_dark_hr};
    
    plot_individual_lines(time_align_1st_dark_hr, foods, food_params);
    %set(gcf,'position',[100,150 + screen(4)/3,screen(3) - 200,screen(4)/3]);
    title('Food with no down sampling or smoothing');
    %ylim([0, 10]);
    savefig('food_individuals.fig');
   input('Press enter to confirm');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all
    gray = [0.94902, 0.94902, 0.94902];
    %% Plot average smoothed and downsampled curves
    figure
    y1 = ones(1, size(time_align_1st_dark_hr, 2))*highest_y;
    y2 = ones(1, size(time_align_1st_dark_hr, 2))*lowest_y;
    con_err = std(smooth_ds_smooth_con_RERs')/sqrt(size(smooth_ds_smooth_con_RERs,2));
    AreaPlot(ds_time_align_1st_dark_hr,avg_smooth_ds_smooth_con_RERs',smooth(con_err),gray,0.1,1);
    hold on
    con_line = plot(ds_time_align_1st_dark_hr, mean(smoothed_down_samp_smoothed_RERs(:, con), 2), 'k');
    
    exp_err = std(smooth_ds_smooth_exp_RERs')/sqrt(size(smooth_ds_smooth_exp_RERs,2));
    AreaPlot(ds_time_align_1st_dark_hr,avg_smooth_ds_smooth_exp_RERs',smooth(exp_err),[0.80000, 1.00000, 0.80000],0.1,1);
    
    % Plot ds dark cycles
    for day = 1:size(dark_ends, 2)
        first_patch = zeros(1, size(time_align_1st_dark_hr, 2));
        first_patch(dark_starts(day):dark_ends(day)) = 1;
        first_patch = first_patch == 1;
        P = patch([time_align_1st_dark_hr(first_patch) fliplr(time_align_1st_dark_hr(first_patch))], [y1(first_patch), fliplr(y2(first_patch))], [0.90196, 0.90196, 1.00000], 'FaceAlpha', 0.5, 'LineStyle', 'none');    
    end
    
    exp_line = plot(ds_time_align_1st_dark_hr, mean(smoothed_down_samp_smoothed_RERs(:, exp), 2), 'g');
    set(gcf,'position',[100,175 + screen(4)/3,screen(3) - 200,screen(4)/3.2]);
    title(['Averaged smoothed RERs, with downsampling factor: ', num2str(rf), ' and standard error bars']);
    ylim([lowest_y, highest_y]);
    legend([exp_line, con_line, P], {'Exp', 'Con', 'Dark'});
    savefig('Avg_RERs.fig')
    
    %% Plot averages with no downsampling
    figure
    plot(time_align_1st_dark_hr, avg_smooth_exp_RERs, 'g')
    hold on
    plot(time_align_1st_dark_hr, avg_smooth_con_RERs, 'k')
    set(gcf,'position',[100,100,screen(3) - 200,screen(4)/3.2]);
    title('Averaged smoothed RERs, without downsampling');
    % Plot ds dark cycles
    for day = 1:size(dark_ends, 2)
        first_patch = zeros(1, size(time_align_1st_dark_hr, 2));
        first_patch(dark_starts(day):dark_ends(day)) = 1;
        first_patch = first_patch == 1;
        P = patch([time_align_1st_dark_hr(first_patch) fliplr(time_align_1st_dark_hr(first_patch))], [y1(first_patch), fliplr(y2(first_patch))], [0.90196, 0.90196, 1.00000], 'FaceAlpha', 0.5, 'LineStyle', 'none');    
    end
    ylim([lowest_y, highest_y]);
    legend([exp_line, con_line, P], {'Exp', 'Con', 'Dark'});
    savefig('Avg_RERs_no_ds.fig')
   
    save('full_trace_data.mat', 'heat', 'vo2', 'vco2', 'do2', 'dco2', 'foods', ...
        'drinks', 'time_align_1st_dark_hr', 'dark_cycles', 'exp', 'con')
 
end



