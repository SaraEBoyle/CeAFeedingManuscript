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
%First choose the file from running cell reg
%Commented is for SST6
load('/Users/sboyle/Documents/Li Lab/Image analysis/SST5/cell_registration/cellRegistered.mat');
%load('/Users/sboyle/Documents/Li Lab/Image analysis/SST5and6onlywatertrials/Cell_reg/Testing_SST5/cellRegistered_20180905_105057');
%load('/Users/sboyle/Documents/Li Lab/Image analysis/SST5and6onlywatertrials/Cell_reg/Testing_SST6/cellRegistered_20180905_192938');
%% Choose the correct response pool- these report US and CS responses
response_pools_1 = load('/Users/sboyle/Documents/Li Lab/Image analysis/SST5/onlysaline/response_pools.mat');
response_pools_2 = load('/Users/sboyle/Documents/Li Lab/Image analysis/SST5/amph/response_pools.mat');
%response_pools_3 = load('/Users/sboyle/Documents/Li Lab/Image analysis/SST5and6onlywatertrials/Cell_reg/Testing_response_pool/response_pools_3.mat');
%response_pools_4 = load('/Users/sboyle/Documents/Li Lab/Image analysis/SST5and6onlywatertrials/Cell_reg/Testing_response_pool/response_pools_4.mat');
pools = [response_pools_1, response_pools_2];

%% Define variables
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
%3/4 amph water US neurons also respond to shock
%5 water US sal neurons are in both conditions, and 
disp(horzcat('cells that respond to shock CS in amph: ', num2str(shock_cs_responders_amph)));
disp(horzcat('cells that respond to shock US in amph: ', num2str(shock_us_responders_amph)));

disp(horzcat(['Cells in both saline and amphetamine trials: ' newline 'Total: ', num2str(size(cells_in_both, 2)) newline], num2str(cells_in_both)));
disp(horzcat(['Cells that respond to water CS in both: ' newline 'Total: ', num2str(size(water_cs_responders_both, 2)), ' out of ', num2str(size(cells_in_both, 2)), ' total cells in both' newline], num2str(water_cs_responders_both)));
disp(horzcat(['Cells that respond to water US in both: ' newline 'Total: ', num2str(size(water_us_responders_both, 2))]));
for n=water_us_responders_both
    disp(horzcat('Saline and amphetamine indeces of US responders: ', '[', num2str(cell_to_index_map(n, [1,2])), ']'));    
end


%First row- Is it a shock cs responder?
both_cell_stats_1 = ismember(cells_in_both, shock_cs_responders_amph);
%Second row- Is it a shock us responder?
both_cell_stats_2 = ismember(cells_in_both, shock_us_responders_amph);
%Third row- Is it a water cs responder?
both_cell_stats_3 = ismember(cells_in_both, water_cs_responders_amph);
%Fourth row- Is it a water us responder?
both_cell_stats_4 = ismember(cells_in_both, water_us_responders_amph);
%fifth row- water cs sal
both_cell_stats_5 = ismember(cells_in_both, water_cs_responders_sal);
%Sixth row- Is it a water us sal responder?
both_cell_stats_6 = ismember(cells_in_both, water_us_responders_sal);

both_cell_stats = vertcat(both_cell_stats_1, ...
    both_cell_stats_2, ...
    both_cell_stats_3, ...
    both_cell_stats_4, ...
    both_cell_stats_5, ...
    both_cell_stats_6);

%2 neurons respond to both shock US and water US, but one of the two only
%responds in saline trials, not amph...

%% Graph the neurons responsive to US in both conditions
files = {'/Users/sboyle/Documents/Li Lab/Image analysis/SST5/onlysaline/', ...
    '/Users/sboyle/Documents/Li Lab/Image analysis/SST5/amph/', ...
    '/Users/sboyle/Documents/Li Lab/Image analysis/SST6/onlysaline/', ...
    '/Users/sboyle/Documents/Li Lab/Image analysis/SST6/amph/'
    };

% TODO: MAKE THIS AUTOMATICALLY PRINT THE CONSISTANT US RESPONSIVE ONES
%% Graph all neurons that respond consistently to US or CS in both conditions
f1 = (files(1));
f1 = load(fullfile(cell2mat(f1), 'Sorted.mat'));
f2 = (files(2));
f2 = load(fullfile(cell2mat(f2), 'Sorted.mat'));
if ~isempty(water_us_responders_both)       %If there are US responders in both
    for cells = water_us_responders_both
        sal_ind = cell_to_index_map(n, 1);
        amph_ind = cell_to_index_map(n, 2);
        neuron_s = f1.DF(sal_ind,:); %This is a US responsive neuron from both conditions
        %You need to average over the 20 trials and then plot that graph
        b = 1;
        e = 200;
        av = zeros(1,200);
        for i = (1:20)
            av = av + neuron_s(b:e);
            b = b + 200;
            e = e + 200;
        end
        %graph1 = plot(av);

        neuron_a = f2.DF(amph_ind,:);
        %You need to average over the 20 trials and then plot that graph
        b = 1;
        e = 200;
        av2 = zeros(1,200);
        for i = (1:20)
            av2 = av2 + neuron_a(b:e);
            b = b + 200;
            e = e + 200;
        end

        f = figure;
        p = uipanel('Parent',f,'BorderType','none'); 
        p.Title = horzcat('Neuron ', num2str(sal_ind), ' US Response'); 
        p.TitlePosition = 'centertop'; 
        p.FontSize = 12;
        p.FontWeight = 'bold';

        %Code for graphing normalized neuron
        %combo = vertcat(av/max(av), av2/max(av2));
        %title(horzcat('Neuron ', num2str(us_responders_both)));
        %subplot(2,1,1, 'Parent',p)       % add first plot in 2 x 1 grid
        %plot(combo.');
        %legend({'Saline', 'Amphetamine'});
        %title('Normalized')

        subplot(1,1,1, 'Parent',p)       
        plot(vertcat(av, av2).')      
        legend({'Saline', 'Amphetamine'});
        title('dF/F');
        line([100 100], [0 1.8]);
        line([125 125], [0 1.8]);
    end %Graph sal + amph response
end
if ~isempty(water_cs_responders_both)       %If there are CS responders in both
    for cells = water_cs_responders_both
        sal_ind = cell_to_index_map(n, 1);
        amph_ind = cell_to_index_map(n, 2);
        neuron_s = f1.DF(sal_ind,:); %This is a CS responsive neuron from both conditions
        %You need to average over the 20 trials and then plot that graph
        b = 1;
        e = 200;
        av = zeros(1,200);
        for i = (1:20)
            av = av + neuron_s(b:e);
            b = b + 200;
            e = e + 200;
        end
        %graph1 = plot(av);

        neuron_a = f2.DF(amph_ind,:);
        %You need to average over the 20 trials and then plot that graph
        b = 1;
        e = 200;
        av2 = zeros(1,200);
        for i = (1:20)
            av2 = av2 + neuron_a(b:e);
            b = b + 200;
            e = e + 200;
        end

        f = figure;
        p = uipanel('Parent',f,'BorderType','none'); 
        p.Title = horzcat('Neuron ', num2str(sal_ind), ' CS Response'); 
        p.TitlePosition = 'centertop'; 
        p.FontSize = 12;
        p.FontWeight = 'bold';

        %Code for graphing normalized neuron
        %combo = vertcat(av/max(av), av2/max(av2));
        %title(horzcat('Neuron ', num2str(us_responders_both)));
        %subplot(2,1,1, 'Parent',p)       % add first plot in 2 x 1 grid
        %plot(combo.');
        %legend({'Saline', 'Amphetamine'});
        %title('Normalized')

        subplot(1,1,1, 'Parent',p)       
        plot(vertcat(av, av2).')      
        legend({'Saline', 'Amphetamine'});
        title('dF/F');
        line([100 100], [0 1.8]);
        line([125 125], [0 1.8]);
    end %Graph sal + amph response
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

% Calculates average dF/F over 3 second baseline period in both conditions
%Record the average dF/F per neuron per trial
sal_per_trial = [];
amph_per_trial_w = [];
amph_per_trial_s = [];
co = 1;
%{
for y=cells_in_both
    sal_neuron = f1.DF(y,:); %saline
    %You need to average over the 20 trials and then plot that graph
    b = 1;
    e = 200;
    av = zeros(1,200);
    for i = (1:20)
        sal_per_trial(co, i) = mean(sal_neuron(b:e));
        av = av + sal_neuron(b:e);
        b = b + 200;
        e = e + 200;
    end
    average_baselines_sal_n = average_baselines_sal_n + av(69:99)/max(av); %Normalized by dividing by largest response from each neuron
    average_baselines_sal = average_baselines_sal + av(69:99);
    average_neuron_baseline_sal = horzcat(average_neuron_baseline_sal, mean(av));
    amph_neuron_w = f2.DF(cell_to_index_map(y, 2),1:4000); %amphetamine
    amph_neuron_s = f2.DF(cell_to_index_map(y, 2),4001:7000);
    %You need to average over the 20 trials and then plot that graph
    b = 1;
    e = 200;
    av2 = zeros(1,200);
    for i = (1:20)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        amph_per_trial_w(co, i) = mean(amph_neuron_w(b:e));
        av2 = av2 + amph_neuron_w(b:e);
        b = b + 200;
        e = e + 200;
    end
    
    b = 1;
    e = 200;
    av2_s = zeros(1,200);
    for i = (1:15)
        amph_per_trial_s(co, i) = mean(amph_neuron_s(b:e));
        av2_s = av2_s + amph_neuron_s(b:e);
        b = b + 200;
        e = e + 200;
    end
    average_baselines_amph_n = average_baselines_amph_n + av2(69:99)/max(av2); %Normalized by dividing by largest response from each neuron
    average_baselines_amph_w = average_baselines_amph_w + av2(69:99);
    average_baselines_amph_s = average_baselines_amph_s + av2_s(69:99);
    average_neuron_baseline_amph_w = horzcat(average_neuron_baseline_amph_w, mean(av2));
    average_neuron_baseline_amph_s = horzcat(average_neuron_baseline_amph_s, mean(av2_s));
    co = co + 1;
end %Calculates average baseline in frames 69-99
%}
%UNCOMMENT IF YOU DON'T WANT TO PLOT TRIAL BY TRIAL
counter = 1;

for y=horzcat(cells_in_both,cells_only_in_sal)
    sal_neuron = f1.DF(y,:); %saline
    %You need to average over the 20 trials and then plot that graph
    b = 1;
    e = 200;
    av = zeros(1,200);
    for i = (1:20)
        trial = sal_neuron(b:e);
        av = av + trial;
        sal_per_trial(counter, i) = mean(trial(69:99));
        b = b + 200;
        e = e + 200;
    end
    %average_baselines_sal_n = average_baselines_sal_n + av(69:99)/max(av); %Normalized by dividing by largest response from each neuron
    average_baselines_sal = average_baselines_sal + av(69:99);
    average_neuron_baseline_sal = horzcat(average_neuron_baseline_sal, mean(av(69:99)));
    counter = counter + 1;
    %You need to average over the 20 trials and then plot that graph
    % UNCOMMENT
end

counter = 1;
for y=horzcat(cells_in_both,cells_only_in_amph)
    amph_neuron_w = f2.DF(cell_to_index_map(y, 2),1:4000); %amphetamine
    amph_neuron_s = f2.DF(cell_to_index_map(y, 2),4001:7000);
    b = 1;
    e = 200;
    av2 = zeros(1,200);
    for i = (1:20)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        trial = amph_neuron_w(b:e);
        amph_per_trial_w(counter, i) = mean(trial(69:99));
        av2 = av2 + trial;
        b = b + 200;
        e = e + 200;
    end
    b = 1;
    e = 200;
    av2_s = zeros(1,200);
    for i = (1:15)
        trial = amph_neuron_s(b:e);
        amph_per_trial_s(counter, i) = mean(trial(69:99));
        av2_s = av2_s + trial;
        b = b + 200;
        e = e + 200;
    end
    average_baselines_amph_n = average_baselines_amph_n + av2(69:99)/max(av2); %Normalized by dividing by largest response from each neuron
    average_baselines_amph_w = average_baselines_amph_w + av2(69:99);
    average_baselines_amph_s = average_baselines_amph_s + av2_s(69:99);
    average_neuron_baseline_amph_w = horzcat(average_neuron_baseline_amph_w, mean(av2(69:99)));
    average_neuron_baseline_amph_s = horzcat(average_neuron_baseline_amph_s, mean(av2_s(69:99)));
    counter = counter + 1;
end 
% 20 columns, 48 rows

%average_baselines_sal = average_baselines_sal/20;
%average_baselines_amph_w = average_baselines_amph_w/20;
%average_baselines_amph_s = average_baselines_amph_s/15;
%average_baselines_sal_n = average_baselines_sal_n/20;
%average_baselines_amph_n = average_baselines_amph_n/20;
%sal_per_trial is 38 neurons by 20 trials
disp(mean(mean(sal_per_trial)));
disp(mean(mean(amph_per_trial_w)));
disp(mean(mean(amph_per_trial_s)));
figure;
hold on
plot(mean(sal_per_trial));
plot(mean(amph_per_trial_w));
plot(mean(amph_per_trial_s));
title('Average dF/F during baseline per trial');
xlabel('dF/F per Trial');
ylabel('Trial Number');
legend({'Saline Water','Amph Water', 'Amph Shock'},'Location','southwest')
%plot(horzcat(mean(sal_per_trial), mean(amph_per_trial_w), mean(amph_per_trial_s))); 