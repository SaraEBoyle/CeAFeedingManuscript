%% define parameter %%
%% Use this to find which neurons are active during CS/US delivery%%
%% Saves variables into response_pools.mat.
clear
files = { %It's important to put them in this order. sal amph sal amph
    '/Users/sboyle/Documents/Li Lab/Image analysis/SST16/20181115/saline1/', ...
    '/Users/sboyle/Documents/Li Lab/Image analysis/SST16/20181115/saline2/'%, ...
    %'/Users/sboyle/Documents/Li Lab/Image analysis/SST6/onlysaline/', ...
    %'/Users/sboyle/Documents/Li Lab/Image analysis/SST6/amph/'
    };

%Set these
nWaterTrials = 20;
nShockTrials = 0;%20 for water, 15 for shock
seconds_per_trial = 18;
baseline = [49 79]; %og= 69 99 3 seconds right before the stimulus
shock_us = [106 140]; %100-125 is the CS-- tone. 126-150,160 is the US-- water/shock
shock_cs = [80 105]; %100-125 is the CS-- tone. 126-150,160 is the US-- water/shock
water_us = [106 140]; %100-125 is the CS-- tone. 126-150,160 is the US-- water/shock
water_cs = [80 105]; %100-125 is the CS-- tone. 126-150,160 is the US-- water/shock

NUM_ITER = 184756;
framerate = 10;
saline_frames_per_session = framerate*seconds_per_trial*nWaterTrials;
amph_frames_per_session = framerate*seconds_per_trial*(nWaterTrials + nShockTrials);
amph_water_frames_per_session = framerate*seconds_per_trial*nWaterTrials;
amph_shock_frames_per_session = framerate*seconds_per_trial*nShockTrials;
FrameLength_sal = saline_frames_per_session/nWaterTrials;
FrameLength_amph_water = amph_water_frames_per_session/nWaterTrials;
FrameLength_amph_shock = amph_shock_frames_per_session/nShockTrials;

q = 1;
for x = files%([2 4])
    response_pool = [];
    p_val = [];
    disp(x);
    load(fullfile(cell2mat(x), 'Sorted.mat'));
    %% First do water cs responses %%
    DF_new = DF(:,1:saline_frames_per_session);%Change this value to control which frames. 4001-7000 shock, 1-4000 water
    fpath = sprintf('%sresponses', cell2mat(x));
    mkdir (fpath);
    total_neurons = size(DF_new,1);
    % calculate permutation test and p value %%
    for o=1:size(DF_new,1)   % neuron number
          p_temp = DF_new(o,:); %pulls DF/F for each frame for 1 neuron
          p_temp2 = reshape(p_temp,FrameLength_sal,nWaterTrials); %splits frames into trials
          p_temp3 = p_temp2';  % trial X frames
          p_temp4 = p_temp3(:,water_cs(1):water_cs(2));   % calcualte average vuale within stimulus period
          p_temp5 = mean(p_temp4,2);
          positive = p_temp5;
          p_temp7 = p_temp3(:,baseline(1):baseline(2));     % calcualte average vuale within baseline period
          p_temp8 = mean(p_temp7,2);
          negtive = p_temp8;
          mean_val = mean(positive) - mean(negtive);
          lenx_p = length(positive);
          leny_p = length(negtive);
          p = [negtive;positive];
          p_average = zeros (NUM_ITER,1);
          for s = 1: NUM_ITER
                 inds= randperm(lenx_p+leny_p);
                shuff_po = p(inds(1:lenx_p));
                shuff_ne = p(inds(1+lenx_p:end));
                p_average (s,1) = mean(shuff_po)-mean(shuff_ne);
          end
          p_sort = sort(p_average);
          p_val(o) = sum(p_sort>mean_val)/ NUM_ITER ;
    end
    response_number = sum(p_val<0.05);
    w=1;
    for i=1:size(p_val,2)
        if(p_val(i)<0.05)
            response_pool(w) = i;
            pool(w) = p_val(i);
            w=w+1;
        end
    end
    for i=1: response_number
        j = response_pool(i);
        neuron_data = DF_new(j,:);
        neuron_data = neuron_data';
        neuron_data1 = reshape(neuron_data, FrameLength_sal,nWaterTrials);
        neuron_data1 = neuron_data1';
        h(i) = figure, imagesc(neuron_data1);
        title(['Neuron # ',num2str(j),' C']);
        answer = questdlg('Does this look real?', ...
        'Choose CS responsive neurons', ...
        'yes','no','halp', 'yes');
        if strcmpi(answer, 'yes')
            saveas(h(i),fullfile(fpath, sprintf('water_cs_set_%s_neuron_%s.jpg',num2str(q), num2str(j))));
        elseif strcmpi(answer, 'no') 
            response_pool(i) = 0;
        else
            disp('Sorry buddy. I quit the program for you.')
            close all;
            return
        end
        close(h(i));
    end
    response_pool(response_pool == 0) = [];
    water_cs_response_pool = response_pool; 
    Test1_p_val = p_val;
    %% Then do the water US
    %Reset some variables
    response_pool = [];
    p_val = [];
    
    % calculate permutation test and p value %%
    for o=1:size(DF_new,1)   % neuron number
          p_temp = DF_new(o,:); %pulls DF/F for each frame for 1 neuron
          p_temp2 = reshape(p_temp,FrameLength_sal,nWaterTrials); %splits frames into trials
          p_temp3 = p_temp2';  % trial X frames
          p_temp4 = p_temp3(:,water_us(1):water_us(2));   % calcualte average vuale within stimulus period
          p_temp5 = mean(p_temp4,2);
          positive = p_temp5;
          p_temp7 = p_temp3(:,baseline(1):baseline(2));     % calcualte average vuale within baseline period
          p_temp8 = mean(p_temp7,2);
          negtive = p_temp8;
          mean_val = mean(positive) - mean(negtive);
          lenx_p = length(positive);
          leny_p = length(negtive);
          p = [negtive;positive];
          p_average = zeros (NUM_ITER,1);
          for s = 1: NUM_ITER
                 inds= randperm(lenx_p+leny_p);
                shuff_po = p(inds(1:lenx_p));
                shuff_ne = p(inds(1+lenx_p:end));
                p_average (s,1) = mean(shuff_po)-mean(shuff_ne);
          end
          p_sort = sort(p_average);
          p_val(o) = sum(p_sort>mean_val)/ NUM_ITER ;
    end

    response_number = sum(p_val<0.05);
    w=1;
    for i=1:size(p_val,2)
        if(p_val(i)<0.05)
            response_pool(w) = i;
            pool(w) = p_val(i);
            w=w+1;
        end
    end
    for i=1: response_number
        j = response_pool(i);
        neuron_data = DF_new(j,:);
        neuron_data = neuron_data';
        neuron_data1 = reshape(neuron_data, FrameLength_sal,nWaterTrials);
        neuron_data1 = neuron_data1';
        h(i) = figure, imagesc(neuron_data1);
        title(['Neuron # ',num2str(j),' C']);
        answer = questdlg('Does this look real?', ...
        'Choose US responsive neurons', ...
        'yes','no','halp', 'yes');
% Handle response
        if strcmpi(answer, 'yes')
            saveas(h(i),fullfile(fpath, sprintf('water_us_set_%s_neuron_%s.jpg',num2str(q), num2str(j))));
        elseif strcmpi(answer, 'halp')
            disp('You just quit the program.')
            close all;
            return
        else
            response_pool(i) = 0;
        end
        close(h(i));
    end
    response_pool(response_pool == 0) = [];
    water_us_response_pool = response_pool; 
    Test1_p_val = p_val; 
    %% Then do the shock CS
    %Reset some variables
    if nShockTrials > 0
        response_pool = [];
        p_val = [];
        DF_new = DF(:,amph_water_frames_per_session + 1:amph_frames_per_session);%Change this value to control which frames. 4001-7000 shock, 1-4000 water
        % calculate permutation test and p value %%
        for o=1:size(DF_new,1)   % neuron number
              p_temp = DF_new(o,:); %pulls DF/F for each frame for 1 neuron
              p_temp2 = reshape(p_temp,FrameLength_amph_shock,nShockTrials); %splits frames into trials
              p_temp3 = p_temp2';  % trial X frames
              p_temp4 = p_temp3(:,shock_cs(1):shock_cs(2));   % calcualte average vuale within stimulus period
              p_temp5 = mean(p_temp4,2);
              positive = p_temp5;
              p_temp7 = p_temp3(:,baseline(1):baseline(2));     % calcualte average vuale within baseline period
              p_temp8 = mean(p_temp7,2);
              negtive = p_temp8;
              mean_val = mean(positive) - mean(negtive);
              lenx_p = length(positive);
              leny_p = length(negtive);
              p = [negtive;positive];
              p_average = zeros (NUM_ITER,1);
              for s = 1: NUM_ITER
                     inds= randperm(lenx_p+leny_p);
                    shuff_po = p(inds(1:lenx_p));
                    shuff_ne = p(inds(1+lenx_p:end));
                    p_average (s,1) = mean(shuff_po)-mean(shuff_ne);
              end
              p_sort = sort(p_average);
              p_val(o) = sum(p_sort>mean_val)/ NUM_ITER ;
        end
        response_number = sum(p_val<0.05);
        w=1;
        for i=1:size(p_val,2)
            if(p_val(i)<0.05)
                response_pool(w) = i;
                pool(w) = p_val(i);
                w=w+1;
            end
        end
        for i=1: response_number
            j = response_pool(i);
            neuron_data = DF_new(j,:);
            neuron_data = neuron_data';
            neuron_data1 = reshape(neuron_data, FrameLength_amph_shock,nShockTrials);
            neuron_data1 = neuron_data1';
            h(i) = figure, imagesc(neuron_data1);
            title(['Neuron # ',num2str(j),' C']);
            answer = questdlg('Does this look real?', ...
            'Choose CS responsive neurons', ...
            'yes','no','halp', 'yes');
    % Handle response
            if strcmpi(answer, 'yes')
                saveas(h(i),fullfile(fpath, sprintf('shock_cs_set_%s_neuron_%s.jpg',num2str(q), num2str(j))));
            else
                response_pool(i) = 0;
            end
            close(h(i));
        end
        response_pool(response_pool == 0) = [];
        shock_cs_response_pool = response_pool; 
       %% Then do the shock US
        %Reset some variables
        response_pool = [];
        p_val = [];

        % calculate permutation test and p value %%
        for o=1:size(DF_new,1)   % neuron number
              p_temp = DF_new(o,:); %pulls DF/F for each frame for 1 neuron
              p_temp2 = reshape(p_temp,FrameLength_amph_shock,nShockTrials); %splits frames into trials
              p_temp3 = p_temp2';  % trial X frames
              p_temp4 = p_temp3(:,shock_us(1):shock_us(2));   % calcualte average vuale within stimulus period
              p_temp5 = mean(p_temp4,2);
              positive = p_temp5;
              p_temp7 = p_temp3(:,baseline(1):baseline(2));     % calcualte average vuale within baseline period
              p_temp8 = mean(p_temp7,2);
              negtive = p_temp8;
              mean_val = mean(positive) - mean(negtive);
              lenx_p = length(positive);
              leny_p = length(negtive);
              p = [negtive;positive];
              p_average = zeros (NUM_ITER,1);
              for s = 1: NUM_ITER
                     inds= randperm(lenx_p+leny_p);
                    shuff_po = p(inds(1:lenx_p));
                    shuff_ne = p(inds(1+lenx_p:end));
                    p_average (s,1) = mean(shuff_po)-mean(shuff_ne);
              end
              p_sort = sort(p_average);
              p_val(o) = sum(p_sort>mean_val)/ NUM_ITER ;
        end
        response_number = sum(p_val<0.05);
        w=1;
        for i=1:size(p_val,2)
            if(p_val(i)<0.05)
                response_pool(w) = i;
                pool(w) = p_val(i);
                w=w+1;
            end
        end
        for i=1: response_number
            j = response_pool(i);
            neuron_data = DF_new(j,:);
            neuron_data = neuron_data';
            neuron_data1 = reshape(neuron_data, FrameLength_amph_shock,nShockTrials);
            neuron_data1 = neuron_data1';
            h(i) = figure, imagesc(neuron_data1);
            title(['Neuron # ',num2str(j),' C']);
            answer = questdlg('Does this look real?', ...
            'Choose US responsive neurons', ...
            'yes','no','I have no idea what I am doing', 'yes');
    % Handle response
            if strcmpi(answer, 'yes')
                saveas(h(i),fullfile(fpath, sprintf('shock_us_set_%s_neuron_%s.jpg',num2str(q), num2str(j))));
            else
                response_pool(i) = 0;
            end
            close(h(i));
        end
        response_pool(response_pool == 0) = [];
        shock_us_response_pool = response_pool; 
    else
        shock_us_response_pool = [];
        shock_cs_response_pool = [];
    end
    %% Save response pools
    o = sprintf('%sresponse_pools', cell2mat(x));
    save(o, 'water_cs_response_pool', 'water_us_response_pool', ...
        'shock_us_response_pool', 'shock_cs_response_pool', 'total_neurons');
    q = q + 1;
end

