function [] = calculate_z_score_w_consistent_baselinesV2()
% Previously I calculate z score using baselines that are different for
% each stimulus. Now I will pool all calculated baselines to have more
% comparable z scores.
    start_AUC = 0;
    end_AUC = 5;
    delete_negative_AUC = 1;
    %% Load signal events
    left_trials = [];
    right_trials = [];
    file = horzcat(pwd, '\');
    tem = dir(horzcat(file, '*', 'signal')); 
    for signal_folder_ind = 1:size(tem, 1)
        signal_folder = tem(signal_folder_ind, :);
        name = signal_folder.name;
        children = dir(horzcat(file, name, '\')); 
        all_baselines = [];
        for child_ind = 3:size(children, 1)
            child = children(child_ind);
            child_name = child.name;
            if or(contains(child_name, 'Entire'), contains(child_name, 'pearson'))
                continue
            elseif contains(child_name, 'full_session')
                continue
            elseif contains(child_name, 'raw_data')
                continue
            end
            if exist(horzcat(file, name, '\', child_name, '\dFF_data.mat'))
                load(horzcat(file, name, '\', child_name, '\dFF_data.mat'));
                all_baselines = horzcat(all_baselines, reshape(baselines, [1, (size(baselines, 1)*size(baselines, 2))]));
            else
                child_children = dir(horzcat(file, name, '\', child_name));
                for u = 3:4
                    grandkid = child_children(u);
                    og_file = horzcat(file, name, '\', child_name, '\');
                    load(horzcat(og_file, grandkid.name, '\dFF_data.mat'));
                    all_baselines = horzcat(all_baselines, reshape(baselines, [1, (size(baselines, 1)*size(baselines, 2))]));
                end
            end
            
        end
        %Calculate std_avg for all data
        std_avg = std(all_baselines);
        avg = mean(all_baselines);
        
        grandkids_included = children;
        for child_ind = 3:size(children, 1)
            grandkid = dir(horzcat(file, name, '\', children(child_ind).name));
            grandkids_included = vertcat(grandkids_included, grandkid);
        end

        for child_ind = 3:size(grandkids_included, 1)
            US_trials = [];
            child = grandkids_included(child_ind);
            child_name = child.name;
            child_folder = child.folder;
            if or(contains(child_name, 'Entire'), contains(child_name, 'pearson'))
                continue
            elseif contains(child_name, 'full_session')
                continue
            elseif contains(child_name, 'raw_data')
                continue
            elseif strcmp(child_name, '.')
                continue
            elseif strcmp(child_name, '..')
                continue
            elseif ~exist(horzcat(child_folder, '\', child_name, '\dFF_data.mat'))
                continue
            end
            load(horzcat(child_folder, '\', child_name, '\dFF_data.mat'));
            if ~isempty(left_trials)
                US_trials_percent = left_trials;
                left_trials = [];
            elseif ~isempty(right_trials)
                US_trials_percent = right_trials;
                right_trials = [];
            elseif ~isempty(US_trials)
                US_trials_percent = US_trials;
                US_trials = [];
            end
            diffs = (US_trials_percent - avg);
            z_US_trials = (diffs ./ std_avg);
            cur_file = horzcat(child_folder, '\', child_name);
            if exist(horzcat(cur_file, '\artifacts_removed'), 'dir')
                load(horzcat(cur_file, '\', 'artifacts_removed\deleted_trials.mat'));
                z_US_trials(to_remove, :) = [];
            end
            %mean_z_US_trials = mean(z_US_trials);
            
            if size(z_US_trials, 1) > 1
                mean_z_US_trials = mean(z_US_trials);
            else
                mean_z_US_trials = z_US_trials; 
            end
            lefts = find(x_axis >= start_AUC);
            rights = find(x_axis <= end_AUC);
            zero_thru_five = intersect(lefts, rights);
            zero_thru_five_axis = x_axis(zero_thru_five);
            
            AUCs = [];
            for m = 1:size(z_US_trials, 1)
              
                current_trial = z_US_trials(m, :);
                current_AUC = current_trial(zero_thru_five);
                current_axis = zero_thru_five_axis;
                secs = 1;
                if delete_negative_AUC
                    zero_inds = find(current_AUC < 0); % delete the ones below zero
                    current_AUC(zero_inds) = [];
                    current_axis(zero_inds) = [];
                    secs = size(current_AUC, 2)/10;
                end
                if and(~isempty(current_axis), size(current_axis, 2) > 1)
                    AUCs = horzcat(AUCs, trapz(current_axis, current_AUC/secs));
                else
                    AUCs = horzcat(AUCs, 0);
                end
            end
            %% TODO FIGURE SOMETHING OUT
            mean_AUC_trial_by_trial = mean(AUCs);
            zero_thru_five_mean_z = mean_z_US_trials(zero_thru_five);
            if delete_negative_AUC
                zero_inds = find(zero_thru_five_mean_z < 0); % delete the ones below zero
                zero_thru_five_axis(zero_inds) = [];
                zero_thru_five_mean_z(zero_inds) = [];
            end
            zero_thru_five_AUC = trapz(zero_thru_five_axis, zero_thru_five_mean_z);
            
            save(horzcat(child_folder, '\', child_name, '\consistent_z_scores.mat'), 'z_US_trials', 'mean_z_US_trials', 'zero_thru_five_AUC', 'mean_AUC_trial_by_trial', 'AUCs');
        end
    end
    
    %% Now for control
    tem = dir(horzcat(file, '*', 'control')); 
    for signal_folder_ind = 1:size(tem, 1)
        signal_folder = tem(signal_folder_ind, :);
        name = signal_folder.name;
        children = dir(horzcat(file, name, '\')); 
        all_baselines = [];
        for child_ind = 3:size(children, 1)
            child = children(child_ind);
            child_name = child.name;
            if or(contains(child_name, 'Entire'), contains(child_name, 'pearson'))
                continue
            end
            %load(horzcat(file, name, '\', child_name, '\dFF_data.mat'));
            
            
            
            
            if exist(horzcat(file, name, '\', child_name, '\dFF_data.mat'))
                load(horzcat(file, name, '\', child_name, '\dFF_data.mat'));
                all_baselines = horzcat(all_baselines, reshape(baselines, [1, (size(baselines, 1)*size(baselines, 2))]));
            else
                child_children = dir(horzcat(file, name, '\', child_name));
                for u = 3:4
                    grandkid = child_children(u);
                    og_file = horzcat(file, name, '\', child_name, '\');
                    load(horzcat(og_file, grandkid.name, '\dFF_data.mat'));
                    all_baselines = horzcat(all_baselines, reshape(baselines, [1, (size(baselines, 1)*size(baselines, 2))]));
                end
            end
            
            
            
            
            
            
            
            
           % all_baselines = horzcat(all_baselines, reshape(baselines, [1, (size(baselines, 1)*size(baselines, 2))]));
        end
        %Calculate std_avg for all data
        std_avg = std(all_baselines);
        avg = mean(all_baselines);
        
        grandkids_included = children;
        for child_ind = 3:size(children, 1)
            grandkid = dir(horzcat(file, name, '\', children(child_ind).name));
            grandkids_included = vertcat(grandkids_included, grandkid);
        end

        for child_ind = 3:size(grandkids_included, 1)
            US_trials = [];
            child = grandkids_included(child_ind);
            child_name = child.name;
            child_folder = child.folder;
            if or(contains(child_name, 'Entire'), contains(child_name, 'pearson'))
                continue
            elseif contains(child_name, 'full_session')
                continue
            elseif contains(child_name, 'raw_data')
                continue
            elseif strcmp(child_name, '.')
                continue
            elseif strcmp(child_name, '..')
                continue
            elseif ~exist(horzcat(child_folder, '\', child_name, '\dFF_data.mat'))
                continue
            end
            load(horzcat(child_folder, '\', child_name, '\dFF_data.mat'));
            if ~isempty(left_trials)
                US_trials_percent = left_trials;
                left_trials = [];
            elseif ~isempty(right_trials)
                US_trials_percent = right_trials;
                right_trials = [];
            elseif ~isempty(US_trials)
                US_trials_percent = US_trials;
                US_trials = [];
            end
            diffs = (US_trials_percent - avg);
            z_US_trials = (diffs ./ std_avg);
            fake_name = strrep(child_folder, 'control', 'signal');
            cur_file = horzcat(fake_name, '\', child_name);
            if exist(horzcat(cur_file, '\artifacts_removed'), 'dir')
                load(horzcat(cur_file, '\', 'artifacts_removed\deleted_trials.mat'));
                z_US_trials(to_remove, :) = [];
            end
            if size(z_US_trials, 1) > 1
                mean_z_US_trials = mean(z_US_trials);
            else
                mean_z_US_trials = z_US_trials; 
            end
            lefts = find(x_axis >= start_AUC);
            rights = find(x_axis <= end_AUC);
            zero_thru_five = intersect(lefts, rights);
            zero_thru_five_axis = x_axis(zero_thru_five);
            %zero_thru_five_mean_z = mean_z_US_trials(zero_thru_five);
            AUCs = [];
            
            for m = 1:size(z_US_trials, 1)
              
                current_trial = z_US_trials(m, :);
                current_AUC = current_trial(zero_thru_five);
                current_axis = zero_thru_five_axis;
                
                if delete_negative_AUC
                    zero_inds = find(current_AUC < 0); % delete the ones below zero
                    current_AUC(zero_inds) = [];
                    current_axis(zero_inds) = [];
                end
                if and(~isempty(current_axis), (size(current_axis, 2) ~= 1))
                    AUCs = horzcat(AUCs, trapz(current_axis, current_AUC));
                else 
                    AUCS = horzcat(AUCs, 0);
                end
            end
          
            
%            zero_thru_five_AUC = trapz(zero_thru_five_axis, zero_thru_five_mean_z);
            save(horzcat(child_folder, '\', child_name, '\consistent_z_scores.mat'), 'z_US_trials', 'mean_z_US_trials', 'zero_thru_five_AUC');
        end
    end
    
end
