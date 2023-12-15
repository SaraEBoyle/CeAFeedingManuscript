function [] = calculate_z_score_w_consistent_baselinesV2()
% Previously I calculate z score using baselines that are different for
% each stimulus. Now I will pool all calculated baselines to have more
% comparable z scores.

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
            elseif contains(child_name, 'STOP')
                continue
            elseif contains(child_name, 'drops')
                continue
            elseif contains(child_name, 'z_scores_whole_session.mat')
                continue
            end
            load(horzcat(file, name, '\', child_name, '\dFF_data.mat'));
            all_baselines = horzcat(all_baselines, reshape(baselines, [1, (size(baselines, 1)*size(baselines, 2))]));
        end
        %Calculate std_avg for all data
        std_avg = std(all_baselines);
        avg = mean(all_baselines);

        for child_ind = 3:size(children, 1)
            child = children(child_ind);
            child_name = child.name;
            if or(contains(child_name, 'Entire'), contains(child_name, 'pearson'))
                continue
            elseif contains(child_name, 'full_session')
                continue
            elseif contains(child_name, 'STOP')
                continue
            elseif contains(child_name, 'z_scores_whole_session.mat')
                continue
            elseif contains(child_name, 'drops')
                continue
            end
            load(horzcat(file, name, '\', child_name, '\dFF_data.mat'));
            if ~isempty(left_trials)
                US_trials_percent = left_trials;
                left_trials = [];
            elseif ~isempty(right_trials)
                US_trials_percent = right_trials;
                right_trials = [];
            end
            diffs = (US_trials_percent - avg);
            z_US_trials = (diffs ./ std_avg);
            cur_file = horzcat(file, name, '\', child_name);
            if exist(horzcat(cur_file, '\artifacts_removed'), 'dir')
                load(horzcat(cur_file, '\', 'artifacts_removed\deleted_trials.mat'));
                if ~isempty(to_remove)
                if to_remove(end) <= size(z_US_trials, 1)
                    z_US_trials(to_remove, :) = [];
                end
                end
                to_remove = [];
                
            end
            mean_z_US_trials = mean(z_US_trials);
            save(horzcat(file, name, '\', child_name, '\consistent_z_scores.mat'), 'z_US_trials', 'mean_z_US_trials');
        end
        load(horzcat(file, name, '\Corrected_data_full_session.mat'));
        %gcamp = cell2mat(raw_signal_dFFs);
        new_z_score = gcamp(:, 2);
        new_z_score = (new_z_score' - (avg/100))/(std_avg/100);
        gcamp(:, 2) = new_z_score';
        z_scores_whole_session = gcamp;
        plot(z_scores_whole_session(:, 1), z_scores_whole_session(:, 2));
        save(horzcat(file, name, '\z_scores_whole_session.mat'), 'z_scores_whole_session');
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
            elseif contains(child_name, 'STOP')
                continue
            elseif contains(child_name, 'z_scores_whole_session.mat')
                continue
            elseif contains(child_name, 'drops')
                continue
            end
            load(horzcat(file, name, '\', child_name, '\dFF_data.mat'));
            all_baselines = horzcat(all_baselines, reshape(baselines, [1, (size(baselines, 1)*size(baselines, 2))]));
        end
        %Calculate std_avg for all data
        std_avg = std(all_baselines);
        avg = mean(all_baselines);

        for child_ind = 3:size(children, 1)
            child = children(child_ind);
            child_name = child.name;
            if or(contains(child_name, 'Entire'), contains(child_name, 'pearson'))
                continue
            elseif contains(child_name, 'full_session')
                continue
            elseif contains(child_name, 'STOP')
                continue
            elseif contains(child_name, 'z_scores_whole_session.mat')
                continue
            elseif contains(child_name, 'drops')
                continue
            end
            load(horzcat(file, name, '\', child_name, '\dFF_data.mat'));
            if ~isempty(left_trials)
                US_trials_percent = left_trials;
                left_trials = [];
            elseif ~isempty(right_trials)
                US_trials_percent = right_trials;
                right_trials = [];
            end
            diffs = (US_trials_percent - avg);
            z_US_trials = (diffs ./ std_avg);
            fake_name = strrep(name, 'control', 'signal');
            cur_file = horzcat(file, fake_name, '\', child_name);
            if exist(horzcat(cur_file, '\artifacts_removed'), 'dir')
                load(horzcat(cur_file, '\', 'artifacts_removed\deleted_trials.mat'));
                if ~isempty(to_remove)
                if to_remove(end) <= size(z_US_trials, 1)
                z_US_trials(to_remove, :) = [];
                end
                end
                to_remove = [];
            end
            mean_z_US_trials = mean(z_US_trials);
            save(horzcat(file, name, '\', child_name, '\consistent_z_scores.mat'), 'z_US_trials', 'mean_z_US_trials');
        end
        load(horzcat(file, fake_name, '\Corrected_data_full_session.mat'));
        new_z_score = iso(:, 2);
        new_z_score = (new_z_score' - (avg/100))/(std_avg/100);
        iso(:, 2) = new_z_score';
        control_z_scores_whole_session = iso;
        save(horzcat(file, fake_name, '\control_z_scores_whole_session.mat'), 'control_z_scores_whole_session');
    end
    
end
