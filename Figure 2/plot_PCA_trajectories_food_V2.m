function [] = plot_PCA_trajectories_food_V2()
%% Loads z scores and food contact times, deletes sketchy neurons, then runs PCA on 10 seconds post contact
    load('z_scores_session_1.mat');
    NeuronZScores = NeuronZScores_aligned; % use info aligned to food pellet
    NeuronZScores(neurons_to_delete, :) = [];
    load('time_of_food_contact.mat');
    close all
    response_length = 10;
    response = [target_frame - 35 target_frame + response_length*10]; % time post-delivery to measure
    baseline = [2 30]; % baseline to include
    x_axis = 0.1:.1:trial_length;
    record_avg_time_trajectories = 0;
    plot_trial_by_trial_trajectories = 0;
    plot_trial_by_trial_points = 1;
    
    % holds appended averaged high fat + sucrose + missed responses X n neurons
    % keeps timing intact at expense of trial by trial resolution
    trajectory_response_matrix = [];
    
    % trial by trial trajectory responses. PCA over entire dataset, no
    % averaging
    trajectory_response_matrix_all_trials = [];
    
    % average trial by trial responses, averaged over 10 seconds
    avg_response_matrix_trial_by_trial = [];
    
    
    % This particular session requires trunctation
    %trial_inds = trial_inds(1:50);
    %relative_sorted_food_times = relative_sorted_food_times(1:50);
    
    %Need to get the food specific index of the trial
    first_food_trials = find(trial_inds == 1);
    second_food_trials = find(trial_inds == 2);
    third_food_trials = find(trial_inds == 3);
    response_matrix_trial_by_trial = [];
    first_response_all = [];
    second_response_all = [];
    third_response_all = [];
    %for each neuron do this. Redundant but no one pays me for speed
    for neuron_ind = 1:size(NeuronZScores, 1)
        % For each food contact time, get together z scores response_length secs after
        
        %this is to keep track of trial inds
        trial_ind_mat = [1, 1, 1];
        
        % keep track of single trial responses
        first_response_mat = [];
        first_baseline_mat = [];
        second_response_mat = [];
        second_baseline_mat = [];
        third_response_mat = [];
        third_baseline_mat = [];
        
        % smush together single trial responses for ultimate PCA
        first_response_mat_all = [];
        second_response_mat_all = [];
        missed_response_mat_all = [];
        
        current_neuron = NeuronZScores(neuron_ind, :);
        first_trials = current_neuron{:, 1};
        second_trials = current_neuron{:, 2};
        third_trials = current_neuron{:, 3};
        
        first_response_mat = [first_response_mat; first_trials(:, response(1):response(2))];
        first_baseline_mat = [first_baseline_mat, first_trials(:, baseline(1):baseline(2))];
        
        second_response_mat = [second_response_mat; second_trials(:, response(1):response(2))];
        second_baseline_mat = [second_baseline_mat, second_trials(:, baseline(1):baseline(2))];
        
        third_response_mat = [third_response_mat; third_trials(:, response(1):response(2))];
        third_baseline_mat = [third_baseline_mat, third_trials(:, baseline(1):baseline(2))];
        
        % Average the trial responses by neuron
        trajectory_response_matrix(neuron_ind, :) = horzcat(mean(first_response_mat), mean(second_response_mat), mean(third_response_mat));
        
        %neuron x avg response
        avg_response_matrix_trial_by_trial(neuron_ind, :) = [mean(first_response_mat, 2)', mean(second_response_mat, 2)', mean(third_response_mat, 2)'];
        %response_matrix_trial_by_trial = [response_matrix_trial_by_trial;first_response_mat;second_response_mat;third_response_mat];
        % (# fat trials per individual neuron 19*52) x frames in response
        first_response_all = [first_response_all;reshape(first_response_mat, 1, size(first_response_mat, 1) * size(first_response_mat, 2))];
        second_response_all = [second_response_all;reshape(second_response_mat, 1, size(second_response_mat, 1) * size(second_response_mat, 2))];
        third_response_all = [third_response_all;reshape(third_response_mat, 1, size(third_response_mat, 1) * size(third_response_mat, 2))];
    end
    response_matrix_trial_by_trial = [first_response_all,second_response_all];%,third_response_all];
    trial_ind_mat = [size(NeuronZScores{1, 1}, 1), size(NeuronZScores{1, 2}, 1), size(NeuronZScores{1, 3}, 1)];
    %% Make video of PCA component traces
    % Now run PCA on the response matrix. This matrix is average trial
    % response, averaged over 10 seconds. This keeps timing in tact but
    % sacrifices trial by trial resolution
    % Also saves video of the PCA traces
    if 1 %GOOD
    %% PCA on averaged data
        [COEFF,score,latent] = pca(trajectory_response_matrix'); %1/21 changed to ' and switched from plotting coeff to score 
        bigness = size(first_response_mat, 2);
        %score, which is
        %the representation of X in the principal component space. Rows of SCORE
        %correspond to observations, columns to components.
    
        % Plot mineral oil
        plot3(score(1:bigness, 1)', score(1:bigness, 2)', score(1:bigness,3), 'Color', 'r', 'LineWidth', 2);
        hold on
        plot3(score(1, 1)', score(1, 2)', score(1,3), 'Color', 'r', 'Marker', '*', 'MarkerSize', 12);

        % Plot Xanthan Gum
        plot3(score((bigness + 1):(bigness*2), 1)', score((bigness + 1):(bigness*2), 2)', score((bigness + 1):(bigness*2),3), 'Color', 'b', 'LineWidth', 2);
        plot3(score(bigness + 1, 1)', score(bigness + 1, 2)', score(bigness + 1,3), 'Color', 'b', 'Marker', '*', 'MarkerSize', 12);

        % Plot water
        plot3(score((bigness*2 + 1):(bigness*3), 1)', score((bigness*2 + 1):(bigness*3), 2)', score((bigness*2 + 1):(bigness*3),3), 'Color', 'c', 'LineWidth', 2);
        plot3(score((bigness*3), 1)', score((bigness*3), 2)', score((bigness*3),3), 'Color', 'c', 'Marker', '*', 'MarkerSize', 12);

    
        grid on
        axis equal
        set(gca,'YTickLabel',[]);
        set(gca,'ZTickLabel',[]);
        set(gca,'xTickLabel',[]);
        for AZ = 0:1:360
            view(AZ, 5);
            F(AZ + 1) = getframe(gcf); 
            pause(0.1);
        end
        writerObj = VideoWriter('Principal Components with baseline.mp4', 'MPEG-4');
        writerObj.FrameRate = 10;
        open(writerObj);
        writeVideo(writerObj, F)
        close(writerObj);
    end
    
    %% Plots average single neuron trajectory in PCA 3D space for each trial (not just the average over 10 seconds)
    if 0 % has some major issues
        %% PCA on single trials
        [COEFF_all,score_all,~,~,explainedVar_all] = pca(trajectory_response_matrix_all_trials');%1/21 changed to ' and switched from plotting coeff to score
        COEFF_all = score_all;
        % Plot first pellet trials
        frame_num = response_length*10;

        for first_ind = 1:(trial_ind_mat(1) - 1)
            frame_start = ((first_ind - 1)*100 + 1);
            frame_end = (first_ind*100 + 1);
            plot3(COEFF_all(frame_start:frame_end, 1)', COEFF_all(frame_start:frame_end, 2)', COEFF_all(frame_start:frame_end,3), 'Color', 'r', 'LineWidth', 2);
            if first_ind == 1
                hold on
            end
            plot3(COEFF_all(frame_start, 1)', COEFF_all(frame_start, 2)', COEFF_all(frame_start,3), 'Color', 'r', 'Marker', '*', 'MarkerSize', 12);
        end
        % Plot second pellet trials
        for second_ind = 1:(trial_ind_mat(2) - 1)
            previous_pellet_frame = size(first_response_mat_all, 2);
            frame_start = previous_pellet_frame + ((second_ind - 1)*100 + 1);
            frame_end = previous_pellet_frame + (second_ind*100 + 1);
            plot3(COEFF_all(frame_start:frame_end, 1)', COEFF_all(frame_start:frame_end, 2)', COEFF_all(frame_start:frame_end,3), 'Color', 'b', 'LineWidth', 2);
            plot3(COEFF_all(frame_start, 1)', COEFF_all(frame_start, 2)', COEFF_all(frame_start,3), 'Color', 'b', 'Marker', '*', 'MarkerSize', 12);
        end

        % Plot missed trials
        for third_ind = 1:(trial_ind_mat(3) - 1)
            previous_pellet_frame = size(second_response_mat_all, 2) + size(first_response_mat_all, 2);
            frame_start = previous_pellet_frame + ((third_ind - 1)*100 + 1);
            frame_end = previous_pellet_frame + (third_ind*100 + 1);
            plot3(COEFF_all(frame_start:frame_end, 1)', COEFF_all(frame_start:frame_end, 2)', COEFF_all(frame_start:frame_end,3), 'Color', 'c', 'LineWidth', 2);
            plot3(COEFF_all(frame_start, 1)', COEFF_all(frame_start, 2)', COEFF_all(frame_start,3), 'Color', 'c', 'Marker', '*', 'MarkerSize', 12);
        end
    end
    
    %% Plot PCA of average trial by trial responses. Then train an svm classifier
    % First trains by reserving all combos of 1 and 2 trials for test
    % datasets, keeping track of how often predictions are correct. Then
    % does the same with shuffled labels. Prints % of correct predictions
    close all
    neuron_num = size(NeuronZScores,1);
    if 1
   %     [COEFF,score,latent] = pca(avg_response_matrix_trial_by_trial'); %1/21 changed to ' and switched from plotting coeff to score, this is just avg responses
        [COEFF,score,latent] = pca(response_matrix_trial_by_trial'); %this is mean trace responses
        frame_num = 1;
        COEFF_all = score; %avg trial response x reduced dim
        response_length = response(2) - response(1) + 1;
        fat_coords = [];
        sucrose_coords = [];
        for first_ind = 1:trial_ind_mat(1)
            %frame_start = ((first_ind - 1)*1 + 1);
            
            % 19 fat trials repeated for each neuron
            frame_start = (1 + response_length*(first_ind - 1)):(response_length*first_ind);
            %plot3(COEFF_all(frame_start, 1)', COEFF_all(frame_start, 2)', COEFF_all(frame_start,3), 'Color', 'r', 'LineWidth', 2);
            %plot3(mean(COEFF_all(frame_start, 1))', mean(COEFF_all(frame_start, 2))', mean(COEFF_all(frame_start,3))', 'Color', 'r', 'LineWidth', 2);
            if first_ind == 1
                hold on
            end
            fat_coords = [fat_coords;[mean(COEFF_all(frame_start, 1))', mean(COEFF_all(frame_start, 2))', mean(COEFF_all(frame_start,3))']];
            plot3(mean(COEFF_all(frame_start, 1))', mean(COEFF_all(frame_start, 2))', mean(COEFF_all(frame_start,3))', 'Color', 'r', 'Marker', '*', 'MarkerSize', 12);
            %plot3(mean(COEFF_all(frame_start, 1))', mean(COEFF_all(frame_start, 2))', mean(COEFF_all(frame_start,3))', 'Color', 'b', 'Marker', '*', 'MarkerSize', 12);
            % plot3(COEFF_all(frame_start, 1)', COEFF_all(frame_start, 2)', COEFF_all(frame_start,3), 'Color', 'r', 'Marker', '*', 'MarkerSize', 12);
        end
        % Plot second pellet trials
        for second_ind = 1:(trial_ind_mat(2))
            response_length = response(2) - response(1) + 1;
            % 19 fat trials repeated for each neuron
            frame_start = ((1 + response_length*(second_ind - 1)):(response_length*second_ind)) + first_ind*response_length;
            %plot3(COEFF_all(frame_start, 1)', COEFF_all(frame_start, 2)', COEFF_all(frame_start,3), 'Color', 'r', 'LineWidth', 2);
            %plot3(mean(COEFF_all(frame_start, 1))', mean(COEFF_all(frame_start, 2))', mean(COEFF_all(frame_start,3))', 'Color', 'r', 'LineWidth', 2);
            if second_ind == 1
                hold on
            end
            sucrose_coords = [sucrose_coords;[mean(COEFF_all(frame_start, 1))', mean(COEFF_all(frame_start, 2))', mean(COEFF_all(frame_start,3))']];
            plot3(mean(COEFF_all(frame_start, 1))', mean(COEFF_all(frame_start, 2))', mean(COEFF_all(frame_start,3))', 'Color', 'b', 'Marker', '*', 'MarkerSize', 12);
        end
        % You want to trim off missed trials
        %X = [score(1:(end - (trial_ind_mat(3)*response_length - 1)), [1,2,3])];
        X = [fat_coords;sucrose_coords];
        first_ones = ones(1, trial_ind_mat(1)); %label first type 1
        second_zeros = zeros(1, trial_ind_mat(2)); %label second type 0
        group = [first_ones, second_zeros];
        %Plot the graph and the plane
        %fit simple classifier
        mdl = fitcsvm(X,group);
        svm_3d_plot(mdl,X,group)
        
        %go through whole dataset and train on all but 1, then test on the
        %1 left out.
        corrects = 0;
        one_left_out = 1:size(X, 1);
        two_left_out = nchoosek(one_left_out, 2);
        one_left_out = [one_left_out', zeros(size(one_left_out, 2), 1)];
        one_and_two_left = vertcat(one_left_out, two_left_out);
        for left_out_ind = 1:size(one_and_two_left, 1)
            % get all permutations of groups size 1 and 2 to use for
            % testing
            left_out = one_and_two_left(left_out_ind, :);
            left_out = nonzeros(left_out);
            test_X = X(left_out, :);
            test_label = group(left_out);
            train_X = X;
            train_X(left_out, :) = [];
            train_labels = group;
            train_labels(left_out) = [];
            train_mdl = fitcsvm(train_X,train_labels);
            prediction = predict(train_mdl,test_X);
            % if first test was right
            if prediction(1) == test_label(1)
                corrects = corrects + 1;
            else
                %disp('got one wrong')
            end
            %if the test is 2 items, not 1
            if size(test_label, 2) == 2
                if prediction(2) == test_label(2)
                    corrects = corrects + 1;
                else
                    %disp('got one wrong')
                end
            end
        end
        percent_correct = corrects/size(nonzeros(one_and_two_left), 1);
        disp(horzcat('Percent correct: ', num2str(percent_correct)))
        
        %mdl = fitclinear(X,[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
        %    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]);
        
        %% Now do the same but shuffle your labels
        corrects = 0;
        for left_out_ind = 1:size(one_and_two_left, 1)
            % get all permutations of groups size 1 and 2 to use for
            % testing
            left_out = one_and_two_left(left_out_ind, :);
            left_out = nonzeros(left_out);
            
            rand_inds = randperm(size(X, 1));
            new_group = group(rand_inds);
            
            test_X = X(left_out, :);
            test_label = new_group(left_out);
            train_X = X;
            train_X(left_out, :) = [];
            train_labels = new_group;
            train_labels(left_out) = [];
            train_mdl = fitcsvm(train_X,train_labels);
            prediction = predict(train_mdl,test_X);
            
            
            if prediction(1) == test_label(1)
                corrects = corrects + 1;
            else
                %disp('got one wrong')
            end
            %if the test is 2 items, not 1
            if size(test_label, 2) == 2
                if prediction(2) == test_label(2)
                    corrects = corrects + 1;
                else
                    %disp('got one wrong')
                end
            end
            
            
            %if prediction == test_label
            %    corrects = corrects + 1;
            %end
        end
        
        
        
        percent_correct_shuffled = corrects/size(nonzeros(one_and_two_left), 1);
        disp(horzcat('Percent shuffled correct: ', num2str(percent_correct_shuffled)))
        save('SVM_classifier', 'percent_correct_shuffled', 'percent_correct');
        
        %{
        x = linspace(0,5);
        y = linspace(0,5);
        z = linspace(0,5);
        [XX,YY,ZZ] = meshgrid(x,y,z);
        pred = [XX(:),YY(:),ZZ(:)];
        p = predict(mdl,pred);
        f = @(x) -(x*mdl.Beta(1) + mdl.Bias)/mdl.Beta(2);
        y = f(x);
        % Plot line dividing
        plot(x,y,'g--','LineWidth',2,'DisplayName','Boundary')
        %scatter3(pred(:,1),pred(:,2),pred(:,3),p)
        %}
        %plot(mdl)
        % Plot missed trials
        %{
        for third_ind = 1:(trial_ind_mat(3) - 1)
            previous_pellet_frame = trial_ind_mat(1) - 1 + trial_ind_mat(2) - 1;
            frame_start = previous_pellet_frame + ((third_ind - 1)*1 + 1);
            plot3(COEFF_all(frame_start, 1)', COEFF_all(frame_start, 2)', COEFF_all(frame_start,3), 'Color', 'c', 'LineWidth', 2);
            plot3(COEFF_all(frame_start, 1)', COEFF_all(frame_start, 2)', COEFF_all(frame_start,3), 'Color', 'c', 'Marker', '*', 'MarkerSize', 12);
        end
        %}
    end
    
   close all 
end