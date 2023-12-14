function [] = plot_individual_lines(time_axis, y_axis, params)
    %% Plot each cage, with whichever data is inputed. Not the main function.
    [exp, con, lowest_y, highest_y, dark_starts, dark_ends, time_from, time_to, ...
        screen, smooth_span, time_align_1st_dark_hr] =  params{[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]};
    figure
    h1 = plot(time_axis, y_axis(:, exp), 'g');
    hold on
    h2 = plot(time_axis, y_axis(:, con), 'k');
    y1 = ones(1, size(time_align_1st_dark_hr, 2))*highest_y;
    y2 = ones(1, size(time_align_1st_dark_hr, 2))*lowest_y;
    % Plot dark cycles
    for day = 1:size(dark_ends, 2)
        first_patch = zeros(1, size(time_align_1st_dark_hr, 2));
        first_patch(dark_starts(day):(dark_ends(day))) = 1;
        first_patch = first_patch == 1;
        P = patch([time_align_1st_dark_hr(first_patch) fliplr(time_align_1st_dark_hr(first_patch))], [y1(first_patch), fliplr(y2(first_patch))], [0.90196, 0.90196, 1.00000], 'FaceAlpha', 0.5, 'LineStyle', 'none');    
    end
    
    legend([h1(1) h2(2) P],'Exp', 'Con', 'Dark');
    xlim([time_from time_to]);
    ylim([lowest_y, highest_y]);
    set(gcf,'position',[100,100,screen(3) - 200,screen(4)/3]);
    
    
  
    
    %{
    figure
    h1 = plot(time_align_1st_dark_hr, smoothed_RERs(:, exp), 'g');
    hold on
    h2 = plot(time_align_1st_dark_hr, smoothed_RERs(:, con), 'k');
    y1 = ones(1, size(time_align_1st_dark_hr, 2))*highest_y;
    y2 = ones(1, size(time_align_1st_dark_hr, 2))*lowest_y;
    % Plot dark cycles
    for day = 1:size(dark_ends, 2)
        first_patch = zeros(1, size(time_align_1st_dark_hr, 2));
        first_patch(dark_starts(day):(dark_ends(day))) = 1;
        first_patch = first_patch == 1;
        P = patch([time_align_1st_dark_hr(first_patch) fliplr(time_align_1st_dark_hr(first_patch))], [y1(first_patch), fliplr(y2(first_patch))], [0.90196, 0.90196, 1.00000], 'FaceAlpha', 0.5, 'LineStyle', 'none');    
    end
    legend([h1(1) h2(2) P],'Exp', 'Con', 'Dark');
    xlim([time_from time_to]);
    ylim([lowest_y, highest_y]);
    set(gcf,'position',[100,100,screen(3) - 200,screen(4)/3]);
    title(['Smoothed RERs with window size: ', num2str(smooth_span), ', no down sampling']);
    savefig('smoothed_RER_individuals.fig');
    
    %% Plot downsampled individual lines
    figure;
    h1 = plot(ds_time_align_1st_dark_hr, smooth_ds_smooth_exp_RERs, 'g');
    hold on
    h2 = plot(ds_time_align_1st_dark_hr, smooth_ds_smooth_con_RERs, 'k');
    for day = 1:size(dark_ends, 2)
        first_patch = zeros(1, size(time_align_1st_dark_hr, 2));
        first_patch(dark_starts(day):dark_ends(day)) = 1;
        first_patch = first_patch == 1;
        P = patch([time_align_1st_dark_hr(first_patch) fliplr(time_align_1st_dark_hr(first_patch))], [y1(first_patch), fliplr(y2(first_patch))], [0.90196, 0.90196, 1.00000], 'FaceAlpha', 0.5, 'LineStyle', 'none');    
    end
    legend([h1(1) h2(2), P],'Exp', 'Con', 'Dark');
    xlim([time_from time_to]);
    set(gcf,'position',[100,150 + screen(4)/3,screen(3) - 200,screen(4)/3]);
    title(['Smoothed RERs with window size: ' num2str(smooth_span), ', with down sampling factor: ', num2str(rf)]);
    ylim([lowest_y, highest_y]);
    savefig('ds_smoothed_RER_individuals.fig');
    input('press enter to confirm');
    %}
end