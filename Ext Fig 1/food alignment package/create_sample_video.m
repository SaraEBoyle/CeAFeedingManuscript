%% Make example video of dF/F
close all
times = z_scores_whole_session(:, 1);
signal = z_scores_whole_session(:, 2);
control = control_z_scores_whole_session(:, 2);

start_time = 190;
end_time = 270;

bigs = find(times >= start_time);
start_frame = bigs(1);

smalls = find(times <= end_time);
end_frame = smalls(end);

plot(times(start_frame:end_frame), signal((start_frame:end_frame)));
base = gca;
base_x = base.XLim;
base_y = base.YLim;
close all
baseplot = figure;
xlim(base_x);
ylim(base_y);
hold on
title('White Chocolate Trial', 'FontSize', 16);
ylabel('Z-Score(dF/F)', 'FontSize', 12);
for i = 1:(end_frame - start_frame)
      
    %imshow(processo(:,:,1,i))
    cur_frame = (start_frame - 1) + i;
    time = times(cur_frame);
    sig_point = signal(cur_frame);
    con_point = control(cur_frame);
    frame_to_plot = start_frame:cur_frame;
    frame_to_plot = frame_to_plot';
    times_to_plot = times(frame_to_plot);
    sigs_to_plot = signal(frame_to_plot);
    cont_to_plot = control(frame_to_plot);
    
    plot(times_to_plot', sigs_to_plot', 'g');
    plot(times_to_plot',cont_to_plot', 'm');
    %plot(X1,Y1,'o')
    %plot(X2,Y2,'o')
    %plot(X3,Y3,'o')
    %hold off
    F(i) = getframe(gcf);
    drawnow
end
% create the video writer with 1 fps
writerObj = VideoWriter('sampleVideo.avi');
writerObj.FrameRate = 10;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);