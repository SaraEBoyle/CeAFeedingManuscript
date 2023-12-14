function [] = plot_pca_trajectories()
%% Plot the trajectories of the responses in 3D in 3 principal components
% 1 mo, 2 xg, 3 w, 4 fat, 5 sugar, 6 quinine
    
    % Plot mineral oil
    plot3(COEFF(1:61, 1)', COEFF(1:61, 2)', COEFF(1:61,3), 'Color', 'r', 'LineWidth', 2);
    hold on
    plot3(COEFF(1, 1)', COEFF(1, 2)', COEFF(1,3), 'Color', 'r', 'Marker', '*', 'MarkerSize', 12);
    
    % Plot Xanthan Gum
    plot3(COEFF(62:123, 1)', COEFF(62:123, 2)', COEFF(62:123,3), 'Color', 'b', 'LineWidth', 2);
    plot3(COEFF(62, 1)', COEFF(62, 2)', COEFF(62,3), 'Color', 'b', 'Marker', '*', 'MarkerSize', 12);
    
    % Plot water
    plot3(COEFF(123:184, 1)', COEFF(123:184, 2)', COEFF(123:184,3), 'Color', 'c', 'LineWidth', 2);
    plot3(COEFF(123, 1)', COEFF(123, 2)', COEFF(123,3), 'Color', 'c', 'Marker', '*', 'MarkerSize', 12);
    
    % plot fat
    plot3(COEFF(184:245, 1)', COEFF(184:245, 2)', COEFF(184:245,3), 'Color', 'y', 'LineWidth', 2);
    plot3(COEFF(184, 1)', COEFF(184, 2)', COEFF(184,3), 'Color', 'y', 'Marker', '*', 'MarkerSize', 12);
    
    % plot sugar
    plot3(COEFF(245:306, 1)', COEFF(245:306, 2)', COEFF(245:306,3), 'Color', 'm', 'LineWidth', 2);
    plot3(COEFF(245, 1)', COEFF(245, 2)', COEFF(245,3), 'Color', 'm', 'Marker', '*', 'MarkerSize', 12);
    
    % plot quinine
    plot3(COEFF(306:367, 1)', COEFF(306:367, 2)', COEFF(306:367,3), 'Color', 'g', 'LineWidth', 2);
    plot3(COEFF(306, 1)', COEFF(306, 2)', COEFF(306,3), 'Color', 'g', 'Marker', '*', 'MarkerSize', 12);
    
    % plot air
    plot3(COEFF(368:387, 1)', COEFF(368:387, 2)', COEFF(368:387,3), 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
    plot3(COEFF(368, 1)', COEFF(368, 2)', COEFF(368,3), 'Color', [0.6350 0.0780 0.1840], 'Marker', '*', 'MarkerSize', 12);
    
    %plot(COEFF(62:123, 1)', COEFF(62:123, 2)', 'r');
    %plot(COEFF(123:184, 1)', COEFF(123:184, 2)', 'm');
    %plot(COEFF(184:245, 1)', COEFF(184:245, 2)', 'g');
    %plot(COEFF(245:306, 1)', COEFF(245:306, 2)', 'c');
    %plot(COEFF(306:367, 1)', COEFF(306:367, 2)', 'y');
    %plot(COEFF(368:387, 1)', COEFF(368:387, 2)', 'b');
    
    scatter3(COEFF(1:61,1),COEFF(1:61,2),COEFF(1:61,3));
    hold on
    scatter3(COEFF(62:123,1),COEFF(62:123,2),COEFF(62:123,3), 'r');
    scatter3(COEFF(123:184,1),COEFF(123:184,2),COEFF(123:184,3), 'm');
    scatter3(COEFF(184:245,1),COEFF(184:245,2),COEFF(184:245,3), 'g');
    scatter3(COEFF(245:306,1),COEFF(245:306,2),COEFF(245:306,3), 'c');
    scatter3(COEFF(306:367,1),COEFF(306:367,2),COEFF(306:367,3), 'y');
    scatter3(COEFF(368:387,1),COEFF(368:387,2),COEFF(368:387,3), 'b');
    axis equal
    set(gca,'YTickLabel',[]);
    set(gca,'ZTickLabel',[]);
    set(gca,'xTickLabel',[]);
    for AZ = 0:1:360
        view(AZ, 5);
        F(AZ + 1) = getframe(gcf); 
        pause(0.1);
    end
    writerObj = VideoWriter('Principal Components.mp4', 'MPEG-4');
    writerObj.FrameRate = 10;
    open(writerObj);
    writeVideo(writerObj, F)
    close(writerObj);
    
    xlabel('1st Principal Component')
    ylabel('2nd Principal Component')
    zlabel('3rd Principal Component')
    
    
    
    
    %Train a discriminant anlysis classifier
    %1 row per observation, 1 column per predictor in response matrix
    tree = fitcdiscr(responses, labels);
    
    %Use it to predict which liquid you give
    meanclass = predict(MdlLinear,meanmeas);
    
    
end