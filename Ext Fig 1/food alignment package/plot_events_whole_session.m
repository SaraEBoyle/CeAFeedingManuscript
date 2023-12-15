function [] = plot_events_whole_session()
    file = horzcat(pwd, '/');
    tem = dir(horzcat(file, '*.xlsx'));
    loaded_events = {}; %Names are row 1, times row 2
    colors = {'black', 'red', 'blue', 'yellow', 'cyan', 'green', 'magenta'};
    %{
    for event_ind = 1:size(tem, 1)
        event_data = tem(event_ind);
        name = event_data.name;
        event_time = xlsread(name);
        name = name(1:(end - 5));
        loaded_events{1, event_ind} = name;
        loaded_events{2, event_ind} = event_time;
    end
    %}
    % Read csv from Boris
    T = readtable('boris.csv');
    
    %Convert weird table to cell
    X = table2cell(T);
    
    % Search for the row that starts with "Time"
    for w = 1:size(X, 1)
        current = X{w, 1};
        if strcmp(current, 'Time')
            start_row = w;
            break
        end
    end
    %Remove junk from beginning
    X(1:(start_row - 1), :) = [];
    
    %Find unique Behavior events
    for w = 1:size(X, 2)
        current = X{1, w};
        if strcmp(current, 'Behavior')
            behavior_ind = w;
            break
        end
    end
    
    %Find whether the event is a start or stop or point
    for w = 1:size(X, 2)
        current = X{1, w};
        if strcmp(current, 'Status')
            status_ind = w;
            break
        end
    end
    
    behavior_column = X(:, behavior_ind);
    behavior_names = unique(behavior_column);
    
    %Remove "Behavior"
    behavior_names(1) = [];
    behavior_names = behavior_names';
    behavior_column(1) = [];
    
    %Find status
    status_column = X(:, status_ind);
    status_column(1) = [];
    
    %Combine behavior and status to get names
    event_names = {};
    for q = 1:size(status_column, 1)
        current_behavior = behavior_column{q};
        current_status = status_column{q};
        event_names{1, q} = horzcat(current_behavior, ' ', current_status);
    end
    
    behavior_names = unique(event_names);
    
    behavior_times = X(2:end,1);
    loaded_events = behavior_names;
    
    %for each event name, find times, create loaded_events cell
    for name_ind = 1:size(behavior_names, 2)
        name = behavior_names{name_ind};
        cur_times = [];
        for cur_name_ind = 1:size(event_names, 2)
            cur_name = event_names{cur_name_ind};
            if strcmp(cur_name, name)
                cur_time = behavior_times{cur_name_ind};
                cur_times = horzcat(cur_times, str2double(cur_time));
            end
        end
        loaded_events{2, name_ind} = cur_times;
    end
    
    %Make top row name, bottom row list of times
    
    
    signal_folders = dir(horzcat(file, '* signal'));
    for sig = 1:size(signal_folders, 1)
        current_folder = signal_folders(sig);
        name = current_folder.name;
        openfig(horzcat(file, name, '/Entire_Session_Plot_no_fit.fig'));
        ax = gca;
        y_limits = ax.YLim;
        
        %plot all the events and then rerun legend with gcamp iso plus the
        %names
        %legend([h1 h2 h3 h5],{'h1','h2','h3','h5})
        %set(hLeg,'visible','off')
        handles = {};
        for ev_ind = 1:size(loaded_events, 2)
            name_to_plot = loaded_events{1, ev_ind};
            times_to_plot = loaded_events{2, ev_ind};
            color_to_plot = colors{ev_ind};
            first = 1;
            for points = 1:size(times_to_plot, 2)
                %Remember the first handle of each group for the legend
                hand = plot([times_to_plot(points), times_to_plot(points)], y_limits, 'color', color_to_plot);
                if first
                    eval(horzcat('hand', num2str(ev_ind), ' = hand'));
                    handles{ev_ind} = eval(horzcat('hand', num2str(ev_ind)));
                end
                first = 0;
            end
        end
        %legend(handles(1,:), loaded_events(1, 1:end));%{'dark chocolate', 'eraser', 'water', 'white chocolate drop', 'white chocolate'});
        if size(handles, 2) == 5
            legend([hand1 hand2 hand3 hand4 hand5], loaded_events(1, 1:end));%{'dark chocolate', 'eraser', 'water', 'white chocolate drop', 'white chocolate'});
        else
            legend([hand1 hand2 hand3 hand4 hand5 hand6], loaded_events(1, 1:end));%{'dark chocolate', 'eraser', 'water', 'white chocolate drop', 'white chocolate'});
        end
        savefig(horzcat(file, name, '/Entire_Session_w_events.fig'))
        close all
    end
end