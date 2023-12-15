function [] = calculate_auc_of_bouts()
    %Find the bout length
    load('loaded_events.mat');
    tem = dir('*signal');
    for folder_ind = 1:size(tem, 1)
        folder = tem(folder_ind);
        load(horzcat(folder.folder, '/', folder.name, '/', 'Corrected_data_full_session.mat'));
        
    end
end