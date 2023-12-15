function[data] = load_data(signal_key_word, dir_nm, photometry_data)
    %% select EthoVision data from excel
    % clear nam
    %% Load in data
    %Retrieve the bpod behavior file
    file = dir_nm;
    if photometry_data
    %Load photometry data
        signal_file = dir(horzcat(file, '*', signal_key_word, '*.csv'));
        data = dlmread(horzcat(file, signal_file.name));
    end
    
    if photometry_data
        %Convert analog in file to seconds
        data(:,1) = (data(:,1) - data(1,1))/1000;
    end
end