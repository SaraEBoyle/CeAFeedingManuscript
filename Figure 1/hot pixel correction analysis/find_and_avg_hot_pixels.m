close all

P = pwd;
S = dir(fullfile(P,'*.tif'));
filePath = matlab.desktop.editor.getActiveFilename;
filePath = filePath(1:(end-25));
pix = dir(fullfile(filePath,'*pixels.csv'));
N = {S.name};
pix_name = pix.name;
X = ~cellfun('isempty',strfind(N,'tif'));
%Coords = ~cellfun('isempty',strfind(P,'csv'));
%load(N{1});

%% Read list of dead pixels and replace them with average of surrounding good pixels
%name = 'D:\D:\GRIN imaging\mouse 1\09.26.21 all liquids\fat sucrose quinine 7 ul\files\fat sucrose quinine 7 u croppedl.tif';
name = N{1};
pixels = readtable([filePath, pix_name]);
pixels = pixels{:,1:2} + 1;
[bad_pixels,I,J] = unique(pixels, 'rows', 'first'); % find unique bad pixels

% for each bad pixel, find the surround pixels and check if they're also
% bad. Take an average of the surrounding good pixels and replace the bad
% one

pixels_for_correction = {}; %cell of bounding pixels to use for average correction. To use with bad_pixels array

for bad_ind = 1:size(bad_pixels, 1)
    bad_pixel = bad_pixels(bad_ind, :);
    edge = bad_pixel < 2;
    if sum(edge) > 0
        disp('Bad pixel on edge. Dont add 0 index into bounding pixels');
    end
    bad_x = bad_pixel(1);
    bad_y = bad_pixel(2);
    good_left = [bad_x - 1, bad_y];
    good_up = [bad_x, bad_y - 1];
    good_right = [bad_x + 1, bad_y];
    good_down = [bad_x, bad_y + 1];
    
    surrounding_pixels = [good_left;good_up;good_right;good_down];
    [trash,bad_surrounds,J] = unique([bad_pixels;surrounding_pixels], 'rows', 'first');
    bad_surrounds = size(bad_surrounds, 1) - size(bad_pixels, 1) - 4;
    if bad_surrounds > 0
        disp('bounding pixel bad. FIGURE IT OUT');
    end
    
    pixels_for_correction{bad_ind} = surrounding_pixels;
    
end


%% read data and convert to double

%V = [99,10,50,10,0,-1,50];
%[U,idx,idy] = unique(V,'first');
%cnt = histc(V,U);
%U(cnt>1) % duplicate values
%{
info = imfinfo(name);
numberOfPages = length(info);
full_file = zeros(472, 510, 21280);
halp = Tiff(name, 'r');
for k = 1:21280
    % Read the kth image in this multipage tiff file.
    %thisPage = imread(name, k);
    %image3D(:, :, 2);
    %why = imread(name, k);
    %imshow(why);
    full_file(:, :, k) = halp.read();
    % Now process thisPage somehow...
end	
halp.close();
%}
info = imfinfo(name);
bigness = info.FileSize;
if bigness < 4000000000
%% For files under 4 Gb
    Yf = read_file(name);
else
    disp('File larger than 4Gb detected. Converting to readable format...');
    numFramesStr = regexp(info.ImageDescription, 'images=(\d*)', 'tokens');
    numFrames = str2double(numFramesStr{1}{1});
    % Use low-level File I/O to read the file
    fp = fopen(name , 'rb');
    % The StripOffsets field provides the offset to the first strip. Based on
    % the INFO for this file, each image consists of 1 strip.
    fseek(fp, info.StripOffsets, 'bof');
    % Assume that the image is 16-bit per pixel and is stored in big-endian format.
    % Also assume that the images are stored one after the other.
    % For instance, read the first 100 frames
    imData = zeros(info.Height, info.Width, numFrames);
    framenum=numFrames;
    %imData=cell(1,framenum);
    for cnt = 1:framenum
        temp = fread(fp, [info.Width info.Height], 'uint16', 0, 'ieee-be')';
        imData(:, :, cnt) = temp;
        %imshow(imData(:, :, cnt));
        %imData{cnt} = fread(fp, [info.Width info.Height], 'uint16', 0, 'ieee-be')';
    end
    fclose(fp);

    %for frame =1:1200
    %    ip(:,:,frame) = imread(name, frame);
    %end

    Yf = single(imData); %only when normal method doesnt work
    %Yf = single(YF); %this is for normal files
    disp('File converted. Initiating pixel correction.')
end

%% For every bad pixel, average bounding pixels through length of video
for z_ind = 1:size(Yf, 3)
    % for each bad pixel % for some reason this tif is stored Y x X x Z
    %imshow(Yf(:, :, z_ind), [34000 36000]);
    for bad_ind = 1:size(bad_pixels, 1)
        bad_pixel = bad_pixels(bad_ind, :);
        good_pixels = pixels_for_correction{bad_ind};
        added_surrounds = [];
        % make for loop to go through each good pixel, y = 1, x = 2 :(
        for good_ind = 1:size(good_pixels, 1)
            good_value = Yf(good_pixels(good_ind, 2), good_pixels(good_ind, 1), z_ind);
            added_surrounds = [added_surrounds, good_value];
        end
        mean_surrounding = round(mean(added_surrounds));
        %y = 1, x = 2 :(
        Yf(bad_pixel(1, 2), bad_pixel(1,1), z_ind) = mean_surrounding;
    end
    %imshow(Yf(:, :, z_ind), [34000 36000]);
end


disp('Saving bad pixel corrected .tif...');
FileName_Save = [name(1:(end - 4)),'_removed_bad.tif'];
opts_tiff.append = true;
opts_tiff.big = true;
saveastiff(Yf,FileName_Save,opts_tiff);








