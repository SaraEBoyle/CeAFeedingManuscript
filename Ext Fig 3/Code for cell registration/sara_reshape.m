%% Reshapes sorted into a shape that can be used for cell registration
clear
files = {pwd};

%{%'D:\GRIN imaging\mouse 1\11.13.21 fr high fat sucrose palm fed3\', ...
    %'D:\GRIN imaging\mouse 1\11.12.21 sated high fat palm sugar fed3\', ...
    %'E:\Pkcd GRIN imaging\GRIN 1m\03.16.22 6 hr fr hf su fed3\'
    %};
    %'/Users/sboyle/Documents/Li Lab/Image analysis/SST16/20181114/saline/', ...
    %'/Users/sboyle/Documents/Li Lab/Image analysis/SST16/20181114/amph/'
adapt_sorted_before_cellreg;    
    
i = 0;
for x = files
    disp(x);
    load(fullfile(cell2mat(x), '\Sorted.mat'));
    %% reshape ROIs NeuroNum x Yaxis x Xaxis %%
    i = i + 1;
    Yaxis = ROIsize(1);
    Xaxis = ROIsize(2);
    NeuroNum = size(ROIs,2);
    ROIs_reshape = reshape(ROIs,[Yaxis,Xaxis,NeuroNum]);
    allFiltersMat = permute(ROIs_reshape,[3,1,2]);
    o = strcat(string(x), '\spatial_footprints'); %, num2str(i));
    disp(o);
    save(o, 'allFiltersMat');
    
end
