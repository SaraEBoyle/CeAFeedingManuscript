%%
%% TODO DELETE NEURONS THAT NEED DELETING 
clear,clc;close all
%% load the data and combine them
tem = dir('z_scores*.mat');
auROC_matrix = [];
aaMatrix = [];
NeuronSum1 = [];
NeuronSum2 = [];
NeuronSum3 = [];
NeuronSum4 = [];
NeuronSum5 = [];
NeuronSum6 = [];
NeuronSum7 = [];
use_dFF = 1;
PCA_low_bound = 4;
PCA_high_bound = 10;
cluster_num = 12; %choose # clusters
do_kmeans = 1;
time_interest_win = [4,10];
time = (1:140)/10;
time_sel = (time>=0 & time<=14);
%first use Wilcoxon test to find responses different from baseline
%Gives auROC_matrix [liquid 1 dif from baseline, p, liquid 2 dif, p]
all_neurons_to_delete = [];
for NN = 1:numel(tem)
    cur_name = tem(NN).name;
    load(tem(NN).name);
    NeuronZScores = NeurondFF;
    NeurondFF_mean = {};
    for x = 1:size(NeurondFF, 1)
        for y = 1:size(NeurondFF, 2)
            if isempty(NeurondFF{x, y})
                continue
            end
            NeurondFF_mean{x, y} = mean(NeurondFF{x, y});
        end
    end
    NeuronZScores_mean = NeurondFF_mean;
    
    all_neurons_to_delete = horzcat(all_neurons_to_delete, neurons_to_delete);
    if isempty(NeuronZScores{1, 1})
        pallet_cleansers = 1;
    else
        pallet_cleansers = 0;
    end
    if contains(cur_name, '1')
        for k=1:size(NeuronZScores,1)
            data1 = NeuronZScores{k,1 + pallet_cleansers};
            data2 = NeuronZScores{k,2 + pallet_cleansers};
            data3 = NeuronZScores{k,3 + pallet_cleansers};
            
            d1 = data1(:,time>=time_interest_win(1) & time<=time_interest_win(2));
            d2 = data2(:,time>=time_interest_win(1) & time<=time_interest_win(2));        
            d3 = data3(:,time>=time_interest_win(1) & time<=time_interest_win(2));
            
            d1_base = data1(:,time<time_interest_win(1));
            d2_base = data2(:,time<time_interest_win(1));
            d3_base = data3(:,time<time_interest_win(1));
            
            d1 = mean(d1,2);
            d2 = mean(d2,2);
            d3 = mean(d3,2);
            
            d1_base = mean(d1_base,2);
            d2_base = mean(d2_base,2);
            d3_base = mean(d3_base,2);
            
            p1 = signrank(d1, d1_base); % Wilcoxon signed rank test for zero median 
            %paired, two-sided test of the hypothesis
            %that the difference between the matched samples in the vectors X and Y
            %comes from a distribution whose median is zero
            au1 = mean(d1)-mean(d1_base); %dif between mean and baseline    
            p2 = signrank(d2, d2_base);
            au2 = mean(d2)-mean(d2_base);
            p3 = signrank(d3, d3_base);
            au3 = mean(d3)-mean(d3_base);
            
            auROC_matrix = [auROC_matrix;[au1, p1,au2, p2, au3, p3]];
            
            a1 = NeuronZScores_mean{k,1 + pallet_cleansers};
            a2 = NeuronZScores_mean{k,2 + pallet_cleansers};
            a3 = NeuronZScores_mean{k,3 + pallet_cleansers};
            NeuronSum1 = [NeuronSum1;a1(time_sel)];
            NeuronSum2 = [NeuronSum2;a2(time_sel)];
            NeuronSum3 = [NeuronSum3;a3(time_sel)];
        end
    elseif contains(cur_name, '2')
        for k=1:size(NeuronZScores,1)
            data4 = NeuronZScores{k,1 + pallet_cleansers};
            data5 = NeuronZScores{k,2 + pallet_cleansers};
            data6 = NeuronZScores{k,3 + pallet_cleansers};
            
            d4 = data4(:,time>=time_interest_win(1) & time<=time_interest_win(2));
            d5 = data5(:,time>=time_interest_win(1) & time<=time_interest_win(2));        
            d6 = data6(:,time>=time_interest_win(1) & time<=time_interest_win(2));
            
            d4_base = data4(:,time<time_interest_win(1));
            d5_base = data5(:,time<time_interest_win(1));
            d6_base = data6(:,time<time_interest_win(1));
            
            d4 = mean(d4,2);
            d5 = mean(d5,2);
            d6 = mean(d6,2);
            
            d4_base = mean(d4_base,2);
            d5_base = mean(d5_base,2);
            d6_base = mean(d6_base,2);
            
            p4 = signrank(d4, d4_base); % Wilcoxon signed rank test for zero median 
            %paired, two-sided test of the hypothesis
            %that the difference between the matched samples in the vectors X and Y
            %comes from a distribution whose median is zero
            au4 = mean(d4)-mean(d4_base); %dif between mean and baseline    
            p5 = signrank(d5, d5_base);
            au5 = mean(d5)-mean(d5_base);
            p6 = signrank(d6, d6_base);
            au6 = mean(d6)-mean(d6_base);
            
            auROC_matrix(k, 7:12) = [au4, p4,au5, p5, au6, p6];
            %auROC_matrix = [auROC_matrix;[au4, p4,au5, p5, au6, p6]];
            
            a4 = NeuronZScores_mean{k,1 + pallet_cleansers};
            a5 = NeuronZScores_mean{k,2 + pallet_cleansers};
            a6 = NeuronZScores_mean{k,3 + pallet_cleansers};
            NeuronSum4 = [NeuronSum4;a4(time_sel)];
            NeuronSum5 = [NeuronSum5;a5(time_sel)];
            NeuronSum6 = [NeuronSum6;a6(time_sel)];
        end
    elseif contains(cur_name, '3')
        for k=1:size(NeuronZScores,1)
            data7 = NeuronZScores{k,1 + pallet_cleansers};
            
            d7 = data7(:,time>=time_interest_win(1) & time<=time_interest_win(2));
           
            d7_base = data7(:,time<time_interest_win(1));
           
            d7 = mean(d7,2);

            d7_base = mean(d7_base,2);
            
            p7 = signrank(d7, d7_base); % Wilcoxon signed rank test for zero median 
            %paired, two-sided test of the hypothesis
            %that the difference between the matched samples in the vectors X and Y
            %comes from a distribution whose median is zero
            au7 = mean(d7)-mean(d7_base); %dif between mean and baseline    
            
            auROC_matrix(k, 13:14) = [au7, p7];
            %auROC_matrix = [auROC_matrix;[au4, p4,au5, p5, au6, p6]];
            
            a7 = NeuronZScores_mean{k,1 + pallet_cleansers};

            NeuronSum7 = [NeuronSum7;a7(time_sel)];
        end
    end
end
all_neurons_to_delete = unique(all_neurons_to_delete);

%% Population activity
cc = [0.2 0.8 0 %green
    0.8 0.2 0 %red
    .6 .8 1 %light blue
    1.0000 0.5000 0.8000 %pink
    0.4900 0.2800 0.5600 %purple
    0.9300 0.6900 0.1300 %orange
    0.8000 0.8000 0.8000 %gray
    0.5 0.5 0.5];
% no idea what cc is
figure('Position', [100 100 400 350],'Name','Summary','numbertitle','off');
hold on
NeuronSum1(all_neurons_to_delete, :) = [];
NeuronSum1 = NeuronSum1 - mean(mean(NeuronSum1(:, 5:39)));
a = NeuronSum1;
Y1 = mean(a,1); %mean all responses
Y1 = Y1 - mean(Y1(5:39));
Y1_err = std(a)/sqrt(size(a,1)); %std error
AreaPlot(time,smooth(Y1,3)',smooth(Y1_err)',cc(1,:),0.4,1);

NeuronSum2(all_neurons_to_delete, :) = [];
NeuronSum2 = NeuronSum2 - mean(mean(NeuronSum2(:, 5:39)));
a = NeuronSum2;
Y2 = mean(a,1);
Y2 = Y2 - mean(Y2(5:39));
Y2_err = std(a)/sqrt(size(a,1));
AreaPlot(time,smooth(Y2,3)',smooth(Y2_err)',cc(2,:),0.4,1);

NeuronSum3(all_neurons_to_delete, :) = [];
NeuronSum3 = NeuronSum3 - mean(mean(NeuronSum3(:, 5:39)));
a = NeuronSum3;
Y3 = mean(a,1);
Y3 = Y3 - mean(Y3(5:39));
Y3_err = std(a)/sqrt(size(a,1));
AreaPlot(time,smooth(Y3,3)',smooth(Y3_err)',cc(3,:),0.4,1);

NeuronSum4(all_neurons_to_delete, :) = [];
NeuronSum4 = NeuronSum4 - mean(mean(NeuronSum4(:, 5:39)));
a = NeuronSum4;
Y4 = mean(a,1);
Y4 = Y4 - mean(Y4(5:39));
Y4_err = std(a)/sqrt(size(a,1));
AreaPlot(time,smooth(Y4,3)',smooth(Y4_err)',cc(4,:),0.4,1);

NeuronSum5(all_neurons_to_delete, :) = [];
NeuronSum5 = NeuronSum5 - mean(mean(NeuronSum5(:, 5:39)));
a = NeuronSum5;
Y5 = mean(a,1);
Y5 = Y5 - mean(Y5(5:39));
Y5_err = std(a)/sqrt(size(a,1));
AreaPlot(time,smooth(Y5,3)',smooth(Y5_err)',cc(5,:),0.4,1);

NeuronSum6(all_neurons_to_delete, :) = [];
NeuronSum6 = NeuronSum6 - mean(mean(NeuronSum6(:, 5:39)));
a = NeuronSum6;
Y6 = mean(a,1);
Y6 = Y6 - mean(Y6(5:39));
Y6_err = std(a)/sqrt(size(a,1));
AreaPlot(time,smooth(Y6,3)',smooth(Y6_err)',cc(6,:),0.4,1);

NeuronSum7(all_neurons_to_delete, :) = [];
NeuronSum7 = NeuronSum7 - mean(mean(NeuronSum7(:, 5:39)));
a = NeuronSum7;
Y7 = mean(a,1);
Y7 = Y7 - mean(Y7(5:39));
Y7_err = std(a)/sqrt(size(a,1));
AreaPlot(time,smooth(Y7,3)',smooth(Y7_err)',cc(7,:),0.4,1);

xlabel('Time (s)','FontSize', 12)
ylabel('z-Score(dF/F)','FontSize', 12)
set(gca,'TickDir', 'out','xlim',[0,14],'xtick',2:2:10,'FontSize', 12,'box','off');
%ylim([-1,1])
print(gcf,['psth_population_activity'],'-dpng','-r0');

%%

cluster_folder_name = 'neuron_cluster';

if ~exist(cluster_folder_name, 'dir')
    mkdir(cluster_folder_name);
end

%% load the data and pre-process the data
time = (0:(size(NeuronSum1,2)-1))/10;

dataS = [NeuronSum1(:,time>=PCA_low_bound & time<=PCA_high_bound),NeuronSum2(:,time>=PCA_low_bound & time<=PCA_high_bound), NeuronSum3(:,time>=PCA_low_bound & time<=PCA_high_bound), NeuronSum4(:,time>=PCA_low_bound & time<=PCA_high_bound), NeuronSum5(:,time>=PCA_low_bound & time<=PCA_high_bound), NeuronSum6(:,time>=PCA_low_bound & time<=PCA_high_bound), NeuronSum7(:,time>=PCA_low_bound & time<=6)];
%not sure why he choose 2 and 8 here
%% PCA (reduce the time)
%Rows of X correspond to observations and columns to variables.
%Xiong appends random chuncks of different trial types, backwards
close all
dataS(1:3, :) = [];
[COEFF,score,~,~,explainedVar] = pca(dataS);
% SCORE*COEFF' to recapitulate the OG data
% 
scatter3(score(:,1),score(:,2),score(:,3))
axis equal
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')
%the coefficients are the eigenvalues, which are how much variance carried
%by each component. score is percentage of all eigenvalues
%score is 93
%COEFF, SCORE, LATENT, TSQUARED, EXPLAINED
explainedVar_3 = sum(explainedVar(1:3));
figure
bar(explainedVar(1:10));
title(['Explained Variance explained by first three PC: ',num2str(explainedVar_3), '%']);
ylabel('PC')
set(gca,'TickDir','Out','box','off');
D = score(:,1:size(score, 2)); % the first three principal components
print(gcf,[cluster_folder_name,'\pca_to_cluster'],'-dpdf','-r0');

%% Hierarchical clustering
% cgo = clustergram(D,'Standardize','Row');
distD = pdist(D, 'euclidean');%'complete'); % 'euclidean' Euclidean distance; 'correlation' One minus the sample correlation between points (treated as sequences of values)
Z = linkage(distD, 'complete');%'complete'); % 'average' 'ward' 'complete'%OG
%Z = linkage(distD);%'complete'); % 'average' 'ward' 'complete'
% cophenetic correlation coefficient to compare the results of clustering the same data set using different distance calculation methods or clustering algorithms
c = cophenet(Z,distD);
cutoff = median([Z(end-cluster_num+1,3) Z(end-cluster_num+2,3)]);

%Try clustering with kmeans
if do_kmeans
    kmeans_inds = kmeans(D, cluster_num);
end
T = cluster(Z,'maxclust',cluster_num);
%T = cluster(Z,'Cutoff',cutoff);
if do_kmeans
    T = kmeans_inds;
end

[sorted_inds, ind] = sort(T);
if do_kmeans
    display_clusters(sorted_inds);
end
D2 = D(ind,:);

DD = [];
DD{1} = NeuronSum1(ind,:);
DD{2} = NeuronSum2(ind,:);
DD{3} = NeuronSum3(ind,:);
DD{4} = NeuronSum4(ind,:);
DD{5} = NeuronSum5(ind,:);
DD{6} = NeuronSum6(ind,:);
DD{7} = NeuronSum7(ind,:);

figure('Position', [200 50 200 500],'Name','Cluster plot','numbertitle','off')
hc = subplot(1,3,1);
ColorMap = [linspace(0,1,64)',linspace(0,1,64)',linspace(0,1,64)'];
imagesc(D2);
colormap(hc,ColorMap)
set(gca,'xtick',[],'ytick',[],'box','off');
if ~do_kmeans
subplot(1,3,[2 3])
dendrogram(Z,0,'ColorThreshold',cutoff)
set(gca,'xtick',[],'ytick',[],'box','off');
view([90 90]) %90 -90 flipped vertically wrong way
end
% -90 90 flipped vertically and horizontally wrong

print(gcf,[cluster_folder_name,'\pca_cluster'],'-dpdf','-r0');

%% plot the pca of dendrogram
indd = [find(T==1);find(T==2);find(T==3);find(T==4);find(T==5);find(T==6);find(T==7)];
[~,denD,~,~,~] = pca(dataS(indd,:));

pca_crange = [-5,5];

figure('Position', [200 50 200 500],'Name','Cluster plot','numbertitle','off')
hc = subplot(1,3,1);
ColorMap = [linspace(0,1,64)',linspace(0,1,64)',linspace(0,1,64)'];
imagesc(denD(:,1:3),pca_crange);
%this plots the first 3 principal components
%colormap(hc,ColorMap)
%set(gca,'xtick',[],'ytick',[],'box','off');
if ~do_kmeans
subplot(1,3,[2 3])
dendrogram(Z,0,'ColorThreshold',cutoff)
set(gca,'xtick',[],'ytick',[],'box','off');
view([90 90])
end

print(gcf,[cluster_folder_name,'\pca_cluster_dendrogram_pca'],'-dpdf','-r0');

%%
figure('Position', [200 50 200 500],'Name','Cluster plot','numbertitle','off')
hc = subplot(1,3,1);
ColorMap = [linspace(0,1,64)',linspace(0,1,64)',linspace(0,1,64)'];
imagesc(denD(:,1:3),pca_crange);
colormap(hc,ColorMap)
set(gca,'xtick',[],'ytick',[],'box','off');
%c = colorbar('northoutside');
print(gcf,[cluster_folder_name,'\pca_cluster_dendrogram_pca_colorbar'],'-dpdf','-r0');

%%
ap_sort = [];
for k=1:cluster_num
    ind_temp = (T==k);
    ap_sort = [ap_sort;find((T==k))];
    
    AA = [];
    AA{1} = NeuronSum1(ind_temp,:);
    AA{2} = NeuronSum2(ind_temp,:);
    AA{3} = NeuronSum3(ind_temp,:);
    AA{4} = NeuronSum4(ind_temp,:);
    AA{5} = NeuronSum5(ind_temp,:);
    AA{6} = NeuronSum6(ind_temp,:);
    AA{7} = NeuronSum7(ind_temp,:);
    n_number = sum(ind_temp);
    
    figure
    hold on
    for ii=1:numel(AA)       
        Y = mean(AA{ii},1);
        Y = Y - mean(Y(5:39));
        Y_err = std(AA{ii})/sqrt(size(AA{ii},1));        
        AreaPlot(time,smooth(Y,3)',smooth(Y_err)',cc(ii,:),0.4,1);
    end
    xlabel('Time (s)'); ylabel('Normalized activity');
    title(['Neuron number ',num2str(n_number)]);    
    set(gca,'TickDir', 'out','xlim',[2,10],'xtick',2:2:10,'FontSize', 12,'box','off');
    print(gcf,[cluster_folder_name,'\psth_cluster',num2str(k)],'-dpdf','-r0');
%     close;
    
end

%%
% Plot heatmap of ranked cluster
R = 28; % for colorspace
C = 28; % for colorspace
ColorMap = zeros(64,3);
ColorMap(64-R+1:64,1) = linspace(0,1,R); % R
ColorMap(1:C,2) = linspace(1,0,C); % G
ColorMap(1:C,3) = linspace(1,0,C); % B

ap_sort = [find(T==1);find(T==2);find(T==3);find(T==4);find(T==5);find(T==6);find(T==7)];

NeuronSum1_sort = NeuronSum1(ap_sort,:);
NeuronSum2_sort = NeuronSum2(ap_sort,:);
NeuronSum3_sort = NeuronSum3(ap_sort,:);
NeuronSum4_sort = NeuronSum4(ap_sort,:);
NeuronSum5_sort = NeuronSum5(ap_sort,:);
NeuronSum6_sort = NeuronSum6(ap_sort,:);
NeuronSum7_sort = NeuronSum7(ap_sort,:);

figure
subplot(1,7,1)
imagesc(time,1:size(NeuronSum1,1),NeuronSum1_sort,[-3 3]);
set(gca,'TickDir', 'out','xlim',[2,10],'xtick',[2:2:10],...
    'FontSize', 16,'box','off');
colormap(ColorMap)

subplot(1,7,2)
imagesc(time,1:size(NeuronSum2,1),NeuronSum2_sort,[-3 3]);
set(gca,'TickDir', 'out','xlim',[2,10],'xtick',[2:2:10],...
    'FontSize', 16,'box','off');
colormap(ColorMap)

subplot(1,7,3)
imagesc(time,1:size(NeuronSum3,1),NeuronSum3_sort,[-3 3]);
set(gca,'TickDir', 'out','xlim',[2,10],'xtick',[2:2:10],...
    'FontSize', 16,'box','off');
colormap(ColorMap)

subplot(1,7,4)
imagesc(time,1:size(NeuronSum4,1),NeuronSum4_sort,[-3 3]);
set(gca,'TickDir', 'out','xlim',[2,10],'xtick',[2:2:10],...
    'FontSize', 16,'box','off');
colormap(ColorMap)

subplot(1,7,5)
imagesc(time,1:size(NeuronSum5,1),NeuronSum5_sort,[-3 3]);
set(gca,'TickDir', 'out','xlim',[2,10],'xtick',[2:2:10],...
    'FontSize', 16,'box','off');
colormap(ColorMap)

subplot(1,7,6)
imagesc(time,1:size(NeuronSum6,1),NeuronSum6_sort,[-3 3]);
set(gca,'TickDir', 'out','xlim',[2,10],'xtick',[2:2:10],...
    'FontSize', 16,'box','off');
colormap(ColorMap)

subplot(1,7,7)
imagesc(time,1:size(NeuronSum7,1),NeuronSum7_sort,[-3 3]);
set(gca,'TickDir', 'out','xlim',[2,10],'xtick',[2:2:10],...
    'FontSize', 16,'box','off');
colormap(ColorMap)
print(gcf,['psth_heatmap'],'-dpdf','-r0');
