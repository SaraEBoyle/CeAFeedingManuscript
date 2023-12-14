function [] = cluster_based_on_correlation_matrix()
%% Takes a response matrix n neurons by m stimuli
%% creates a stimulus correlation matrix and uses 1 - correlation as the distance
load('single_neuron_avg_response_matrices.mat');
titles = response_matrix(1, :);
neuron_reponses = cell2mat(response_matrix(2:end, :));
R = corrcoef(neuron_reponses);
flat_correlation_distances = [];
for row = 1:size(R, 1)
    for col = (row + 1):size(R, 2)
        flat_correlation_distances = horzcat(flat_correlation_distances, 1 - R(row, col));
    end
end
%cosD = pdist(cell2mat(neuron_reponses),'correlation');
%clustTreeCos = linkage(cosD,'complete');
%cophenet(clustTreeCos,cosD)
%[h,nodes] = dendrogram(clustTreeCos,0);
%h_gca = gca;
%h_gca.TickDir = 'out';
%h_gca.TickLength = [.002 0];
%T = cluster(clustTreeCos,'maxclust',4);

% This means the bigger the positive correlation, the less distance
%D rows are neurons, columns are PCAs
%distD = pdist(D, 'euclidean');%'complete'); % 'euclidean' Euclidean distance; 'correlation' One minus the sample correlation between points (treated as sequences of values)
%   D is a M-by-N matrix.
%   1-by-(M*(M-1)/2) row vector, corresponding to the M*(M-1)/2 pairs of
%   observations in X.
%The output D is arranged in the order of ((2,1),(3,1),..., (M,1),
%   (3,2),...(M,2),.....(M,M-1)), i.e. the lower left triangle of the full
%   M-by-M distance matrix in column order.  To get the distance between
%   the Ith and Jth observations (I < J), either use the formula
%   D((I-1)*(M-I/2)+J-I), or use the helper function Z = SQUAREFORM(D),
%   which returns an M-by-M square symmetric matrix, with the (I,J) entry
%   equal to distance between observation I and observation J.
cluster_num = 4;

% fat, sucrose, quinine, mineral oil, xanthan gum, water, air
%cor_distances = [1 - .8284, 1 - .4811, 1 - .3141, 1 - .4857, 1 - .42, 1 - .0935, 1 - .5166, ...
%    1 - .2048, 1 - .4426, 1 - .4080, 1 - .0673, 1 - .2595, 1 - .4663, 1 - .5104, 1 - -.0129, ...
%    1 - .4121, 1 - .3714, 1 - -.0247, 1 - .9606, 1 - -.2332, 1 - -.1931];

cor_distances = flat_correlation_distances;
%distD is every pairs distance without replacement
Z = linkage(cor_distances, 'complete');%'complete'); % 'average' 'ward' 'complete'
% cophenetic correlation coefficient to compare the results of clustering the same data set using different distance calculation methods or clustering algorithms
c = cophenet(Z,cor_distances);
cutoff = median([Z(end-cluster_num+1,3) Z(end-cluster_num+2,3)]);

%T = cluster(Z,'maxclust',cluster_num);

%[~, ind] = sort(T);
%D2 = D(ind,:);

%DD = [];
%DD{1} = NeuronSum1(ind,:);
%DD{2} = NeuronSum2(ind,:);
%DD{3} = NeuronSum3(ind,:);
%DD{4} = NeuronSum4(ind,:);
%DD{5} = NeuronSum5(ind,:);
%DD{6} = NeuronSum6(ind,:);
%DD{7} = NeuronSum7(ind,:);

%figure('Position', [200 50 200 500],'Name','Cluster plot','numbertitle','off')
%hc = subplot(1,3,1);
%ColorMap = [linspace(0,1,64)',linspace(0,1,64)',linspace(0,1,64)'];
%imagesc(D2);
%colormap(hc,ColorMap)
%set(gca,'xtick',[],'ytick',[],'box','off');

%subplot(1,3,[2 3])
%dendrogram(Z,0,'ColorThreshold',cutoff)
%set(gca,'xtick',[],'ytick',[],'box','off');
%view([90 -90])

%print(gcf,[cluster_folder_name,'\pca_cluster'],'-dpdf','-r0');

%% plot the pca of dendrogram
%indd = [find(T==1);find(T==2);find(T==3);find(T==4);find(T==5);find(T==6);find(T==7)];
%[~,denD,~,~,~] = pca(dataS(indd,:));

%pca_crange = [-10,10];

%figure('Position', [200 50 200 500],'Name','Cluster plot','numbertitle','off')
%hc = subplot(1,3,1);
%ColorMap = [linspace(0,1,64)',linspace(0,1,64)',linspace(0,1,64)'];
%imagesc(denD(:,1:3),pca_crange);
%%this plots the first 3 principal components
%colormap(hc,ColorMap)
%set(gca,'xtick',[],'ytick',[],'box','off');

figure;
[H,T,OUTPERM] = dendrogram(Z, 'Orientation','left','ColorThreshold',.55);
%[H,T,OUTPERM] = dendrogram(Z,0,'ColorThreshold',cutoff);
%set(gca,'xtick',[],'ytick',[],'box','off');
%view([90 -90])

end