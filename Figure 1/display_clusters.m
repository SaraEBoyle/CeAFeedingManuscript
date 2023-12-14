% Make a plot of the sorted inds
%vertical colored lines for each cluster proportional to the 
%Size of the cluster
function [] = display_clusters(sorted_inds)
first = size(find(sorted_inds == 6), 1);
second = size(find(sorted_inds == 5), 1);
third = size(find(sorted_inds == 4), 1);
fourth = size(find(sorted_inds == 3), 1);
fifth = size(find(sorted_inds == 2), 1);
sixth = size(find(sorted_inds == 1), 1);

figure;
hold on
plot([0, 0], [0, first], 'color', 'g', 'LineWidth', 8);
plot([0, 0], [first, first + second], 'color', 'b', 'LineWidth', 8);
plot([0, 0], [first + second, first + second + third], 'color', 'm', 'LineWidth', 8);
plot([0, 0], [first + second + third, fourth + first + second + third], 'color', 'c', 'LineWidth', 8);
plot([0, 0], [first + second + third + fourth, fourth + first + second + third + fifth], 'color', 'r', 'LineWidth', 8);
plot([0, 0], [fifth + first + second + third + fourth, fourth + first + second + third + fifth + sixth], 'color', 'y', 'LineWidth', 8);
