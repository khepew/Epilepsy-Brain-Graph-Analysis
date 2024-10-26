clc; clear;
load("LearnedGraphs_AR.mat");
load("preictal_vs_ictal_ar.mat");
%%
preictal_edge_1 = squeeze(learnedGraphs_preictal(:, 15, 6)); 
ictal_edge_1 = squeeze(learnedGraphs_ictal(:, 15, 6));      
preictal_edge_2 = squeeze(learnedGraphs_preictal(:, 17, 4)); 
ictal_edge_2 = squeeze(learnedGraphs_ictal(:, 17, 4));       
if length(preictal_edge_1) ~= length(preictal_edge_2) || ...
   length(ictal_edge_1) ~= length(ictal_edge_2) || ...
   length(preictal_edge_1) ~= length(ictal_edge_1)
    error('Vectors must be of the same length.');
end
figure;
hold on;
scatter3(preictal_edge_1, preictal_edge_2, 1:412, 'b', 'DisplayName', 'Preictal Edges (15,6) and (17,4)');
scatter3(ictal_edge_1, ictal_edge_2, 1:412, 'r', 'DisplayName', 'Ictal Edges (15,6) and (17,4)');
xlabel('Edge Weight (15,6)');
ylabel('Edge Weight (17,4)');
zlabel('Window');
title('3D Visualization of Edges (15,6) and (17,4) in Preictal and Ictal States');
legend('show');
hold off;
view(45, 30); 



