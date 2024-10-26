clear;
clc;
load("channel_Names.mat");
load("X_chLocs.mat");
load("Y_chLocs.mat");
%%
load('preictal_vs_ictal_ar.mat');
load("preictal_vs_postictal_ar.mat");
load("ictal_vs_postictal_ar.mat");
%%
load('preictal_vs_ictal_power.mat');
load("preictal_vs_postictal_power.mat");
load("ictal_vs_postictal_power.mat");
%%
significant_edges_matrices = {
    final_significant_edges12, ...
    final_significant_edges13, ...
    final_significant_edges23
};
significance_directions_matrices = {
    final_significance_direction12, ...
    final_significance_direction13, ...
    final_significance_direction23
};
titles = {
    'Preictal vs. Ictal (AR for k = 8)', ...
    'Preictal vs. Postictal (AR for k = 8)', ...
    'Ictal vs. Postictal (AR for k = 8)'
};
save_names = {
    'graph_ar_pre_ict.fig', ...
    'graph_ar_pre_post.fig', ...
    'graph_ar_ict_post.fig'
};
folder_name = 'k8_ar';
mkdir(folder_name);
for k = 1:length(significant_edges_matrices)
    G = graph(significant_edges_matrices{k}); 
    figure;
    edgeIndices = G.Edges.EndNodes;
    edgeColors = repmat([0.5, 0.5, 0.5], size(edgeIndices, 1), 1); 
    for e = 1:size(edgeIndices, 1)
        node1 = edgeIndices(e, 1);
        node2 = edgeIndices(e, 2);
        if significant_edges_matrices{k}(node1, node2) ~= 0
            direction = significance_directions_matrices{k}{node1, node2};
            if strcmp(direction, 'increased')
                edgeColors(e, :) = [1, 0, 0]; 
            elseif strcmp(direction, 'decreased')
                edgeColors(e, :) = [0, 0, 1]; 
            end
        end
    end
    h = plot(G, 'XData', -yCoords, 'YData', xCoords, 'NodeLabel', chanLabels, 'EdgeColor', edgeColors);
    title(titles{k});
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    full_path = fullfile(folder_name, save_names{k});  
    savefig(full_path);  
end

