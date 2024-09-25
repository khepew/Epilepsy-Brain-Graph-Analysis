clc;
clear all;
freq = 256;
load('Example1.mat');
data_dimensions = size(newX);
numberOfChannels = data_dimensions(1);
numberOfSamples = data_dimensions(2);
time = (1:numberOfSamples) / freq;
channel_names = {'F7', 'T7', 'P7', 'O1', 'O2', 'F8', 'T8', 'P8'};
% extracting coords
fieldsToExtract = {'X', 'Y'};
x_coords = zeros(1, numel(channel_names));
y_coords = zeros(1, numel(channel_names));
for i = 1:numel(channel_names)
    channelName = channel_names{i};
    rowIndex = find(strcmp({ChanLocs.labels}, channelName));
    x_coords(i) = ChanLocs(rowIndex).X;
    y_coords(i) = ChanLocs(rowIndex).Y;
end

%% Part Alef

window_length = 5 * freq;

% number of windows for each interval
wnum_preictals = floor((Seizure_start) / window_length) - 2;
wnum_ictals = floor((Seizure_end - Seizure_start + 1) / window_length) - 2;
wnum_postictals = floor((numberOfSamples - Seizure_end + 1) / window_length) - 2;

% initializing cells
cells_preictals = cell(1, wnum_preictals);
cells_ictals = cell(1, wnum_ictals);
cells_postictals = cell(1, wnum_postictals);

% Extracting Data
% PRE-ICTAL
for i = 1:wnum_preictals
    start_idx = (i) * window_length + 1;
    end_idx = (i + 1) * window_length;
    cells_preictals{i} = newX(:, start_idx:end_idx);
end
% ICTAL
%ictal_newX = newX(:,Seizure_start:Seizure_end);
for i = 1: wnum_ictals
    start_idx = (i) * window_length + 1;
    end_idx = (i + 1) * window_length;
    cells_ictals{i} = newX(:, start_idx:end_idx);
end
% POST-ICTAL
%postictal_newX = newX(:,Seizure_end:end);
for i = 1:wnum_postictals
    start_idx = (i) * window_length + 1;
    end_idx = (i + 1) * window_length;
    cells_postictals{i} = newX(:, start_idx:end_idx);
end

%% Be 
modelOrder = 10;
k = 6;

% Learning Graphs
learnedGraphs_preictals = compute_learningGraphs(cells_preictals, numberOfChannels, modelOrder, k);
learnedGraphs_ictals = compute_learningGraphs(cells_ictals, numberOfChannels, modelOrder, k);
learnedGraphs_postictals = compute_learningGraphs(cells_postictals, numberOfChannels, modelOrder, k);


%% Part Jim

meaningful_edges_pre_ict = find_meaningful_edges(learnedGraphs_preictals, learnedGraphs_ictals);
meaningful_edges_pre_post = find_meaningful_edges(learnedGraphs_preictals, learnedGraphs_postictals);
meaningful_edges_ict_post = find_meaningful_edges(learnedGraphs_ictals, learnedGraphs_postictals);



%% Part Ve

meaningful_edges2_pre_ict = find_meaningful_edges2(learnedGraphs_preictals, learnedGraphs_ictals);
meaningful_edges2_pre_post = find_meaningful_edges2(learnedGraphs_preictals, learnedGraphs_postictals);
meaningful_edges2_ict_post = find_meaningful_edges2(learnedGraphs_ictals, learnedGraphs_postictals);


%%
similarity_percentage = calculate_similarity(meaningful_edges_pre_ict, meaningful_edges2_pre_ict);
disp(['Similarity (Pre - Ict): ', num2str(similarity_percentage), '%']);

similarity_percentage = calculate_similarity(meaningful_edges_pre_post, meaningful_edges2_pre_post);
disp(['Similarity (Pre - Post): ', num2str(similarity_percentage), '%']);

similarity_percentage = calculate_similarity(meaningful_edges_ict_post, meaningful_edges2_ict_post);
disp(['Similarity (Ict - Post): ', num2str(similarity_percentage), '%']);



%%%%%%%%%%%%% Functions %%%%%%%%%%%%%
function learnedGraphs = compute_learningGraphs(selected_data, numberOfChannels, modelOrder, k)
    % AR coefficients for (X)
    modelAR_coefs = zeros(length(selected_data), numberOfChannels, modelOrder);
    for i = 1:length(selected_data)
        wind = selected_data(i);
        coefMatrix = zeros(numberOfChannels, modelOrder);
        % Fit autoregressive models for each feature
        for j = 1:numberOfChannels
            currentChannel = wind{1}(j, :);
            model = ar(currentChannel, modelOrder);
            coefMatrix(j, :) = model.A(:, 2:end);
        end
        modelAR_coefs(i, :, :) = coefMatrix;
    end
    learnedGraphs = zeros(length(selected_data), numberOfChannels, numberOfChannels);
    % Computing matrix Z
    for i = 1:length(selected_data)
        X = modelAR_coefs(i, :, :);
        X = squeeze(X);
        Z = gsp_distanz(X').^2;
        theta = gsp_compute_graph_learning_theta(Z,k);
        [W, ~] = gsp_learn_graph_log_degrees(theta * Z, 1, 1);
        % W(W<1e-4) = 0;
        learnedGraphs(i, :, :) = W;
    end
end

function meaningful_edges = find_meaningful_edges(learnedGraphs_preictals, learnedGraphs_ictals)
    [~, numberOfChannels, ~] = size(learnedGraphs_preictals);
    meaningful_edges = zeros(numberOfChannels, numberOfChannels);
    counter = 0;
    all_pvalues = zeros(1, nchoosek(numberOfChannels, 2)); 

    for i = 1:numberOfChannels
        for j = i+1:numberOfChannels
            counter = counter + 1;
            preictal_weights = squeeze(learnedGraphs_preictals(:, i, j));
            ictal_weights = squeeze(learnedGraphs_ictals(:, i, j));
            
            % T-TEST
            [~, p] = ttest2(preictal_weights, ictal_weights);
            all_pvalues(counter) = p;
        end
    end
    
    % Correction
    p_adjusted = mafdr(all_pvalues, 'BHFDR', true);
    % meaningful_edges_indices = find(p_adjusted<0.05);
    % meaningful_edges(meaningful_edges_indices) = 1;
    
    % meaningful edges based on p-values
    counter = 0;
    for i = 1:numberOfChannels
        for j = i+1:numberOfChannels
            counter = counter + 1;
            if p_adjusted(counter) < 0.05
                meaningful_edges(i, j) = 1;
                meaningful_edges(j, i) = 1;
            end
        end
    end
end

function meaningful_edges = find_meaningful_edges2(learnedGraphs_preictals, learnedGraphs_ictals)
    [~, numberOfChannels, ~] = size(learnedGraphs_preictals);
    meaningful_edges = zeros(numberOfChannels, numberOfChannels);
    counter = 0;
    all_pvalues = zeros(1, nchoosek(numberOfChannels, 2)); 

    for i = 1:numberOfChannels
        for j = i+1:numberOfChannels
            counter = counter + 1;
            preictal_weights = squeeze(learnedGraphs_preictals(:, i, j));
            ictal_weights = squeeze(learnedGraphs_ictals(:, i, j));
            
            % T-TEST
            [~, p] = ranksum(preictal_weights, ictal_weights);
            all_pvalues(counter) = p;
        end
    end
    
    % Correction
    p_adjusted = mafdr(all_pvalues, 'BHFDR', true);
    % meaningful_edges_indices = find(p_adjusted<0.05);
    % meaningful_edges(meaningful_edges_indices) = 1;
    
    % meaningful edges based on p-values
    counter = 0;
    for i = 1:numberOfChannels
        for j = i+1:numberOfChannels
            counter = counter + 1;
            if p_adjusted(counter) < 0.05
                meaningful_edges(i, j) = 1;
                meaningful_edges(j, i) = 1;
            end
        end
    end
end

function similarity_percentage = calculate_similarity(matrix1, matrix2)

    matching_elements = sum(matrix1(:) == matrix2(:));
    total_elements = numel(matrix1);
    similarity_percentage = (matching_elements / total_elements) * 100;
end











