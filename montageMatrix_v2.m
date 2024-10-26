% Specify the directory containing the filtered files
input_dir = 'D:\BSc Project\MAIN\10\filtered_data\';

% Get a list of all filtered .set files in the directory
filtered_files = dir(fullfile(input_dir, '*_filtered.set'));
filtered_files = {filtered_files.name};  % Extract file names into a cell array

% Define the transformation matrix A
A = zeros(23, 21);
A(1, 1) = 1; A(1, 2) = -1;
A(2, 2) = 1; A(2, 3) = -1;
A(3, 3) = 1; A(3, 4) = -1;
A(4, 4) = 1; A(4, 5) = -1;
A(5, 1) = 1; A(5, 6) = -1;
A(6, 6) = 1; A(6, 7) = -1;
A(7, 7) = 1; A(7, 8) = -1;
A(8, 8) = 1; A(8, 5) = -1;
A(9, 9) = 1; A(9, 10) = -1;
A(10, 10) = 1; A(10, 11) = -1;
A(11, 11) = 1; A(11, 12) = -1;
A(12, 12) = 1; A(12, 13) = -1;
A(13, 9) = 1; A(13, 14) = -1;
A(14, 14) = 1; A(14, 15) = -1;
A(15, 15) = 1; A(15, 16) = -1;
A(16, 16) = 1; A(16, 13) = -1;
A(17, 17) = 1; A(17, 18) = -1;
A(18, 18) = 1; A(18, 19) = -1;
A(19, 4) = 1; A(19, 3) = -1;
A(20, 3) = 1; A(20, 20) = -1;
A(21, 20) = 1; A(21, 21) = -1;
A(22, 21) = 1; A(22, 15) = -1;
A(23, 15) = 1; A(23, 16) = -1;

for i = 1:length(filtered_files)
    EEG = pop_loadset('filename', filtered_files{i}, 'filepath', input_dir);
    val = EEG.data;
    newval = val([1:16, 19:end], :);  % Remove Fz, Pz, Cz
    newA = A([1:16, 19:end], [1:16, 20:21]);  
    newX = pinv(newA) * newval;  
    EEG.nbchan = 18;  
    chanlocs = struct('labels', { 'FP1', 'F7', 'T7', 'P7', 'O1', 'F3', 'C3', 'P3', ...
                                   'FP2', 'F4', 'C4', 'P4', 'O2', 'F8', 'T8', 'P8', ...
                                   'FT9', 'FT10' });
    
    EEG.chanlocs = chanlocs;  
    EEG.data = newX; 
    EEG.setname = 'AfterStepsOneTwo'; 
    output_filename = ['CAR_' filtered_files{i}(1:end-4) '.set'];  
    pop_saveset(EEG, 'filename', output_filename, 'filepath', input_dir);
end
