output_dir = 'D:\BSc Project\MAIN\10\filtered_data\CAR data\';  
output_files = dir(fullfile(output_dir, 'CAR_*.set'));  
output_files = {output_files.name};  

for i = 1:length(output_files)
    EEG = pop_loadset('filename', output_files{i}, 'filepath', output_dir);
    EEGclean = pop_cleanline(EEG, 'LineFrequencies', 60, 'Bandwidth', 2, ...
                             'ChanCompIndices', 1:18, 'SignalType', 'Channels', ...
                             'ComputeSpectralPower', false, 'ScanForLines', false, ...
                             'PlotFigures', false, 'VerboseOutput', false);
    output_clean_filename = ['cleanLined_' output_files{i}]; 
    pop_saveset(EEGclean, 'filename', output_clean_filename, 'filepath', output_dir);
end
