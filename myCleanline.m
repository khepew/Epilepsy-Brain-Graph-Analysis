
EEG = pop_loadset('filename', 'CAR_chb07_12_filtered.set', 'filepath', 'D:\BSc Project\MAIN\10\filtered_data\CAR data\');

EEGclean = pop_cleanline(EEG, 'LineFrequencies', 60, 'Bandwidth', 2, ...
                         'ChanCompIndices', 1:3, 'SignalType', 'Channels', ...
                         'ComputeSpectralPower', false, 'ScanForLines', false, ...
                         'PlotFigures', false, 'VerboseOutput', false);


pop_saveset(EEGclean, 'filename', 'MEOWWWWWWWWWWWWWWWWWWWW.set', 'filepath', 'D:\BSc Project\MAIN\10\filtered_data\CAR data\');
