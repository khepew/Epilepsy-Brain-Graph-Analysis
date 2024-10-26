clc;
clear;
Fs = 256;
windowLength = 5 * Fs;
file_names = {'chb01_03', 'chb01_04', 'chb01_15', 'chb01_16', 'chb01_18', 'chb01_21', 'chb01_26', ... 
    'chb07_12', 'chb07_13', 'chb07_19', ...
    'chb08_02', 'chb08_05', 'chb08_11', 'chb08_13', 'chb08_21', ... 
    'chb10_12', 'chb10_20','chb10_27','chb10_30','chb10_31','chb10_38','chb10_89',
    };
globalPreictalWindows = [];
globalIctalWindows = [];
globalPostictalWindows = [];
for i = 1:22
    EEG = pop_loadset([file_names{i} '_afterPrep.set']);
    numWindows = (size(EEG.data,2) / (3 * windowLength));
    preictalWindows = reshape(EEG.data(:, 1:numWindows * windowLength), 18, windowLength, numWindows);
    ictalWindows = reshape(EEG.data(:, ((numWindows * windowLength) + 1): (2 * numWindows * windowLength)), 18, windowLength, numWindows);
    postictalWindows = reshape(EEG.data(:, ((2 * numWindows * windowLength) + 1) : (3 * numWindows * windowLength)), 18, windowLength, numWindows);
    globalPreictalWindows = cat(3, globalPreictalWindows, preictalWindows);
    globalIctalWindows = cat(3, globalIctalWindows, ictalWindows);
    globalPostictalWindows = cat(3, globalPostictalWindows, postictalWindows);
end

save('Global_Phase_Windows.mat', 'globalPreictalWindows', 'globalIctalWindows', 'globalPostictalWindows');

