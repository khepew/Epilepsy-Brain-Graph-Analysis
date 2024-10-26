clc; clear;
seizure_onsets = [
    2996, 1467, 1732, 1015, 1720, 327, 1862, 4920, 3285, 13688, ...
    2670, 2856, 2988, 2417, 2083, 6313, 6888, 2382, 3021, 3801, ...
    4618, 1383
];
seizure_offsets = [
    3036, 1494, 1772, 1066, 1810, 420, 1963, 5006, 3381, 13831, ...
    2841, 3046, 3122, 2577, 2347, 6348, 6958, 2447, 3079, 3877, ...
    4707, 1437
];
seizure_durations = [
    40, 27, 40, 51, 90, 93, 101, 86, 96, 143, ...
    171, 190, 134, 160, 264, 35, 70, 65, 58, 76, ...
    89, 54
];
file_names = {'chb01_03', 'chb01_04', 'chb01_15', 'chb01_16', 'chb01_18', 'chb01_21', 'chb01_26', ... 
    'chb07_12', 'chb07_13', 'chb07_19', ...
    'chb08_02', 'chb08_05', 'chb08_11', 'chb08_13', 'chb08_21', ... 
    'chb10_12', 'chb10_20','chb10_27','chb10_30','chb10_31','chb10_38','chb10_89',
    };
% 
for i = 1:22
    onset = seizure_onsets(i);
    offset = seizure_offsets(i);
    L = seizure_durations(i);  
    EEG = pop_loadset(['cleanLined_CAR_' file_names{i} '_filtered.set']); 
    % Extract the three phases
    L = (L + 1) * EEG.srate;
    onsetSamples = onset * EEG.srate;
    offsetSamples = offset * EEG.srate;
    preictal_data = EEG.data(:, (onsetSamples - L - (EEG.srate - 1)) : (onsetSamples - EEG.srate));
    ictal_data = EEG.data(:, (onsetSamples - (EEG.srate - 1)) : offsetSamples);
    postictal_data = EEG.data(:, (offsetSamples + 1) : (offsetSamples + L));
    durationFloor = floor(seizure_durations(i)/10) * 10;
    combined_data = [preictal_data(:, 1:durationFloor*EEG.srate), ictal_data(:, 1:durationFloor*EEG.srate), postictal_data(:, 1:durationFloor*EEG.srate)];
    EEG_combined = EEG;  
    EEG_combined.data = combined_data;  
    EEG_combined.pnts = size(combined_data, 2);
    EEG_combined.times = linspace(EEG.xmin, EEG.xmin + (EEG_combined.pnts - 1) / EEG.srate, EEG_combined.pnts);
    pop_saveset(EEG_combined, 'filename', ['combined_' file_names{i} '.set']);
end
