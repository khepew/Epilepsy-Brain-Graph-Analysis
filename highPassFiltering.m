files = {'chb10_12.edf', 'chb10_20.edf', 'chb10_27.edf', 'chb10_30.edf', 'chb10_31.edf', 'chb10_38.edf', 'chb10_89.edf', ...
    'chb07_12.edf', 'chb07_13.edf', 'chb07_19.edf', 'chb08_02.edf', 'chb08_05.edf', 'chb08_11.edf', 'chb08_13.edf', 'chb08_21.edf'};

output_dir = './filtered_data/';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

for i = 1:length(files)
    EEG = pop_biosig(files{i});
    EEG = pop_eegfiltnew(EEG, 0.5, []);
    output_filename = [files{i}(1:end-4) '_filtered.set'];
    EEG = pop_saveset(EEG, 'filename', output_filename, 'filepath', output_dir);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
end