clc;
clear;
fs = 256;
WinLength = 10; % 10-sec windows
WinSamples = WinLength * fs;
file_names = {'chb01_03', 'chb01_04', 'chb01_15', 'chb01_16', 'chb01_18', 'chb01_21', 'chb01_26', ... 
    'chb07_12', 'chb07_13', 'chb07_19', ...
    'chb08_02', 'chb08_05', 'chb08_11', 'chb08_13', 'chb08_21', ... 
    'chb10_12', 'chb10_20','chb10_27','chb10_30','chb10_31','chb10_38','chb10_89',
    };
for i = 1:22
    EEG = pop_loadset(['combined_' file_names{i} '.set']);
    NumWin = size(EEG.data,2)/WinSamples ;
    thresh = [0 0;0.8 1; 0.8 1; 0.8 1; 0.8 1; 0.8 1; 0.8 1] ; % Class = { 'Brain'  'Muscle'  'Eye'  'Heart'  'Line Noise'  'Channel Noise'  'Other' }
    EEG1 = EEG.data ;
    for winnum = 1:NumWin
        newEEG = pop_select(EEG,'time',[(winnum-1)*WinSamples,winnum*WinSamples-1]/fs) ;
        [newEEG] = pop_runica( newEEG, 'icatype','runica','pca', 17 );
        [newEEG] = iclabel(newEEG) ;
        % Class = { 'Brain'  'Muscle'  'Eye'  'Heart'  'Line Noise'  'Channel Noise'  'Other' }
        [newEEG] = pop_icflag(newEEG, thresh) ;
        [newEEG,com] = pop_subcomp(newEEG,[],0) ;
        EEG.data(:,(winnum-1)*WinSamples+1:winnum*WinSamples) = newEEG.data ;
    end
    EEG2 = EEG.data ;
    SaveFileName = [file_names{i} '_afterPrep.set']; 
    SavePathName = 'D:\BSc Project\MAIN\10\filtered_data\CAR data\HPF_CAR_CLINED\Extracted_Phases\';
    pop_saveset(EEG, 'filename', SaveFileName, 'filepath', SavePathName, 'savemode', 'onefile');
end


