clear, close all, clc 
% EEG = pop_loadset('after3Steps.set') ;
EEG = pop_loadset('chb01_final_withoutCleanLine.set') ;

fs = EEG.srate ;
WinLength = 30 ; % 30-sec windows
WinSamples = WinLength*fs ;
NumWin = size(EEG.data,2)/WinSamples ;
thresh = [0 0;0.8 1; 0.8 1; 0.8 1; 0.8 1; 0.8 1; 0.8 1] ; % Class = { 'Brain'  'Muscle'  'Eye'  'Heart'  'Line Noise'  'Channel Noise'  'Other' }
EEG1 = EEG.data ; 
for winnum = 1:NumWin
    newEEG = pop_select(EEG,'time',[(winnum-1)*WinSamples,winnum*WinSamples-1]/fs) ;
    % [newEEG] = pop_runica( newEEG, 'icatype','fastica','numOfIC', 17 ) ;
    [newEEG] = pop_runica( newEEG, 'icatype','runica','pca', 17 );
    % newEEG = pop_selectcomps(newEEG) ;
    [newEEG] = iclabel(newEEG) ;
    % Class = { 'Brain'  'Muscle'  'Eye'  'Heart'  'Line Noise'  'Channel Noise'  'Other' }
    [newEEG] = pop_icflag(newEEG, thresh) ;
    [newEEG,com] = pop_subcomp(newEEG,[],0) ;
    EEG.data(:,(winnum-1)*WinSamples+1:winnum*WinSamples) = newEEG.data ;
end
EEG2 = EEG.data ;
SaveFileName = 'chb01_03_main_ICA_denoised' ;
SavePathName = 'D:\' ;
pop_saveset(EEG,'filename',SaveFileName,'filepath',SavePathName,'savemode','onefile');



