% ch1 FP1
% ch2 F7
% ch3 T7
% ch4 P7
% ch5 O1
% ch6 F3
% ch7 C3
% ch8 P3
% ch9 FP2
% ch10 F4
% ch11 C4
% ch12 P4
% ch13 O2
% ch14 F8
% ch15 T8
% ch16 P8
% ch17 Fz
% ch18 Cz
% ch19 Pz
% ch20 FT9
% ch21 FT10

A = zeros(23,21) ;
A(1,1) = 1 ; 
A(1,2) = -1 ;
A(2,2) = 1 ;
A(2,3) = -1 ;
A(3,3) = 1 ;
A(3,4) = -1 ;
A(4,4) = 1 ;
A(4,5) = -1 ;
A(5,1) = 1 ;
A(5,6) = -1 ;
A(6,6) = 1 ;
A(6,7) = -1 ;
A(7,7) = 1 ;
A(7,8) = -1 ;
A(8,8) = 1 ;
A(8,5) = -1 ;
A(9,9) = 1 ;
A(9,10) = -1 ;
A(10,10) = 1 ;
A(10,11) = -1 ;
A(11,11) = 1 ;
A(11,12) = -1 ;
A(12,12) = 1 ;
A(12,13) = -1 ;
A(13,9) = 1 ;
A(13,14) = -1 ;
A(14,14) = 1 ;
A(14,15) = -1 ;
A(15,15) = 1 ;
A(15,16) = -1 ;
A(16,16) = 1 ;
A(16,13) = -1 ;
A(17,17) = 1 ;
A(17,18) = -1 ;
A(18,18) = 1 ;
A(18,19) = -1 ;
A(19,4) = 1 ;
A(19,3) = -1 ;
A(20,3) = 1 ;
A(20,20) = -1 ;
A(21,20) = 1 ;
A(21,21) = -1 ;
A(22,21) = 1 ;
A(22,15) = -1 ;
A(23,15) = 1 ;
A(23,16) = -1 ;

% Re-reference & remove Fz, Pz, Cz
EEG = pop_loadset('filename', 'chb01_03_afterStep1.set', 'filepath', 'D:\BSc Project\chb01_prep\');
val = EEG.data ;
newval = val([1:16,19:end],:) ;
newA = A([1:16,19:end],[1:16,20:21]) ;
newX = pinv(newA)*newval ;

%%% chanlocs
EEG.nbchan = 18; 
% EEG.chanlocs = pop_readlocs('D:\BSc Project\chb01_prep\my18Chan.ced');

chanlocs = struct('labels', { 'FP1', 'F7', 'T7', 'P7', 'O1', 'F3', 'C3', 'P3', 'FP2', 'F4', 'C4', 'P4', 'O2', 'F8', 'T8', 'P8', 'FT9', 'FT10'});
pop_chanedit( chanlocs );

EEG.chanlocs = chanlocs;
EEG.data = newX;  
EEG.setname = 'AfterStepsOneTwo';


pop_saveset(EEG, 'filename', 'new_dataset_CAR.set', 'filepath', 'D:\BSc Project\chb01_prep\');







% % 
% [ seizure_start_time_offset_in_seconds, seizure_length_in_seconds ] = get_seizure_period( 'chb01_03.edf.seizures' ) ;
% 
% Seizure_start = seizure_start_time_offset_in_seconds*256 ;
% Seizure_end =  seizure_length_in_seconds*256 + Seizure_start ;
% ChannelNames = {'FP1','F7','T7','P7','O1','F3','C3','P3','FP2','F4','C4','P4','O2','F8','T8','P8','FT9','FT10'} ;
% save('chb01_03.mat','ChannelNames','newX','Seizure_end','Seizure_start') ;
