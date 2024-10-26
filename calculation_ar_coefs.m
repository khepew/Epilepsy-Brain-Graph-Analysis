clc;
clear;
load("Global_Phase_Windows.mat");
numChannels = size(globalPreictalWindows, 1);
numWindows = size(globalPreictalWindows, 3);
modelOrder = 10;
ARcoefs_preictal = zeros(numChannels, modelOrder, numWindows); 
ARcoefs_ictal = zeros(numChannels, modelOrder, numWindows); 
ARcoefs_postictal = zeros(numChannels, modelOrder, numWindows); 
for win = 1:numWindows
    for ch = 1:numChannels
        data1 = squeeze(globalPreictalWindows(ch, :, win))';
        model1 = ar(data1, modelOrder, 'burg');
        ARcoefs_preictal(ch, :, win) = model1.A(2:end);
        %
        data2 = squeeze(globalIctalWindows(ch, :, win))';
        model2 = ar(data2, modelOrder, 'burg');
        ARcoefs_ictal(ch, :, win) = model2.A(2:end);
        %
        data3 = squeeze(globalPostictalWindows(ch, :, win))';
        model3 = ar(data3, modelOrder, 'burg');
        ARcoefs_postictal(ch, :, win) = model3.A(2:end);
    end
end
save('Global_AR_Coefs.mat',"ARcoefs_preictal", "ARcoefs_ictal","ARcoefs_postictal")


