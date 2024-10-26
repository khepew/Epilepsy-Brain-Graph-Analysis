clc;
clear;
load("Global_Phase_Windows.mat");

% freqBands = [0.5 2; 2 4; 4 6; 6 8; 8 10; 10 12; 12 14; 14 16; ...
%              16 18; 18 20; 20 22; 22 24; 24 26; 26 28; ...
%              28 30; 30 32; 32 34; 34 36; 36 38; 38 40]; 
% 
freqBands = [0.5 2; 
             2 4; 
             4 6; 
             6 8; 
             8 10; 
             10 12; 
             12 14; 
             14 16; 
             16 18; 
             18 20; 
             20 22; 
             22 24; 
             24 26; 
             26 28; 
             28 30; 
             30 32; 
             32 34; 
             34 36; 
             36 38; 
             38 40; 
             40 42; 
             42 44; 
             44 46; 
             46 48; 
             48 50; 
             50 52; 
             52 54; 
             54 56; 
             56 58; 
             58 60; 
             60 62; 
             62 64; 
             64 66; 
             66 68; 
             68 70; 
             70 72; 
             72 74; 
             74 76; 
             76 78; 
             78 80; 
             80 82; 
             82 84; 
             84 86; 
             86 88; 
             88 90; 
             90 92; 
             92 94; 
             94 96; 
             96 98; 
             98 100; 
             100 102; 
             102 104; 
             104 106; 
             106 108; 
             108 110; 
             110 112; 
             112 114; 
             114 116; 
             116 118; 
             118 120; 
             120 122; 
             122 124; 
             124 126; 
             126 128];
% 
%--------------------
% freqBands = [128 130]; 
numBands = size(freqBands, 1);
numChannels = size(globalPreictalWindows, 1);
numWindows = size(globalPreictalWindows, 3);
powerPreictal = zeros(numChannels, numBands, numWindows);
powerIctal = zeros(numChannels, numBands, numWindows);
powerPostictal = zeros(numChannels, numBands, numWindows);
fs = 256;
for win = 1:numWindows
    for ch = 1:numChannels
        % Preictal data
        data1 = squeeze(globalPreictalWindows(ch, :, win))';
        [pxx1, f1] = pwelch(data1, [], [], [], fs); 
        powerPreictal(ch, :, win) = arrayfun(@(i) trapz(f1(f1 >= freqBands(i, 1) & f1 < freqBands(i, 2)), ...
            pxx1(f1 >= freqBands(i, 1) & f1 < freqBands(i, 2))), 1:numBands);

        % Ictal data
        data2 = squeeze(globalIctalWindows(ch, :, win))';
        [pxx2, f2] = pwelch(data2, [], [], [], fs);
        powerIctal(ch, :, win) = arrayfun(@(i) trapz(f2(f2 >= freqBands(i, 1) & f2 < freqBands(i, 2)), ...
            pxx2(f2 >= freqBands(i, 1) & f2 < freqBands(i, 2))), 1:numBands);
        
        % Postictal data
        data3 = squeeze(globalPostictalWindows(ch, :, win))';
        [pxx3, f3] = pwelch(data3, [], [], [], fs);
        powerPostictal(ch, :, win) = arrayfun(@(i) trapz(f3(f3 >= freqBands(i, 1) & f3 < freqBands(i, 2)), ...
            pxx3(f3 >= freqBands(i, 1) & f3 < freqBands(i, 2))), 1:numBands);
    end
end

save('Global_BandPower.mat',"powerPreictal", "powerIctal","powerPostictal")

