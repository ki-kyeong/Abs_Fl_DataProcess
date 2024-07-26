function results = detcal_v2(results)
% UV wavelength만 읽기 시기시작 24/04/30

%% Q12 detuning
% freqQ11 = 417.170602; % Q1(1) F=1 → F'=1
freqQ121 = 834.294485;
% cavity detuning cal
results.det.cavity.Q = (results.freq.cavity-freqQ121)*1e6; % MHz, from Q12(1)


%IR wavemeter detuning cal and analysis
results.det.wm.raw.Q = (results.freq.wm.raw-freqQ121)*1e6; % MHz
results.det.wm.mean.Q = mean(results.det.wm.raw.Q,[2 3]);
results.det.wm.ste.Q = std(results.det.wm.raw.Q,0,[2 3])/sqrt(results.size.iter); 


%% P1 detuning
% IR cavity detuning cal
freqP11 = 834.294356;
results.det.cavity.P = (results.freq.cavity-freqP11)*1e6; % MHz, from P1(1), F=2

%IR wavemeter detuning cal and analysis
results.det.wm.raw.P = (results.freq.wm.raw-freqP11)*1e6; % MHz
results.det.wm.mean.P = mean(results.det.wm.raw.P,[2 3]);
results.det.wm.ste.P = std(results.det.wm.raw.P,0,[2 3])/sqrt(results.size.iter); 

end

