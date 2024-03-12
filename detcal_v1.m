function results = detcal_v1(results)
%% Q12 detuning
% freqQ11 = 417.170602; % Q1(1) F=1 → F'=1
freqQ121 = 417.1472425;
% IR cavity detuning cal
results.det.IR.cavity.Q = (results.freq.IR.cavity-freqQ121)*1e6; % MHz, from Q12(1)
% UV cavity detuning cal
results.det.UV.cavity.Q = 2*results.det.IR.cavity.Q;

%IR wavemeter detuning cal and analysis
results.det.IR.wm.raw.Q = (results.freq.IR.wm.raw-freqQ121)*1e6; % MHz
results.det.IR.wm.mean.Q = mean(results.det.IR.wm.raw.Q,2);
results.det.IR.wm.ste.Q = std(results.det.IR.wm.raw.Q,0,2)/sqrt(results.size.iter); 

%UV wavemeter detuning cal and analysis
results.det.UV.wm.raw.Q = 2*results.det.IR.wm.raw.Q;
results.det.UV.wm.mean.Q = 2*results.det.IR.wm.mean.Q;
results.det.UV.wm.ste.Q = std(results.det.UV.wm.raw.Q,0,2)/sqrt(results.size.iter); 

%% P1 detuning
% IR cavity detuning cal
freqP11 = 417.147178;
results.det.IR.cavity.P = (results.freq.IR.cavity-freqP11)*1e6; % MHz, from P1(1), F=2
% UV cavity detuning cal
results.det.UV.cavity.P = 2*results.det.IR.cavity.P;

%IR wavemeter detuning cal and analysis
results.det.IR.wm.raw.P = (results.freq.IR.wm.raw-freqP11)*1e6; % MHz
results.det.IR.wm.mean.P = mean(results.det.IR.wm.raw.P,2);
results.det.IR.wm.ste.P = std(results.det.IR.wm.raw.P,0,2)/sqrt(results.size.iter); 

%UV wavemeter detuning cal and analysis
results.det.UV.wm.raw.P = 2*results.det.IR.wm.raw.P;
results.det.UV.wm.mean.P = 2*results.det.IR.wm.mean.P;
results.det.UV.wm.ste.P = std(results.det.UV.wm.raw.P,0,2)/sqrt(results.size.iter); 


end

