function results = readMgFPMTData_v1(name)

clf;
close all;

%===========================
results.name = name;
savename = split(name,'_');
results.savename = savename(end-1);


%===========================
results.readme=readlines(results.name+"readme.txt", Encoding="UTF-8");

%===========================
% read run informations

opt2 = detectImportOptions(name+"info.csv");
expPara2 = readmatrix(name+"info.csv",opt2);

for i = 1 : size(opt2.VariableNames,2)
    results.(opt2.VariableNames{i}) = expPara2(i);
end

results.size.iter = results.iteration;

%======================================
% read each iteration wavemeter data 
% first column = cavity freq
% second column = wavemeter freq

for i = 1 : results.size.iter
    results.wavemeter_data{i} = readmatrix(name+string(i-1)+"_wavemeter_data.csv");
end

% first column to cavity freq & det

results.freq.IR.cavity = results.wavemeter_data{:,1}(:,1);
results.det.IR.cavity = (results.freq.IR.cavity-417.1472425)*1e6; % MHz, from Q12(1)
results.freq.UV.cavity = 2*results.freq.IR.cavity;
results.det.UV.cavity = 2*results.det.IR.cavity;

results.size.freq = size(results.freq.IR.cavity,1);

% second column to wavemeter freq & det
results.freq.IR.wm.raw = zeros(results.size.freq,results.size.iter);
for i = 1: results.size.iter
results.freq.IR.wm.raw(:,i) = results.wavemeter_data{:,i}(:,2);
end

results.freq.IR.wm.mean = mean(results.freq.IR.wm.raw,2);
results.freq.IR.wm.ste = std(results.freq.IR.wm.raw,0,2)/sqrt(results.size.iter);

results.det.IR.wm.raw = (results.freq.IR.wm.raw-417.1472425)*1e6; % MHz
results.det.IR.wm.mean = mean(results.det.IR.wm.raw,2);
results.det.IR.wm.ste = std(results.det.IR.wm.raw,0,2)/sqrt(results.size.iter); 

results.freq.UV.wm.raw = 2*results.freq.IR.wm.raw;
results.freq.UV.wm.mean = 2*results.freq.IR.wm.mean;
results.freq.UV.wm.ste = std(results.freq.UV.wm.raw,0,2)/sqrt(results.size.iter); 

results.det.UV.wm.raw = 2*results.det.IR.wm.raw;
results.det.UV.wm.mean = 2*results.det.IR.wm.mean;
results.det.UV.wm.ste = std(results.det.UV.wm.raw,0,2)/sqrt(results.size.iter); 


%======================================
% read parameter data

opt = detectImportOptions(name+"parameter.csv");
expPara = readmatrix(name+"parameter.csv",opt);

for i = 1 : size(opt.VariableNames,2)
    results.(opt.VariableNames{i}) = expPara(i);
end

results.size.rep = results.repititionPerStep;

%======================================
% read data for initializing

DataParams = detectImportOptions(name+'0_'+num2str(results.freq.IR.cavity(1),'%.6f')+".csv");
Data = readmatrix(name+'0_'+num2str(results.freq.IR.cavity(1),'%.6f')+".csv",DataParams);

results.dt = (Data(2,1)-Data(1,1))*1e6; % µs unit

results.size.time = size(Data,1);

TotalAbsDatas = zeros(results.size.time, results.size.rep, results.size.iter, results.size.freq); % abs time signal row data - bacground
TotalAbsPFMDatas = zeros(results.size.time, results.size.rep, results.size.iter, results.size.freq); % abs time signal row data - bacground
% NormAbsDatas = zeros(size(TotalAbsDatas)); % abs time signal normalized with max val after ablation
% SumAbsDatas = zeros(results.size.freq,results.size.rep, results.size.iter);% summation of each results.iteration of normalized abs time signal
TotalFlDatas = zeros(results.size.time,results.size.rep, results.size.iter, results.size.freq); % fls time signal row data
% NormFlDatas = zeros(size(TotalFlDatas)); % normalized with first fls data
% SumFlDatas = zeros(results.size.freq,results.size.rep, results.size.iter); % summatino of each results.iteration of fls time signal



results.baselineidx = round(96/results.dt); % 96 µs에서 ablation 시작 140 µs에서 끝
results.baselinerange = [1 results.baselineidx];
results.t = (Data(:,1)-Data(results.baselineidx,1))*1e6; % µs, 98 µs을 0초로 설정

wb = waitbar(0, ' Getting started');

for j = 1 : results.size.freq
    f = results.freq.IR.cavity(j);

    waitbar(j/results.size.freq, wb, results.savename+newline+...
        num2str(f, '%.6f') + " THz"+newline+...
        num2str(j/results.size.freq*100, '%.1f')+" % done",...
        'WindowStyle','modal');

    for i = 1:results.size.iter
        Data = readmatrix(name+string(i-1)+"_"+num2str(f,'%.6f')+".csv",DataParams);
        TotalAbsDatas(:,1:results.size.rep,i,j) = Data(:,2:2:2*results.size.rep)-results.AbsBgVoltage_mV_ *1e-3;
        TotalAbsPFMDatas(:,1:results.size.rep,i,j) = Data(:,2*(2*results.size.rep+1):2:2*(3*results.size.rep))-results.AbsBgVoltage_mV_*1e-3;
        TotalFlDatas(:,1:results.size.rep,i,j) = Data(:,2*(results.size.rep+1):2:2*(2*results.size.rep));
    end
end

close(wb)

results.abs.raw = TotalAbsDatas;
results.abs.pfm = TotalAbsPFMDatas;
results.fl.raw = TotalFlDatas;



%% rep들을 미리 평균낸것도 계산해주자.


results.abs.tt = reshape(mean(results.abs.raw, [2 3]),results.size.time, results.size.freq);
results.fl.tt = reshape(mean(results.fl.raw, [2 3]), results.size.time, results.size.freq);

results.abs.power.rep = mean(results.abs.raw(results.baselineidx:end,:,:,:));
results.abs.power.mean = mean(results.abs.power.rep,2);
results.abs.power.ste = std(results.abs.power.rep,0,2)/sqrt(results.size.rep*results.size.iter);

results.fl.power.rep = mean(results.fl.raw(results.baselineidx:end,:,:,:));
results.fl.power.mean = mean(results.fl.power.rep,2);
results.fl.power.ste = std(results.fl.power.rep,0,2)/sqrt(results.size.rep*results.size.iter);

results.abs.rms.rep = rms(results.abs.raw-results.abs.power.rep);
results.abs.rms.mean = mean(results.abs.rms.rep);
results.abs.rms.ste = std(results.abs.rms.rep,0)/sqrt(results.size.rep*results.size.iter);

results.fl.rms.rep = rms(results.fl.raw-results.fl.power.rep);
results.fl.rms.mean = mean(results.fl.rms.rep);
results.fl.rms.ste = std(results.fl.rms.rep,0)/sqrt(results.size.rep*results.size.iter);

end
