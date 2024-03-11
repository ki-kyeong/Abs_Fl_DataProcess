function results = readEOMData_v1(name, plotmode)

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

%%

%======================================
% read each iteration wavemeter data 
% first column = cavity freq
% second column = wavemeter freq

for i = 1 : results.size.iter
    results.wavemeter_data{i} = readmatrix(name+string(i-1)+"_wavemeter_data.csv");
end

% first column to cavity freq & det

results.freq.IR.cavity = results.wavemeter_data{:,1}(:,1);

% random array sorting
[results.freq.IR.cavity, results.randidx]=sort(results.freq.IR.cavity);

results.det.IR.cavity = (results.freq.IR.cavity-417.1472425)*1e6; % MHz, from Q12(1)
% results.det.IR.cavity = (results.freq.IR.cavity-417.147178)*1e6; % MHz, from P1(1), F=2

results.freq.UV.cavity = 2*results.freq.IR.cavity;
results.det.UV.cavity = 2*results.det.IR.cavity;

results.size.freq = size(results.freq.IR.cavity,1);

% second column to wavemeter freq & det
results.freq.IR.wm.raw = zeros(results.size.freq,results.size.iter);
for i = 1: results.size.iter
results.freq.IR.wm.raw(:,i) = results.wavemeter_data{:,i}(:,2);
end

results.freq.IR.wm.raw = results.freq.IR.wm.raw(results.randidx,:);

results.freq.IR.wm.mean = mean(results.freq.IR.wm.raw,2);
results.freq.IR.wm.ste = std(results.freq.IR.wm.raw,0,2)/sqrt(results.size.iter);

results.det.IR.wm.raw = (results.freq.IR.wm.raw-417.1472425)*1e6; % MHz
% results.det.IR.wm.raw = (results.freq.IR.wm.raw-417.147178)*1e6; % MHz
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

TotalEOMDatas = zeros(results.size.time, results.size.rep, results.size.iter, results.size.freq); % abs time signal row data - bacground
TotalPBSRDatas = zeros(results.size.time, results.size.rep, results.size.iter, results.size.freq); % abs time signal row data - bacground
TotalPBSTDatas = zeros(results.size.time, results.size.rep, results.size.iter, results.size.freq); % abs time signal row data - bacground

results.t = Data(:,1)*1e6; 

wb = waitbar(0, ' Getting started');

for j = 1 : results.size.freq
    f = results.freq.IR.cavity(j);

    waitbar(j/results.size.freq, wb, results.savename+newline+...
        num2str(f, '%.6f') + " THz"+newline+...
        num2str(j/results.size.freq*100, '%.1f')+" % done",...
        'WindowStyle','modal');

    for i = 1:results.size.iter
        Data = readmatrix(name+string(i-1)+"_"+num2str(f,'%.6f')+".csv",DataParams);
        TotalEOMDatas(:,1:results.size.rep,i,j) = Data(:,2:2:2*results.size.rep)-expPara(7)*1e-3;
        TotalPBSRDatas(:,1:results.size.rep,i,j) = Data(:,2*(2*results.size.rep+1):2:2*(3*results.size.rep));
        TotalPBSTDatas(:,1:results.size.rep,i,j) = Data(:,2*(results.size.rep+1):2:2*(2*results.size.rep));
    end
end

close(wb)

results.data.raw = TotalEOMDatas;
results.data.PBSR = TotalPBSRDatas;
results.data.PBST = TotalPBSTDatas;

results.data.mean = squeeze(mean(TotalEOMDatas,2));

results.PBSR.sum = sum(TotalPBSRDatas,1);
results.PBST.sum = sum(TotalPBSTDatas,1);


    for j = 1: results.size.iter 
        for k = 1: results.size.freq
        % [~,idx] = max(results.data.mean(:,j,k));
        results.fit(j,k) = EOMfit_v1(results.t*0.3965, results.data.mean(:,j,k),115,5,[0.1,7,2500*0.3965,1,0]);
        results.beta(j,k) = results.fit(j,k).beta;
        end
    end

% for i = 1:results.size.rep
%     for j = 1: results.size.iter 
%         for k = 1: results.size.freq
%         [~,idx] = max(results.data.raw(:,i,j,k));
%         results.fit(i,j,k) = EOMfit_v1(results.t*0.3243, results.data.raw(:,i,j,k),115,2,[0.1,7,results.t(idx)*0.3243,1,0]);
%         results.beta(i,j,k) = results.fit(i,j,k).beta;
%         end
%     end
% end

if plotmode == true
    plotfigureset_v4(results,'UV')
end


end
