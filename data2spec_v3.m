function results = data2spec_v3(name, freq_i, freq_f, freq_step, iter)

results.name = name;
det = freq_i : freq_step : freq_f; % 446.xxxx THz
results.det = ((det-8098.7)*1e2).'; % MHz, detuning from 7Li D2 F=2
dfsize = size(det,2);

opt = detectImportOptions(name+"0_"+num2str(freq_i,'%.2f')+".csv",'Range','A1:Z2');
expPara = readmatrix(name+"0_"+num2str(freq_i,'%.2f')+".csv",opt);

for i = 1 : size(opt.VariableNames,2)
    results.(opt.VariableNames{i}) = expPara(i);
end

Data = readmatrix(name+"0_"+num2str(freq_i,'%.2f')+".csv");
time = Data(:,1)*1e6; % µs
results.t = time;
data_num = size(Data,1);

TotalAbsDatas = zeros(data_num, iter, dfsize); % abs time signal row data - bacground
TotalFlDatas = zeros(data_num, iter, dfsize); % fls time signal row data
NormAbsDatas = zeros(size(TotalAbsDatas)); % abs time signal normalized with max val after ablation
NormFlDatas = zeros(size(TotalFlDatas)); % normalized with first fls data
SumAbsDatas = zeros(dfsize,iter);% summation of each iteration of normalized abs time signal
SumFlDatas = zeros(dfsize, iter); % summatino of each iteration of fls time signal


baselinerange = round(96/1.6);
bg = results.DarkVoltage_mV_*1e-3;


for j = 1 : dfsize
    f = det(j)
    for i = 1:iter
        i;
        Data = readmatrix(name+string(i-1)+"_"+num2str(f,'%.2f')+".csv");
        TotalAbsDatas(:,i,j) = Data(:,2)-bg;
        TotalFlDatas(:,i,j) = Data(:,4);
    end
end

results.AT = TotalAbsDatas;
results.FT = TotalFlDatas;


results.AN = 1-(TotalAbsDatas/mean(TotalAbsDatas(1:baselinerange,1,1)));
results.FN = TotalFlDatas/TotalFlDatas(baselinerange,1,1);

for j = 1 : dfsize
    SumAbsDatas(j,:) = sum(results.AN(baselinerange+4:end,:,j));
    SumFlDatas(j,:) = sum(results.FN(baselinerange:end,:,j));
end

results.AS = SumAbsDatas;
results.FS = SumFlDatas;


results.AM = mean(SumAbsDatas,2);
results.ASte = std(SumAbsDatas,0,2)/sqrt(iter);
results.FM = mean(SumFlDatas,2);
results.FSte = std(SumFlDatas,0,2)/sqrt(iter);


end
