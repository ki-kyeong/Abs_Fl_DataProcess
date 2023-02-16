function results = data2spec_v2(name, freq_i, freq_f, freq_step, iter)

det = freq_i : freq_step : freq_f; % 446.xxxx THz
Det = (det-8098.7)*1e2; % MHz, detuning from 7Li D2 F=2
dfsize = size(det,2);

opt = detectImportOptions(name+"0_"+num2str(freq_i,'%.2f')+".csv",'Range','A1:Z2');
expPara = readmatrix(name+"0_"+num2str(freq_i,'%.2f')+".csv",opt);

for i = 1 : size(opt.VariableNames,2)
    results.(opt.VariableNames{i}) = expPara(i);
end

Data = readmatrix(name+"0_"+num2str(freq_i,'%.2f')+".csv");
time = Data(:,1)*1e6; % µs
data_num = size(Data,1);

TotalAbsDatas = zeros(data_num, iter, dfsize); % abs time signal row data - bacground
TotalFlsDatas = zeros(data_num, iter, dfsize); % fls time signal row data
NormAbsDatas = zeros(size(TotalAbsDatas)); % abs time signal normalized with max val after ablation
SumAbsDatas = zeros(dfsize,iter);% summation of each iteration of normalized abs time signal
SumFlsDatas = zeros(dfsize, iter); % summatino of each iteration of fls time signal
NormFlsDatas = zeros(dfsize, iter); % normalized with first fls data
MeanAbsDatas = zeros(dfsize);
MeanFlsDatas = zeros(dfsize);
SteAbsDatas = zeros(dfsize);
SteFlsDatas = zeros(dfsize);

baselinerange = round(96/1.6);
bg = results.DarkVoltage_mV_*1e-3;


for j = 1 : dfsize
    f = det(j)
    for i = 1:iter
        i;
        Data = readmatrix(name+string(i-1)+"_"+num2str(f,'%.2f')+".csv");
        TotalAbsDatas(:,i,j) = Data(:,2)-bg;
        TotalFlsDatas(:,i,j) = Data(:,4);
    end
end


for j = 1: dfsize
    for i = 1: iter
        NormAbsDatas(:,i,j) = 1-(TotalAbsDatas(:,i,j)/mean(TotalAbsDatas(1:baselinerange,i,j)));
            %% abs data는 iter i, det j마다 baseline을 계산해 나눠주었다. 
    end
end

for j = 1 : dfsize
    SumAbsDatas(j,:) = sum(NormAbsDatas(baselinerange:end,:,j));
    SumFlsDatas(j,:) = sum(TotalFlsDatas(baselinerange:end,:,j));
end

NormFlsDatas= SumFlsDatas./SumFlsDatas(1,:); % Fl data는 iter 마다의 첫번째 데이터로 각각 나눠주었다.

MeanAbsDatas = mean(SumAbsDatas,2);
SteAbsDatas = std(SumAbsDatas,0,2)/sqrt(iter);
MeanFlsDatas = mean(NormFlsDatas,2);
SteFlsDatas = std(NormFlsDatas,0,2)/sqrt(iter);


results = struct('t',time, 'freq', det, 'det',Det, 'AT', TotalAbsDatas, ...
    'FT',TotalFlsDatas,'AN',NormAbsDatas,'AS', SumAbsDatas, 'FS', SumFlsDatas,...
    'FN',NormFlsDatas, 'AM',MeanAbsDatas, 'ASte', SteAbsDatas,...
    'FM', MeanFlsDatas, 'FSte', SteFlsDatas);


end
