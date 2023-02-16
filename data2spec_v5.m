function results = data2spec_v5(name, freq_i, freq_f, freq_step, iter, normmode, plotmode)

results.name = name;
savename = split(name,'_');
savename = savename(end-1);
det = freq_i : freq_step : freq_f; % 446.xxxx THz
results.det = ((det-8098.7)*1e2).'; % MHz, detuning from 7Li D2 F=2
results.absfreq = 446+det*1e-4;
dfsize = size(det,2);

opt = detectImportOptions(name+"0_"+num2str(results.absfreq(1),'%.6f')+".csv",'Range','A1:Z2');
expPara = readmatrix(name+"0_"+num2str(results.absfreq(1),'%.6f')+".csv",opt);

for i = 1 : size(opt.VariableNames,2)
    results.(opt.VariableNames{i}) = expPara(i);
end

Data = readmatrix(name+"0_"+num2str(results.absfreq(1),'%.6f')+".csv");

data_num = size(Data,1);

TotalAbsDatas = zeros(data_num, iter, dfsize); % abs time signal row data - bacground
TotalFlDatas = zeros(data_num, iter, dfsize); % fls time signal row data
NormAbsDatas = zeros(size(TotalAbsDatas)); % abs time signal normalized with max val after ablation
NormFlDatas = zeros(size(TotalFlDatas)); % normalized with first fls data
SumAbsDatas = zeros(dfsize,iter);% summation of each iteration of normalized abs time signal
SumFlDatas = zeros(dfsize, iter); % summatino of each iteration of fls time signal


results.baselinerange = round(96/1.6); % 96 µs에서 ablation 시작 104 µs에서 끝
results.t = (Data(:,1)-Data(results.baselinerange,1))*1e6; % µs, 96 µs을 0초로 설정

bg = results.DarkVoltage_mV_*1e-3;


for j = 1 : dfsize
    f = results.absfreq(j);
    det(j)
    for i = 1:iter
        i;
        Data = readmatrix(name+string(i-1)+"_"+num2str(f,'%.6f')+".csv");
        TotalAbsDatas(:,i,j) = Data(:,2)-bg;
        TotalFlDatas(:,i,j) = Data(:,4);
    end
end

results.AT = TotalAbsDatas;
results.FT = TotalFlDatas;

switch normmode
    case 'basic'
        
        results.AN = 1-(TotalAbsDatas/mean(TotalAbsDatas(1:results.baselinerange,1,1)));
        results.FN = TotalFlDatas/mean(TotalFlDatas(7000/1.6:end,1,1))-1;
        
    case 'powernorm'
        for j = 1 : dfsize
            for i = 1 : iter
                nn(:,i,j) = 1-(TotalAbsDatas(:,i,j)/mean(TotalAbsDatas(1:results.baselinerange,i,j)));
                nnn(:,i,j) = (TotalFlDatas(:,i,j)/mean(TotalFlDatas(7000/1.6:end,i,j)))-1;
            end
        end
        results.AN = nn;
        results.FN = nnn;
        
    otherwise
        error('try basic of powernorm')
end



for j = 1 : dfsize
    SumAbsDatas(j,:) = sum(results.AN(results.baselinerange+4:end,:,j));
    SumFlDatas(j,:) = sum(results.FN(results.baselinerange+4:end,:,j));
end

results.AS = SumAbsDatas;
results.FS = SumFlDatas;


results.AM = mean(SumAbsDatas,2);
results.ASte = std(SumAbsDatas,0,2)/sqrt(iter);
results.FM = mean(SumFlDatas,2);
results.FSte = std(SumFlDatas,0,2)/sqrt(iter);

results.absfit = GF_Li(results.det, results.AM); % abs spectrum gaussian fit

results.Det = results.det-results.absfit.b1; % b1 = F=2 detuning
results.Absfreq = 446.80987+results.Det*1e-4;

results.v = -results.Det/(446809870)*299792458/cos(pi/4+atan(1.5/10)); % for XY beam

if plotmode == true
    plotfigureset_v1(results, savename)
end


end
