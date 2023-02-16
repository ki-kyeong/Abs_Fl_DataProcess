function results = data2spec_v4(name, freq_i, freq_f, freq_step, iter, savename, normmode)

results.name = name;
det = freq_i : freq_step : freq_f; % 446.xxxx THz
results.det = ((det-8098.7)*1e2).'; % MHz, detuning from 7Li D2 F=2
results.v = -results.det/(446809870)*2997924581/cos(pi/180*45); % for XY beam
results.absfreq = 446+det*1e-4;
dfsize = size(det,2);

opt = detectImportOptions(name+"0_"+num2str(det(1),'%.2f')+".csv",'Range','A1:Z2');
expPara = readmatrix(name+"0_"+num2str(det(1),'%.2f')+".csv",opt);

for i = 1 : size(opt.VariableNames,2)
    results.(opt.VariableNames{i}) = expPara(i);
end

Data = readmatrix(name+"0_"+num2str(det(1),'%.2f')+".csv");
time = Data(:,1)*1e6; % µs
results.t = time;
data_num = size(Data,1);

TotalAbsDatas = zeros(data_num, iter, dfsize); % abs time signal row data - bacground
TotalFlDatas = zeros(data_num, iter, dfsize); % fls time signal row data
NormAbsDatas = zeros(size(TotalAbsDatas)); % abs time signal normalized with max val after ablation
NormFlDatas = zeros(size(TotalFlDatas)); % normalized with first fls data
SumAbsDatas = zeros(dfsize,iter);% summation of each iteration of normalized abs time signal
SumFlDatas = zeros(dfsize, iter); % summatino of each iteration of fls time signal


results.baselinerange = round(96/1.6);
bg = results.DarkVoltage_mV_*1e-3;


for j = 1 : dfsize
    f = det(j);
    det(j)
    for i = 1:iter
        i;
        Data = readmatrix(name+string(i-1)+"_"+num2str(f,'%.2f')+".csv");
        TotalAbsDatas(:,i,j) = Data(:,2)-bg;
        TotalFlDatas(:,i,j) = Data(:,4);
    end
end

results.AT = TotalAbsDatas;
results.FT = TotalFlDatas;

switch normmode
    case 'basic'
        
        results.AN = 1-(TotalAbsDatas/mean(TotalAbsDatas(1:results.baselinerange,1,1)));
        results.FN = TotalFlDatas/mean(TotalFlDatas(1:results.baselinerange,1,1))-1;
        
    case 'powernorm'
        for j = 1 : dfsize
            for i = 1 : iter
                nn(:,i,j) = 1-(TotalAbsDatas(:,i,j)/mean(TotalAbsDatas(1:results.baselinerange,i,j)));
                nnn(:,i,j) = (TotalFlDatas(:,i,j)/mean(TotalFlDatas(1:results.baselinerange,i,j)))-1;
                                
%                 nnn(:,i,j) = (TotalFlDatas(:,i,j)/mean(TotalFlDatas(7000/1.6:end,i,j)))-1;
            end
        end
        results.AN = nn;
        results.FN = nnn;
        
    otherwise
        error('try basic of powernorm')
end



for j = 1 : dfsize
    SumAbsDatas(j,:) = sum(results.AN(results.baselinerange+4:end,:,j));
    SumFlDatas(j,:) = sum(results.FN(results.baselinerange+4:6000/1.6,:,j));
end

results.AS = SumAbsDatas;
results.FS = SumFlDatas;


results.AM = mean(SumAbsDatas,2);
results.ASte = std(SumAbsDatas,0,2)/sqrt(iter);
results.FM = mean(SumFlDatas,2);
results.FSte = std(SumFlDatas,0,2)/sqrt(iter);


clf;
figure('Name',savename+" abs spec")
plotspec_v1(gca, results, 'abs');
saveas(gcf, './K_results/'+savename+'_abs_spectrum.png');

figure('Name',savename+" fl spec")
plotspec_v1(gca, results, 'fl');
saveas(gcf, './K_results/'+savename+'_fl_spectrum.png');

fig=figure('Name',savename+" absfl spec")
fig.Position = [120 130 2*560 1*420]
tiledlayout(1,2)
plotspec_v1(nexttile, results, 'abs')
plotspec_v1(nexttile, results, 'fl')
saveas(gcf, './K_results/'+savename+'_absfl_spectrum.png');

figure('Name',savename+" abs 2D")
plot2Dtimedet_v1(gca, results, 'abs');
saveas(gcf, './K_results/'+savename+'_abs_2D.png');
%
figure('Name',savename+" fl 2D")
plot2Dtimedet_v1(gca, results, 'fl');
saveas(gcf, './K_results/'+savename+'_fl_2D.png');


end
