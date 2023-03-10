function results = readLiSpecData_v1(name, freq_i, freq_f, freq_step, iter, normmode, plotmode)

clf;
close all;

results.name = name;
savename = split(name,'_');
results.savename = savename(end-1);
det = freq_i : freq_step*1e-3 : freq_f; % GHz, 446.XXX
results.det = ((det-809.87)*1e3).'; % MHz, detuning from 7Li D2 F=2
results.absfreq = 446+det*1e-3; % THz
dfsize = size(det,2);

opt = detectImportOptions(name+"parameter.csv");
expPara = readmatrix(name+"parameter.csv",opt);

for i = 1 : size(opt.VariableNames,2)
    results.(opt.VariableNames{i}) = expPara(i);
end

DataParams = detectImportOptions(name+'0_'+num2str(results.absfreq(1),'%.6f')+".csv");
Data = readmatrix(name+'0_'+num2str(results.absfreq(1),'%.6f')+".csv",DataParams);

results.dt = (Data(2,1)- Data(1,1))*1e6;

data_num = size(Data,1);

TotalAbsDatas = zeros(data_num, iter, dfsize); % abs time signal row data - bacground
TotalFlDatas = zeros(data_num, iter, dfsize); % fls time signal row data
% NormAbsDatas = zeros(size(TotalAbsDatas)); % abs time signal normalized with max val after ablation
% NormFlDatas = zeros(size(TotalFlDatas)); % normalized with first fls data
SumAbsDatas = zeros(dfsize,iter);% summation of each iteration of normalized abs time signal
SumFlDatas = zeros(dfsize, iter); % summatino of each iteration of fls time signal

results.baselinerange = round(96/results.dt); % 96 µs에서 ablation 시작 104 µs에서 끝
results.t = (Data(:,1)-Data(results.baselinerange,1))*1e6; % µs, 96 µs을 0초로 설정

bg = results.AbsorptionPower_final__uW_*1e-3;

wb = waitbar(0, ' Getting started');

for j = 1 : dfsize
    f = results.absfreq(j);

waitbar(j/dfsize, wb, results.savename+newline+...
        num2str(f, '%.6f') + " THz"+newline+...
        num2str(j/dfsize*100, '%.1f')+" % done",...
        'WindowStyle','modal');
    
    for i = 1:iter
        Data = readmatrix(name+string(i-1)+"_"+num2str(f,'%.6f')+".csv", DataParams);
        TotalAbsDatas(:,i,j) = Data(:,2)-bg;
        TotalFlDatas(:,i,j) = Data(:,4);
    end
end

close(wb)

results.AT = TotalAbsDatas;
results.FT = TotalFlDatas;

switch normmode
    case 'basic'
        results.AN = 1-(TotalAbsDatas/mean(TotalAbsDatas(1:results.baselinerange,1,1)));
        results.FN = TotalFlDatas/mean(TotalFlDatas(7000/results.dt:end,1,1))-1;
        
    case 'powernorm'
        for j = 1 : dfsize
            for i = 1 : iter
                nn(:,i,j) = 1-(TotalAbsDatas(:,i,j)/mean(TotalAbsDatas(1:results.baselinerange,i,j)));
                nnn(:,i,j) = (TotalFlDatas(:,i,j)/mean(TotalFlDatas(5000/results.dt:end,i,j)))-1;
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

% abs spectrum fit으로 Det를 정할때
results.Det = results.det;
results.absfit = GF_Li_v2(results, 'double',1); % abs spectrum gaussian fit
results.Det = results.det-results.absfit.b1; % b1 = F=2 detuning

% z fl spectrum 으로 Det를 정할때
% results.df = results.det(results.FM == max(results.FM));
% results.df= results.det(end);
% results.Det = results.det-results.df;

results.Absfreq = 446.80987+results.Det*1e-6; %THz

% results.v = -results.Det/(446809870)*299792458/cos(pi/4); % for XY beam
% results.v = -results.Det/(446809870)*299792458/cos(pi/4+atan(1.5/10)); % for XYt beam
% % results.v = -results.Det/(446809870)*299792458; % for Z beam
% 
% results.maxv = 0.373./(results.t(results.baselinerange:end)*1e-6); % possible maximum velocity
% results.mint = 0.373./results.v*1e6; % possible minimum arrival time, µs
% 
% NormFlDatas2 = results.FN; % normalized with first fls data
% SumFlDatas2 = zeros(dfsize, iter); % summatino of each iteration of fls time signal

% for j = 1 : dfsize
%    if results.mint(j) <=0
%        NormFlDatas2(:,:,j) = 0;
%    else
%        NormFlDatas2((results.t< results.mint(j)),:,j) = 0;
%    end
%     
% end
% 
% results.FN2 = NormFlDatas2;


% for j = 1 : dfsize
%     SumFlDatas2(j,:) = sum(results.FN2(:,:,j));
% end
% results.FS2 = SumFlDatas2;
% 
% results.FM2 = mean(results.FS2,2);
% results.FSte2 = std(results.FS2,0,2)/sqrt(iter);    

% results.absfit = GF_Li_v2(results, 1); % abs spectrum gaussian fit

if plotmode == true
    plotfigureset_v2(results, results.savename)
end

end
