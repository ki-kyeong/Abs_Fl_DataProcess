function results = readMgFSpecData_v5(name, normmode, plotmode)

clf;
close all;

results.name = name;
savename = split(name,'_');
results.readme=readlines("./data/"+savename(end-1)+"_readme.txt", Encoding="UTF-8");
results.savename = savename(end-1);

opt2 = detectImportOptions(name+"info.csv");
expPara2 = readmatrix(name+"info.csv",opt2);

for i = 1 : size(opt2.VariableNames,2)
    results.(opt2.VariableNames{i}) = expPara2(i);
end


for i = 1 : results.iteration
results.wavemeter_data{:,i} = readmatrix(name+string(i-1)+"_wavemeter_data.csv"); 
end

results.IRcavityfreq = results.wavemeter_data{:,1}(:,1); % THz
results.IRcavitydet = (results.IRcavityfreq-417.147178)*1e6; % MHz
results.UVcavityfreq = 2*results.IRcavityfreq;
results.UVcavitydet = 2*results.IRcavitydet;

dfsize = size(results.IRcavityfreq,1);

results.IRrealfreq = zeros(dfsize,results.iteration);
for i = 1: results.iteration
results.IRrealfreq(:,i) = results.wavemeter_data{i}(:,2);
end

results.IRrealdet = (results.IRrealfreq-417.147178)*1e6; % MHz
results.IRrealdetM = mean(results.IRrealdet,2);
results.IRrealdetSte = std(results.IRrealdet,0,2)/sqrt(results.iteration); 
results.UVrealfreq = 2*results.IRrealfreq;
results.UVrealdet = 2*results.IRrealdet;
results.UVrealdetM = 2*results.IRrealdetM;
results.UVrealdetSte = 2*results.IRrealdetSte;



opt = detectImportOptions(name+"parameter.csv");
expPara = readmatrix(name+"parameter.csv",opt);

for i = 1 : size(opt.VariableNames,2)
    results.(opt.VariableNames{i}) = expPara(i);
end


DataParams = detectImportOptions(name+'0_'+num2str(results.IRcavityfreq(1),'%.6f')+".csv");
Data = readmatrix(name+'0_'+num2str(results.IRcavityfreq(1),'%.6f')+".csv",DataParams);

results.dt = (Data(2,1)-Data(1,1))*1e6; % µs unit

data_num = size(Data,1);

TotalAbsDatas = zeros(data_num, results.repititionPerStep, results.iteration, dfsize); % abs time signal row data - bacground
TotalAbsPFMDatas = zeros(data_num, results.repititionPerStep, results.iteration, dfsize); % abs time signal row data - bacground
% NormAbsDatas = zeros(size(TotalAbsDatas)); % abs time signal normalized with max val after ablation
SumAbsDatas = zeros(dfsize,results.repititionPerStep, results.iteration);% summation of each results.iteration of normalized abs time signal

% TotalFlDatas = zeros(data_num, results.iteration, dfsize); % fls time signal row data
% % NormFlDatas = zeros(size(TotalFlDatas)); % normalized with first fls data
% SumFlDatas = zeros(dfsize, results.iteration); % summatino of each results.iteration of fls time signal



results.baselinerange = round(96/results.dt); % 96 µs에서 ablation 시작 104 µs에서 끝
results.t = (Data(:,1)-Data(results.baselinerange,1))*1e6; % µs, 96 µs을 0초로 설정

wb = waitbar(0, ' Getting started');

for j = 1 : dfsize
    f = results.IRcavityfreq(j);

    waitbar(j/dfsize, wb, results.savename+newline+...
        num2str(f, '%.6f') + " THz"+newline+...
        num2str(j/dfsize*100, '%.1f')+" % done",...
        'WindowStyle','modal');

    for i = 1:results.iteration
        Data = readmatrix(name+string(i-1)+"_"+num2str(f,'%.6f')+".csv",DataParams);
        TotalAbsDatas(:,1:results.repititionPerStep,i,j) = Data(:,2:2:2*results.repititionPerStep)-results.AbsorptionBackgroundVoltage_mV_ *1e-3;
        TotalAbsPFMDatas(:,1:results.repititionPerStep,i,j) = Data(:,2*(2*results.repititionPerStep+1):2:2*(3*results.repititionPerStep))-results.AbsorptionSamplingBackgroundVoltage_mV_*1e-3;
        %         TotalFlDatas(:,i,j) = Data(:,2*(results.repititionPerStep+1):2:2*(2*results.repititionPerStep));
    end
end

close(wb)

results.AT = TotalAbsDatas;
results.ATPFM = TotalAbsPFMDatas;
% results.FT = TotalFlDatas;

switch normmode
    case 'basic'
        results.AN = 1-(TotalAbsDatas/mean(TotalAbsDatas(1:results.baselinerange,1,1)));
        %         results.FN = TotalFlDatas/mean(TotalFlDatas(7000/1.6:end,1,1))-1;

    case 'powernorm'
        for j = 1 : dfsize
            for i = 1 : results.iteration
                for k = 1 : results.repititionPerStep
                nn(:,k,i,j) = 1-(TotalAbsDatas(:,k,i,j)/mean(TotalAbsDatas(1:results.baselinerange,k,i,j)));
                nnsp(:,k,i,j) = 1-(TotalAbsPFMDatas(:,k,i,j)/mean(TotalAbsPFMDatas(1:results.baselinerange,k,i,j)));
                %                 nnn(:,i,j) = (TotalFlDatas(:,i,j)/mean(TotalFlDatas(5000/1.6:end,i,j)))-1;
            end
        end

        results.ANorig = nn;
        results.ANbg = nnsp;
        results.AN = nn-nnsp;
        % results.AN = nn;
        % results.FN = nnn;
        end
    otherwise
        error('try basic of powernorm')

end



for j = 1 : dfsize
    SumAbsDatas(j,:,:) = sum(results.AN(results.baselinerange+round(80/results.dt):end,:,:,j)); % ablation 이후 80 µs 부터
%     SumFlDatas(j,:,:) = sum(results.FN(results.baselinerange+4:end,:,:,j));
end

results.AS = SumAbsDatas;
% results.FS = SumFlDatas;

results.AM = mean(SumAbsDatas(:,:,:),[2 3]); % results.repititionPerStep, results.iteration 전부 평균
results.ASte = std(SumAbsDatas(:,:,:),0,[2 3])/sqrt(results.iteration*results.repititionPerStep);
% results.FM = mean(SumFlDatas,2);
% results.FSte = std(SumFlDatas,0,2)/sqrt(results.iteration);

% % abs spectrum fit으로 Det를 정할때
% results.Det = results.det;
% results.absfit = GF_Li_v2(results, 'single',1); % abs spectrum gaussian fit
% results.Det = results.det-results.absfit.b1; % b1 = F=2 detuning

% z fl spectrum 으로 Det를 정할때
% results.df = results.det(results.FM == max(results.FM));
% results.df= results.det(end);
% results.Det = results.det-results.df;

% results.IRdf= 0;
% 
% results.IRcavityDet = results.IRcavitydet-results.IRdf;

% results.IRcavityfreq = 417.147178+results.IRDet*1e-6; % R_1(1) F=2 전이
% results.IRAbsfreq = 417.20937+results.IRDet*1e-6; % R_1(1) F=2 전이
% results.IRAbsfreq = 417.186005+results.IRDet*1e-6; % R_1(0) F=1 전이
% results.UVAbsfreq = results.IRAbsfreq*2;

% results.v = -results.Det/(446809870)*299792458/cos(pi/4); % for XY beam
% results.v = -results.Det/(446809870)*299792458/cos(pi/4+atan(1.5/10)); % for XYt beam
% results.v = -results.Det/(446809870)*299792458; % for Z beam

% results.maxv = 0.373./(results.t(results.baselinerange:end)*1e-6); % possible maximum velocity
% results.mint = 0.373./results.v*1e6; % possible minimum arrival time, µs

% NormFlDatas2 = results.FN; % normalized with first fls data
% SumFlDatas2 = zeros(dfsize, results.iteration); % summatino of each results.iteration of fls time signal
%
% for j = 1 : dfsize
%    if results.mint(j) <=0
%        NormFlDatas2(:,:,j) = 0;
%    else
%        NormFlDatas2((results.t< results.mint(j)),:,j) = 0;
%    end
%
% end

% results.FN2 = NormFlDatas2;


% for j = 1 : dfsize
%     SumFlDatas2(j,:) = sum(results.FN2(:,:,j));
% end
% results.FS2 = SumFlDatas2;
%
% results.FM2 = mean(results.FS2,2);
% results.FSte2 = std(results.FS2,0,2)/sqrt(results.iteration);

% results.absfit = GF_Li_v2(results, 1); % abs spectrum gaussian fit

if plotmode == true
    plotfigureset_v3(results, results.savename, 'UV')
end


end
