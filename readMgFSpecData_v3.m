function results = readMgFSpecData_v3(name, freq_i, freq_f, freq_step,...
    rep, iter, normmode, plotmode)

% update log 23/02/14
% v2에서 repetition data를 읽어오게끔 변경. Mg signal이 작아서 한번에 rep 번의 pulse때려야
% signal이 나오기 때문.. 현재는 rep번째 pulse의 data만 볼 것임.

clf;
close all;

results.name = name;
savename = split(name,'_');
results.savename = savename(end-1);
det = freq_i : freq_step : freq_f; % 417.xxxx THz
% results.IRdet = ((det-2093.71)*1e2).'; % MHz, detuning from 24MgF R_1(1) F=2 전이
results.IRdet = ((det-1860.055)*1e2).'; % MHz, detuning from 24MgF R_1(0) F=1 전이
results.UVdet = results.IRdet*2;
results.IRabsfreq = 417+det*1e-4;
results.UVabsfreq = results.IRabsfreq*2;
dfsize = size(det,2);
results.rep = rep;

opt = detectImportOptions(name+"parameter.csv");
expPara = readmatrix(name+"parameter.csv",opt);


for i = 1 : size(opt.VariableNames,2)
    results.(opt.VariableNames{i}) = expPara(i);
end

DataParams = detectImportOptions(name+'0_'+num2str(results.IRabsfreq(1),'%.6f')+".csv");
Data = readmatrix(name+'0_'+num2str(results.IRabsfreq(1),'%.6f')+".csv",DataParams);

results.dt = (Data(2,1)-Data(1,1))*1e6; % µs unit


data_num = size(Data,1);

TotalAbsDatas = zeros(data_num, iter, dfsize); % abs time signal row data - bacground
TotalAbsBgDatas = zeros(data_num, iter, dfsize); % abs time signal row data - bacground
% NormAbsDatas = zeros(size(TotalAbsDatas)); % abs time signal normalized with max val after ablation
SumAbsDatas = zeros(dfsize,iter);% summation of each iteration of normalized abs time signal

% TotalFlDatas = zeros(data_num, iter, dfsize); % fls time signal row data
% % NormFlDatas = zeros(size(TotalFlDatas)); % normalized with first fls data
% SumFlDatas = zeros(dfsize, iter); % summatino of each iteration of fls time signal


results.baselinerange = round(96/results.dt); % 96 µs에서 ablation 시작 104 µs에서 끝
results.t = (Data(:,1)-Data(results.baselinerange,1))*1e6; % µs, 96 µs을 0초로 설정

bg = results.AbsorptionBackgroundVoltage_mV_ *1e-3;
bg2 = results.AbsorptionSamplingBackgroundVoltage_mV_*1e-3;
% bg = 0*1e-3;

wb = waitbar(0, ' Getting started');

for j = 1 : dfsize
    f = results.IRabsfreq(j);

    waitbar(j/dfsize, wb, results.savename+newline+...
        num2str(f, '%.6f') + " THz"+newline+...
        num2str(j/dfsize*100, '%.1f')+" % done",...
        'WindowStyle','modal');

    for i = 1:iter
        %         i
        Data = readmatrix(name+string(i-1)+"_"+num2str(f,'%.6f')+".csv",DataParams);

        % Data = readmatrix(name+string(i-1)+"_"+num2str(f,'%.2f')+".csv");
        TotalAbsDatas(:,i,j) = Data(:,2*rep)-bg;
        TotalAbsBgDatas(:,i,j) = Data(:,2*rep+4*rep)-bg2;
        %         TotalFlDatas(:,i,j) = Data(:,2*rep+2*rep);
    end
end
close(wb)

results.AT = TotalAbsDatas;
results.ATbg = TotalAbsBgDatas;
% results.FT = TotalFlDatas;

switch normmode
    case 'basic'

        results.AN = 1-(TotalAbsDatas/mean(TotalAbsDatas(1:results.baselinerange,1,1)));
        %         results.FN = TotalFlDatas/mean(TotalFlDatas(7000/1.6:end,1,1))-1;

    case 'powernorm'
        for j = 1 : dfsize
            for i = 1 : iter
                nn(:,i,j) = 1-(TotalAbsDatas(:,i,j)/mean(TotalAbsDatas(1:results.baselinerange,i,j)));
                nnbg(:,i,j) = 1-(TotalAbsBgDatas(:,i,j)/mean(TotalAbsBgDatas(1:results.baselinerange,i,j)));

                %                 nnn(:,i,j) = (TotalFlDatas(:,i,j)/mean(TotalFlDatas(5000/1.6:end,i,j)))-1;
            end
        end
%         results.ANorig = nn;
        results.ANbg = nnbg;
%                 results.AN = nn-nnbg;
results.AN = nn;
        %         results.FN = nnn;

    otherwise
        error('try basic of powernorm')
end



for j = 1 : dfsize
    SumAbsDatas(j,:) = sum(results.AN(results.baselinerange+round(80/results.dt):end,:,j)); % ablation 이후 10 µs 부터
    %     SumFlDatas(j,:) = sum(results.FN(results.baselinerange+4:end,:,j));
end

results.AS = SumAbsDatas;
% results.FS = SumFlDatas;


results.AM = mean(SumAbsDatas,2);
results.ASte = std(SumAbsDatas,0,2)/sqrt(iter);
% results.FM = mean(SumFlDatas,2);
% results.FSte = std(SumFlDatas,0,2)/sqrt(iter);

% % abs spectrum fit으로 Det를 정할때
% results.Det = results.det;
% results.absfit = GF_Li_v2(results, 'single',1); % abs spectrum gaussian fit
% results.Det = results.det-results.absfit.b1; % b1 = F=2 detuning

% z fl spectrum 으로 Det를 정할때
% results.df = results.det(results.FM == max(results.FM));
% results.df= results.det(end);
% results.Det = results.det-results.df;

results.IRdf= 0;
results.UVdf = results.IRdf*2;
results.IRDet = results.IRdet-results.IRdf;
results.UVDet = results.UVdet-results.UVdf;
% results.IRAbsfreq = 417.20937+results.IRDet*1e-4; % R_1(1) F=2 전이
results.IRAbsfreq = 417.186005+results.IRDet*1e-4; % R_1(0) F=1 전이
results.UVAbsfreq = results.IRAbsfreq*2;

% results.v = -results.Det/(446809870)*299792458/cos(pi/4); % for XY beam
% results.v = -results.Det/(446809870)*299792458/cos(pi/4+atan(1.5/10)); % for XYt beam
% results.v = -results.Det/(446809870)*299792458; % for Z beam

% results.maxv = 0.373./(results.t(results.baselinerange:end)*1e-6); % possible maximum velocity
% results.mint = 0.373./results.v*1e6; % possible minimum arrival time, µs

% NormFlDatas2 = results.FN; % normalized with first fls data
% SumFlDatas2 = zeros(dfsize, iter); % summatino of each iteration of fls time signal
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
% results.FSte2 = std(results.FS2,0,2)/sqrt(iter);

% results.absfit = GF_Li_v2(results, 1); % abs spectrum gaussian fit

if plotmode == true
    plotfigureset_v3(results, results.savename, 'UV')
end


end