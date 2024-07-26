function result = readMgFSpecData_v10(name, plotmode)

% rep의 첫번째 데이터는 버리도록 수정 24/03/18
% UV wavelength만 읽기 시기시작 24/04/30

clf;
close all;

%===========================
result.normtime.start = 10;

%===========================
% results.sumtime.end = results.t(end)*1e-3;
result.abs.sumendtime = 18; % ms
result.fl.sumendtime = 8; % ms

%===========================
beamaxisopt.Default = 'z';
beamaxisopt.Interpreter = 'none';
result.beamaxis = questdlg('choose beamaxis', 'beam axis',...
    'z', 'xy',beamaxisopt);

% result.beamaxis = 'z';

%===========================
result.name = name;
savename = split(name,'_');
result.savename = savename(end-1);

%===========================
result.readme=readlines(result.name+"readme.txt", Encoding="UTF-8");

%======================================
% read parameter data

opt = detectImportOptions(name+"parameter.csv");
expPara = readmatrix(name+"parameter.csv",opt);

for i = 1 : size(opt.VariableNames,2)
    result.(opt.VariableNames{i}) = expPara(i);
end

%===========================
% read run informations

opt2 = detectImportOptions(name+"info.csv");
expPara2 = readmatrix(name+"info.csv",opt2);

for i = 1 : size(opt2.VariableNames,2)
    result.(opt2.VariableNames{i}) = expPara2(i);
end

%===========================
result.size.iter = result.iteration;
result.size.rep = result.repititionPerStep-1;

%======================================
% read each iteration wavemeter data
% first column = cavity freq
% second column = wavemeter freq

for i = 1 : result.size.iter
    result.wavemeter_data{i} = readmatrix(name+string(i-1)+"_wavemeter_data.csv");
end

%===========================
result.freq.cavity = result.wavemeter_data{:,1}(:,1);

for i = 1: result.size.iter
    result.freq.wm.raw(:,:,i) = result.wavemeter_data{:,i}(:,3:end);
end

% freq size
result.size.freq = size(result.freq.cavity,1);

% random array sorting
[result.freq.cavity, result.randidx]=sort(result.freq.cavity);


% wavemeter freq sorting and analysis
result.freq.wm.raw = result.freq.wm.raw(result.randidx,3:end,:);
result.freq.wm.mean = mean(result.freq.wm.raw,[2 3]);
result.freq.wm.ste = std(result.freq.wm.raw,0,[2 3])/sqrt(result.size.iter*result.size.rep);

% detuning cal
result = detcal_v2(result);


%======================================
% read data for initializing

DataParams = detectImportOptions(name+'0_'+num2str(result.freq.cavity(1),'%.6f')+".csv");
Data = readmatrix(name+'0_'+num2str(result.freq.cavity(1),'%.6f')+".csv",DataParams);

result.dt = (Data(2,1)-Data(1,1))*1e6; % µs unit

result.size.time = size(Data,1);

TotalAbsDatas = zeros(result.size.time, result.size.rep, result.size.iter, result.size.freq); % abs time signal row data - bacground
TotalAbsPFMDatas = zeros(result.size.time, result.size.rep, result.size.iter, result.size.freq); % abs time signal row data - bacground
TotalFlDatas = zeros(result.size.time,result.size.rep, result.size.iter, result.size.freq); % fls time signal row data


result.baselineidx = round(40/result.dt); % 40 us에서 ablation peak끝
result.t = (Data(:,1)-Data(result.baselineidx,1))*1e6; % µs, 40 µs을 0초로 설정

wb = waitbar(0, ' Getting started');

for j = 1 : result.size.freq
    f = result.freq.cavity(j);

    waitbar(j/result.size.freq, wb, result.savename+newline+...
        num2str(f, '%.6f') + " THz"+newline+...
        num2str(j/result.size.freq*100, '%.1f')+" % done",...
        'WindowStyle','modal');

    for i = 1:result.size.iter
        Data = readmatrix(name+string(i-1)+"_"+num2str(f,'%.6f')+".csv",DataParams);
        TotalAbsDatas(:,1:result.size.rep,i,j) = Data(:,4:2:2*(result.size.rep+1))-expPara(7)*1e-3;
        TotalAbsPFMDatas(:,1:result.size.rep,i,j) = Data(:,2*(2*(result.size.rep+1)+1)+2:2:2*(3*(result.size.rep+1)))-expPara(8)*1e-3;
        TotalFlDatas(:,1:result.size.rep,i,j) = Data(:,2*((result.size.rep+1)+1)+2:2:2*(2*(result.size.rep+1)));
        % TotalFlDatas(:,1:result.size.rep,i,j) = lowpass(Data(:,2*((result.size.rep+1)+1)+2:2:2*(2*(result.size.rep+1))),50e-4/(result.dt*1e-6), 1/(result.dt*1e-6)); % LPF 2 kHz
    end
end

close(wb)

result.abs.raw = TotalAbsDatas;
result.abs.pfm = TotalAbsPFMDatas;
result.fl.raw = TotalFlDatas;

%% normalizing time trace signal

[result.abs.rawnorm, result.abs.rawbase] = TimeTraceNormalization_v1(result,'abs',result.normtime.start);
[result.abs.pfmnorm, result.abs.pfmbase] = TimeTraceNormalization_v1(result,'pfm',result.normtime.start);
result.abs.norm = result.abs.rawnorm-result.abs.pfmnorm;

[result.fl.norm, result.fl.base] = TimeTraceNormalization_v1(result,'fl', result.normtime.start); % norm value is calculated after 8 ms

%% time trace 평균

result.abs.tt = TimeTraceMean_v1(result, 'abs');
result.fl.tt = TimeTraceMean_v1(result, 'fl');

%% summing up the time trace data

result.abs.sum = TimeTraceSum_v1(result,'abs', result.abs.sumendtime);
result.fl.sum = TimeTraceSum_v1(result,'fl', result.fl.sumendtime);

%% Making spectrum data

result.abs.mean = squeeze(mean(result.abs.sum,[2 3])); % results.repititionPerStep, results.iteration 전부 평균
result.abs.ste = std(result.abs.sum,0,[2 3])/sqrt(result.size.iter*result.size.rep);
result.fl.mean = squeeze(mean(result.fl.sum,[2 3]));
result.fl.ste = std(result.fl.sum,0,[2 3])/sqrt(result.size.iter*result.size.rep);

result.abs.mean_iter = squeeze(mean(result.abs.sum,2)); % results.repititionPerStep 평균
result.abs.ste_iter = std(result.abs.sum,0,2)/sqrt(result.size.rep);
result.fl.mean_iter = squeeze(mean(result.fl.sum,2));
result.fl.ste_iter = squeeze(std(result.fl.sum,0,2)/sqrt(result.size.rep));

%% Frequency compasation

% abs spectrum fit으로 Det를 정할때

result.abs.fit = GF_MgF_v3(result); % abs spectrum gaussian fit
result.det.fit.raw.Q = result.det.wm.raw.Q -result.abs.fit.b3; % b3 = F=0 detuning
result.det.fit.raw.P = result.det.wm.raw.P -result.abs.fit.b3; % b3 = F=0 detuning

% result.det.fit.raw.P = result.det.wm.raw.P; % 그냥 적어줄 때
% result.det.fit.raw.Q = result.det.wm.raw.Q; % 그냥 적어줄 때

result.det.fit.mean.P = mean(result.det.fit.raw.P, [2 3]);
result.det.fit.ste.P = std(result.det.fit.raw.P,0,[2 3])/sqrt(result.size.iter*result.size.rep);
result.det.fit.mean.Q = mean(result.det.fit.raw.Q, [2 3]);
result.det.fit.ste.Q = std(result.det.fit.raw.Q,0,[2 3])/sqrt(result.size.iter*result.size.rep);


% z fl spectrum 으로 Det를 정할때
% results.df = results.det(results.FM == max(results.FM));
% results.df= results.det(end);
% results.Det = results.det-results.df;

result.v  = zeros(size(result.freq));
result.maxv = zeros(size(result.t(result.baselineidx:end))); % possible maximum velocity
result.mint = zeros(size(result.v)); % possible minimum arrival time, µs
% results.mint2 = zeros(size(results.v)); % possible minimum arrival time, µs
result.fl.norm2 = result.fl.norm; % normalized with first fls data
result.fl.sum2 = zeros(size(result.fl.sum)); % summatino of each results.iteration of fls time signal
result.fl.tt2 = zeros(result.size.time, result.size.freq);
result.fl.mean2 = zeros(result.size.freq,1);
result.fl.ste2 = zeros(result.size.freq,1);

if result.beamaxis == 'xy'
    result.v = -result.det.fit.mean.P/(834294356)*299792458/cos(pi/4); % for XY beam
    result.maxv = 0.373./(result.t(result.baselineidx:end)*1e-6); % possible maximum velocity
    result.mint = 0.373./result.v*1e6; % possible minimum arrival time, µs
    % results.mint2 = 0.373./(results.v+results.vshift)*1e6; % extract eom freq signal
    for j = 1 : result.size.freq
        if result.mint(j) <=0
            result.fl.norm2(:,:,:,j) = 0;
        else
            result.fl.norm2(result.t<= result.mint(j),:,:,j) = 0;
            % if results.mint2(j)>=0
            %     results.fl.norm2(results.t>= results.mint2(j),:,:,j) = 0;
            % end
        end
    end

    result.fl.tt2 = reshape(mean(result.fl.norm2, [2 3]), result.size.time, result.size.freq);

    for j = 1 : result.size.freq
        result.fl.sum2(j,:,:) = sum(result.fl.norm2(result.baselineidx+round(80/result.dt):result.baselineidx+6000/result.dt,:,:,j));
    end

    result.fl.mean2 = mean(result.fl.sum2,[2 3]);
    result.fl.ste2 = std(result.fl.sum2,0,[2 3])/sqrt(result.size.iter*result.size.rep);

end

% % results.vshift = -120/(2*417147242.5)*299792458/cos(pi/4); % for XY beam

if plotmode == true
    plotfigureset_v6(result)
end

% load Handel
% sound(y,1.3*Fs)

% save('DataFile'+result.name+'.m',"-v7.3");

end
