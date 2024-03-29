function results = readMgFSpecData_v6(name, normmode, plotmode)

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
SumAbsDatas = zeros(results.size.freq,results.size.rep, results.size.iter);% summation of each results.iteration of normalized abs time signal
TotalFlDatas = zeros(results.size.time,results.size.rep, results.size.iter, results.size.freq); % fls time signal row data
% NormFlDatas = zeros(size(TotalFlDatas)); % normalized with first fls data
SumFlDatas = zeros(results.size.freq,results.size.rep, results.size.iter); % summatino of each results.iteration of fls time signal



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
        TotalAbsDatas(:,1:results.size.rep,i,j) = Data(:,2:2:2*results.size.rep)-expPara(7)*1e-3;
        TotalAbsPFMDatas(:,1:results.size.rep,i,j) = Data(:,2*(2*results.size.rep+1):2:2*(3*results.size.rep))-expPara(7)*1e-3;
        TotalFlDatas(:,1:results.size.rep,i,j) = Data(:,2*(results.size.rep+1):2:2*(2*results.size.rep));
    end
end

close(wb)

results.abs.raw = TotalAbsDatas;
results.abs.pfm = TotalAbsPFMDatas;
results.fl.raw = TotalFlDatas;

switch normmode
    case 'basic'
        results.abs.norm = 1-(results.abs.raw/mean(results.abs.raw(results.baselinerange,1,1)));
        results.fl.norm = results.fl.raw/mean(results.fl.raw(7000/results.dt:end,1,1))-1;

    case 'powernorm'
        for j = 1 : results.size.freq
            for i = 1 : results.size.iter
                for k = 1 : results.size.rep
                nn(:,k,i,j) = 1-(results.abs.raw(:,k,i,j)/mean(results.abs.raw(results.baselinerange,k,i,j)));
                nnsp(:,k,i,j) = 1-(results.abs.pfm(:,k,i,j)/mean(results.abs.pfm(results.baselinerange,k,i,j)));
                % nnn(:,k, i,j) = (results.fl.raw(:,k, i,j)/mean(results.fl.raw(results.baselineidx+11000/results.dt:end,k,i,j)))-1; % time trace 2D image를 그려보고 8 ms으로 정햇음...
                % nnn(:,k, i,j) = (results.fl.raw(:,k, i,j)/mean(results.fl.raw(results.baselinerange,k,i,j)))-1;
                % nnn(:,k, i,j) = results.fl.raw(:,k, i, j)-mean(results.fl.raw(results.baselinerange,k,i,j));
                nnn(:,k, i,j) = (results.fl.raw(:,k, i,j)-mean(results.fl.raw(results.baselineidx+8000/results.dt:end,k,i,j))); % time trace 2D image를 그려보고 8 ms으로 정햇음...
                % nnnsp(:, k,i,j) = results.abs.pfm(:,k,i,j)/mean(results.abs.pfm(results.baselinerange,k,i,j));
                end
            end

        results.abs.rawnorm = nn;
        results.abs.pfmnorm = nnsp;
        results.abs.norm = nn-nnsp;

        results.fl.rawnorm = nnn;
        % results.fl.pfmnorm = nnnsp;
        % results.fl.norm = nnn+nnnsp-1;
        results.fl.norm = nnn;

        end
    otherwise
        error('try basic of powernorm')
end


%% rep들을 미리 평균낸것도 계산해주자.

    results.abs.tt = reshape(mean(results.abs.norm, [2 3]),results.size.time, results.size.freq);
    results.fl.tt = reshape(mean(results.fl.norm, [2 3]), results.size.time, results.size.freq);



for j = 1 : results.size.freq
    SumAbsDatas(j,:,:) = sum(results.abs.norm(results.baselineidx+round(80/results.dt):end,:,:,j)); % ablation 이후 80 µs 부터
    % SumAbsDatas(j,:,:) = sum(results.AN(results.baselinerange+round(80/results.dt):results.baselinerange+10000/results.dt,:,:,j)); % ablation 이후 5 ms 까지
    % SumFlDatas(j,:,:) = sum(results.FN(results.baselinerange+4:end,:,:,j));
    % SumFlDatas(j,:,:) = sum(results.FN(results.baselinerange+80/results.dt:end,:,:,j)); % ablation이 한 40 µs뒤에 끝남
    SumFlDatas(j,:,:) = sum(results.fl.norm(results.baselineidx+1000/results.dt:results.baselineidx+8000/results.dt,:,:,j)); % 한 10 ms까지만 더해보자
    % SumFlDatas(j,:,:) = sum(results.fl.norm(results.baselineidx+round(2000/results.dt):results.baselineidx+6000/results.dt,:,:,j)); % ablation 이후 5 ms 까지
end

results.abs.sum = SumAbsDatas;
results.fl.sum = SumFlDatas;

results.abs.mean = mean(results.abs.sum,[2 3]); % results.repititionPerStep, results.iteration 전부 평균
results.abs.ste = std(results.abs.sum,0,[2 3])/sqrt(results.size.iter*results.size.rep);
results.fl.mean = mean(results.fl.sum,[2 3]);
results.fl.ste = std(results.fl.sum,0,[2 3])/sqrt(results.size.iter*results.size.rep);

% abs spectrum fit으로 Det를 정할때

% results.abs.fit = GF_MgF_v2(results); % abs spectrum gaussian fit
% results.det.UV.fit.raw = results.det.UV.wm.raw -results.abs.fit.b3; % b3 = F=0 detuning

results.det.UV.fit.raw = results.det.UV.wm.raw; % 그냥 적어줄 때

results.det.UV.fit.mean = mean(results.det.UV.fit.raw, 2);
results.det.UV.fit.ste = std(results.det.UV.fit.raw,0,2)/sqrt(results.size.iter);

% z fl spectrum 으로 Det를 정할때
% results.df = results.det(results.FM == max(results.FM));
% results.df= results.det(end);
% results.Det = results.det-results.df;

% 
% results.v = -results.det.UV.fit.mean/(2*417147242.5)*299792458/cos(pi/4); % for XY beam
% % results.vshift = -120/(2*417147242.5)*299792458/cos(pi/4); % for XY beam
% % results.v = -results.det.UV.fit.mean/(2*417147242.5)*299792458/cos(pi/4+atan(1.5/10)); % for XYt beam
% results.v = -results.det.UV.fit.mean/(2*417147242.5)*299792458; % for Z beam
% 
% results.maxv = 0.373./(results.t(results.baselineidx:end)*1e-6); % possible maximum velocity
% results.mint = 0.373./results.v*1e6; % possible minimum arrival time, µs
% % results.mint2 = 0.373./(results.v+results.vshift)*1e6; % extract eom freq signal
% % 
% results.fl.norm2 = results.fl.norm; % normalized with first fls data
% results.fl.sum2 = zeros(size(results.fl.sum)); % summatino of each results.iteration of fls time signal
% 
% for j = 1 : results.size.freq
%    if results.mint(j) <=0
%        results.fl.norm2(:,:,:,j) = 0;
%    else
%        results.fl.norm2(results.t<= results.mint(j),:,:,j) = 0;
%        % if results.mint2(j)>=0
%        %     results.fl.norm2(results.t>= results.mint2(j),:,:,j) = 0;
%        % end
%    end
% 
% end
% 
% results.fl.tt2 = reshape(mean(results.fl.norm2, [2 3]), results.size.time, results.size.freq);
% 
% for j = 1 : results.size.freq
%     results.fl.sum2(j,:,:) = sum(results.fl.norm2(results.baselineidx+round(80/results.dt):results.baselineidx+6000/results.dt,:,:,j));
% end
% 
% 
% 
% results.fl.mean2 = mean(results.fl.sum2,[2 3]);
% results.fl.ste2 = std(results.fl.sum2,0,[2 3])/sqrt(results.size.iter*results.size.rep);


if plotmode == true
    plotfigureset_v4(results,'UV')
end


end
