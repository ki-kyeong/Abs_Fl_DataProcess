function result = AverageRuns_v2(data, runs)

for i = 1: length(runs)
    d(i) = load('./K_results/'+data+string(runs(i))+'.mat');
end


result.savename = d(1).savename+"_merge";
result.baselineidx = d(1).baselineidx;
result.dt = d(1).dt;
result.freq = d(1).freq;
result.det = d(1).det;
result.t = d(1).t;
result.beamaxis = d(1).beamaxis;
result.abs.sumendtime = d(1).abs.sumendtime;
result.fl.sumendtime = d(1).fl.sumendtime;


% abs data processing
% results.abs.norm = cell2mat(arrayfun( @(x) x.abs.norm, data(runs), UniformOutput=false));

result.abs.norm = [];
for i = 1: length(runs)
    if i == 1
        result.abs.norm(:,:,(1:d(i).size.iter),:) = d(i).abs.norm(:,:,:,:);
    else
        result.abs.norm(:,:,end+(1:d(i).size.iter),:) = d(i).abs.norm(:,:,:,:);
    end
end


[result.size.time result.size.rep result.size.iter result.size.freq]=size(result.abs.norm);

result.abs.tt = TimeTraceMean_v1(result,'abs');
result.abs.sum = TimeTraceSum_v1(result, 'abs', result.abs.sumendtime);

result.abs.mean = squeeze(mean(result.abs.sum,[2 3])); % results.repititionPerStep, results.iteration 전부 평균
result.abs.ste = std(result.abs.sum,0,[2 3])/sqrt(result.size.iter*result.size.rep);

% fl data processing
% results.fl.norm = cell2mat(arrayfun( @(x) x.fl.norm, data(runs), UniformOutput=false));

result.fl.norm = [];
for i = 1:length(runs)
    if i == 1
        result.fl.norm(:,:,(1:d(i).size.iter),:) = d(i).fl.norm(:,:,:,:);
    else
        result.fl.norm(:,:,end+(1:d(i).size.iter),:) = d(i).fl.norm(:,:,:,:);
    end
end


result.fl.tt = TimeTraceMean_v1(result,'fl');
result.fl.sum = TimeTraceSum_v1(result, 'fl', result.fl.sumendtime);

result.fl.mean = squeeze(mean(result.fl.sum,[2 3])); % results.repititionPerStep, results.iteration 전부 평균
result.fl.ste = std(result.fl.sum,0,[2 3])/sqrt(result.size.iter*result.size.rep);

%% Frequency compasation

% abs spectrum fit으로 Det를 정할때

% results.abs.fit = GF_MgF_v2(results); % abs spectrum gaussian fit
% results.det.UV.fit.raw = results.det.UV.wm.raw -results.abs.fit.b3; % b3 = F=0 detuning

result.det.UV.fit.raw.P = result.det.UV.wm.raw.P; % 그냥 적어줄 때
result.det.UV.fit.raw.Q = result.det.UV.wm.raw.Q; % 그냥 적어줄 때

result.det.UV.fit.mean.P = mean(result.det.UV.fit.raw.P, 2);
result.det.UV.fit.ste.P = std(result.det.UV.fit.raw.P,0,2)/sqrt(result.size.iter);

result.det.UV.fit.mean.Q = mean(result.det.UV.fit.raw.Q, 2);
result.det.UV.fit.ste.Q = std(result.det.UV.fit.raw.Q,0,2)/sqrt(result.size.iter);

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
    result.v = -result.det.UV.fit.mean.P/(2*417147178)*299792458/cos(pi/4); % for XY beam
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

    plotfigureset_v5(result,'UV')
end