function results = AverageRuns_v1(data, runs)

results.savename = data.savename;
results.baselineidx = data(runs(1)).baselineidx;
results.dt = data(runs(1)).dt;
results.freq = data(runs(1)).freq;
results.det = data(runs(1)).det;
results.t = data(runs(1)).t;


% abs data processing
% results.abs.norm = cell2mat(arrayfun( @(x) x.abs.norm, data(runs), UniformOutput=false));

results.abs.norm = [];
for i = 1: length(runs)
    if i == 1
        results.abs.norm(:,:,(1:data(runs(i)).size.iter),:) = data(runs(i)).abs.norm(:,:,:,:);
    else
        results.abs.norm(:,:,end+(1:data(runs(i)).size.iter),:) = data(runs(i)).abs.norm(:,:,:,:);
    end
end


[results.size.time results.size.rep results.size.iter results.size.freq]=size(results.abs.norm);

results.abs.tt = TimeTraceMean_v1(results,'abs');
results.abs.sum = TimeTraceSum_v1(results, 'abs', 0.08, results.t(end)*1e-3);

results.abs.mean = squeeze(mean(results.abs.sum,[2 3])); % results.repititionPerStep, results.iteration 전부 평균
results.abs.ste = std(results.abs.sum,0,[2 3])/sqrt(results.size.iter*results.size.rep);

% fl data processing
% results.fl.norm = cell2mat(arrayfun( @(x) x.fl.norm, data(runs), UniformOutput=false));

results.fl.norm = [];
for i = 1:length(runs)
    if i == 1
        results.fl.norm(:,:,(1:data(runs(i)).size.iter),:) = data(runs(i)).fl.norm(:,:,:,:);
    else
        results.fl.norm(:,:,end+(1:data(runs(i)).size.iter),:) = data(runs(i)).fl.norm(:,:,:,:);
    end
end


results.fl.tt = TimeTraceMean_v1(results,'fl');
results.fl.sum = TimeTraceSum_v1(results, 'fl', 0.08, results.t(end)*1e-3);

results.fl.mean = squeeze(mean(results.fl.sum,[2 3])); % results.repititionPerStep, results.iteration 전부 평균
results.fl.ste = std(results.fl.sum,0,[2 3])/sqrt(results.size.iter*results.size.rep);

%% Frequency compasation

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
results.v = -results.det.UV.fit.mean/(2*417147178)*299792458/cos(pi/4); % for XY beam
% % results.vshift = -120/(2*417147242.5)*299792458/cos(pi/4); % for XY beam
% % results.v = -results.det.UV.fit.mean/(2*417147242.5)*299792458/cos(pi/4+atan(1.5/10)); % for XYt beam
% results.v = -results.det.UV.fit.mean/(2*417147242.5)*299792458; % for Z beam
% results.v = -results.det.UV.fit.mean/(2*417147178)*299792458; % for Z beam
% 
results.maxv = 0.373./(results.t(results.baselineidx:end)*1e-6); % possible maximum velocity
results.mint = 0.373./results.v*1e6; % possible minimum arrival time, µs
% % results.mint2 = 0.373./(results.v+results.vshift)*1e6; % extract eom freq signal
% % 
results.fl.norm2 = results.fl.norm; % normalized with first fls data
results.fl.sum2 = zeros(size(results.fl.sum)); % summatino of each results.iteration of fls time signal
% 
for j = 1 : results.size.freq
   if results.mint(j) <=0
       results.fl.norm2(:,:,:,j) = 0;
   else
       results.fl.norm2(results.t<= results.mint(j),:,:,j) = 0;
       % if results.mint2(j)>=0
       %     results.fl.norm2(results.t>= results.mint2(j),:,:,j) = 0;
       % end
   end

end
% 
results.fl.tt2 = reshape(mean(results.fl.norm2, [2 3]), results.size.time, results.size.freq);
% 
for j = 1 : results.size.freq
    results.fl.sum2(j,:,:) = sum(results.fl.norm2(results.baselineidx+round(80/results.dt):results.baselineidx+6000/results.dt,:,:,j));
end



results.fl.mean2 = mean(results.fl.sum2,[2 3]);
results.fl.ste2 = std(results.fl.sum2,0,[2 3])/sqrt(results.size.iter*results.size.rep);



end