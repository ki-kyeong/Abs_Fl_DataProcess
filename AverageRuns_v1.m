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
for i = runs
        results.abs.norm(:,:,end+(1:data(i).size.iter),:) = data(i).abs.norm(:,:,:,:);
end


[results.size.time results.size.rep results.size.iter results.size.freq]=size(results.abs.norm);

results.abs.tt = TimeTraceMean_v1(results,'abs');
results.abs.sum = TimeTraceSum_v1(results, 'abs', 0.08, results.t(end)*1e-3);

results.abs.mean = squeeze(mean(results.abs.sum,[2 3])); % results.repititionPerStep, results.iteration 전부 평균
results.abs.ste = std(results.abs.sum,0,[2 3])/sqrt(results.size.iter*results.size.rep);

% fl data processing
% results.fl.norm = cell2mat(arrayfun( @(x) x.fl.norm, data(runs), UniformOutput=false));

results.fl.norm = [];
for i = runs
        results.fl.norm(:,:,end+(1:data(i).size.iter),:) = data(i).fl.norm(:,:,:,:);
end


results.fl.tt = TimeTraceMean_v1(results,'fl');
results.fl.sum = TimeTraceSum_v1(results, 'fl', 0.08, results.t(end)*1e-3);

results.fl.mean = squeeze(mean(results.fl.sum,[2 3])); % results.repititionPerStep, results.iteration 전부 평균
results.fl.ste = std(results.fl.sum,0,[2 3])/sqrt(results.size.iter*results.size.rep);


end