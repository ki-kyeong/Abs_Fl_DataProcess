function [det, TotalDatas, NormDatas, SumDatas, MeanDatas, SteDatas] = abs2spec(name, freq_i, freq_f, freq_step, trial)
% haha I got new macbook~
det = freq_i : freq_step : freq_f;

TotalDatas = cell(size(det));
NormDatas = cell(size(TotalDatas));
SumDatas = cell(size(det));
MeanDatas = zeros(size(det));
SteDatas = zeros(size(det));

Data = readmatrix(name+"1_"+string(freq_i)+"_MHz.csv");
data_num = size(Data,1);

pedetalrange = 96;
absstartpoint = 100;


for j = 1 : size(det,2)
    f = det(j)
    Datas = zeros(data_num,trial);
    for i = 1:trial
        Data = readmatrix(name+string(i)+"_"+string(f)+"_MHz.csv");
        Datas(:,i) = Data(:,2);
    end
    TotalDatas{j} = Datas;
end

for i = 1: size(det,2)
    d = TotalDatas{i};
    v = [];
    for j = 1: trial
        v = [v d(:,j)/mean(d(1:pedetalrange,j))];
    end
    NormDatas{i} = v;
end

for i = 1 : size(det,2)
    data = NormDatas{i};
    v = [];
    for j = 1 : trial
        v = [v sum(data(absstartpoint:end,j))];
    end
    SumDatas{i} = v;
end

for i = 1: size(det,2)
    MeanDatas(i) = mean(SumDatas{i});
    SteDatas(i) = std(SumDatas{i})/sqrt(trial);
end

end
