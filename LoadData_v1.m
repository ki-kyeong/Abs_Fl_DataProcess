function [result, logline] = LoadData_v1(date, runs, figuremode)

for i = runs
    name = "./data/2023-"+date+"_run"+string(i)+"_";
    result(i)= readMgFSpecData_v7(name,figuremode);

    if size(result(i).readme,1) ==0
        logline(i,1) = "empty";
    else
        logline(i,1) = result(i).readme;
    end

end

logline = arrayfun(@(x) x.readme, result, UniformOutput=false)';

% if savemode == 1
%     save('DataFile.mat',"-v7.3");
% end

end
