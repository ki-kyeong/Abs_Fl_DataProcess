function result = TimeTraceSum_v1(data, type, SumEndTime)

switch type
    case 'abs'
        for j = 1 : data.size.freq
            SumData(j,:,:) = sum(data.abs.norm(data.baselineidx+1:data.baselineidx+round(SumEndTime*1e3/data.dt),:,:,j));
        end
    case 'fl'
        for j = 1 : data.size.freq
            SumData(j,:,:) = sum(data.fl.norm(data.baselineidx+1:data.baselineidx+round(SumEndTime*1e3/data.dt),:,:,j));
        end
end

result = SumData;

end