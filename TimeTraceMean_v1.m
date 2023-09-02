function result = TimeTraceMean_v1(data, type)
switch type
    case 'abs'
        result = reshape(mean(data.abs.norm, [2 3]),data.size.time, data.size.freq);
        
    case 'fl'
        result = reshape(mean(data.fl.norm, [2 3]),data.size.time, data.size.freq);

end