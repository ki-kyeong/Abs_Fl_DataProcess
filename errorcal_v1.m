function result = errorcal_v1(data1, data2, mode)

switch mode
    case "sum"
        result = sqrt(data1.ste.^2+data2.ste.^2);
    case "mul"
        result = abs(data1.mean.*data2.mean).*sqrt((data1.ste./data1.mean).^2+(data2.ste./data2.mean).^2);
end

end