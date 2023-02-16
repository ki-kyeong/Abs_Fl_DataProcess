function result = plotspec_v2(ax, data, kind)

switch kind
    case 'abs'
        result = errorbar(ax, data.Det, data.AM, data.ASte);
        title("abs spectrum");
        
    case 'fl1'
        result = errorbar(ax, data.Det, data.FM, data.FSte);
        title("fl spectrum");
        
            case 'fl2'
%         result = errorbar(ax, data.Det, data.FM2, data.FSte2, 'LineWidth',1);
                result = errorbar(ax, data.v, data.FM2, data.FSte2);

        title("fl spectrum-2");
        
    otherwise
        error('try abs or fl');
end

xlabel("det (MHz)");
% xlabel("velocity (m/s)",'fontsize',16);



end