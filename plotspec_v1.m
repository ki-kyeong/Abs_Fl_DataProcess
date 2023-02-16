function result = plotspec_v1(ax, data, kind)

switch kind
    case 'abs'
        result = errorbar(ax, data.Det, data.AM, data.ASte, 'LineWidth',1);
        title("abs spectrum",'fontsize',16);
        
    case 'fl'
        result = errorbar(ax, data.Det, data.FM, data.FSte, 'LineWidth',1);
        title("fl spectrum",'fontsize',16);
        
    otherwise
        error('try abs or fl');
end

xlabel("det from F=2 (MHz)",'fontsize',16);


end