function result = plotspec_v3(ax, data, kind, band)
switch band
    case 'UV'
        det = data.UVDet;

    case 'IR'
        det = data.IRDet;

    otherwise
        error('try UV or IR');

end


switch kind
    case 'abs'
        result = errorbar(ax, det, data.AM, data.ASte);
        title("abs spectrum");
        xlabel("det (MHz)");

    case 'fl1'
        result = errorbar(ax, det, data.FM, data.FSte);
        title("fl spectrum");
        xlabel("det (MHz)");

    case 'fl2'
        %         result = errorbar(ax, data.Det, data.FM2, data.FSte2, 'LineWidth',1);
        result = errorbar(ax, data.v, data.FM2, data.FSte2);

        title("fl spectrum-2");
        % xlabel("velocity (m/s)",'fontsize',16);


    otherwise
        error('try abs or fl or fl2');
end





end