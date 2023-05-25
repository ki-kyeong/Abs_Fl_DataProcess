function result = plotspec_v3(ax, data, kind, band)
switch band
    case 'UV'
        % det = data.UVrealdetM;
        det = data.UVfittedDetM;

    case 'IR'
        % det = data.IRrealdetM;
        det = data.UVfittedDetM/2;

    otherwise
        error('try UV or IR');

end


switch kind
    case 'abs'
        result = errorbar(ax, det, data.AM, data.ASte);
        ylabel("Abs signa (a.u.)");
        title("abs spectrum");
        xlabel("det (MHz)");

    case 'fl1'
        result = errorbar(ax, det, data.FM, data.FSte);
        ylabel("Fl signa (a.u.)");
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