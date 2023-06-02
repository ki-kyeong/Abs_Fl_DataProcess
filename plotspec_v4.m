function result = plotspec_v4(ax, data, kind, band)
switch band
    case 'UV'
        % det = data.UVrealdetM;
        det = data.det.UV.fit.mean;

    case 'IR'
        % det = data.IRrealdetM;
        det = data.det.IR.fit.mean;

    otherwise
        error('try UV or IR');

end


switch kind
    case 'abs'
        result = errorbar(ax, det, data.abs.mean, data.abs.ste);
        ylabel("Abs signal (a.u.)");
        title("abs spectrum");
        xlabel("det (MHz)");

    case 'fl1'
        result = errorbar(ax, det, data.fl.mean, data.fl.ste);
        ylabel("Fl signal (a.u.)");
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