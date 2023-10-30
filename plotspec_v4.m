function result = plotspec_v4(ax, data, kind, band)
set(groot, 'defaultLineLineWidth',1.5)
set(groot, 'defaultAxesFontSize',18)
set(groot, 'defaultAxesLineWidth',1)
pbaspect([1.5,1,1])




switch band
    case 'UV'
        % det = data.UVrealdetM;
        det = data.det.UV.fit.mean;

    case 'IR'
        % det = data.IRrealdetM;
        det = data.det.UV.fit.mean/2;

    otherwise
        error('try UV or IR');

end


switch kind
    case 'abs'
        result = plot(ax, det, data.abs.mean);
        % result = errorbar(ax, det, data.abs.mean, data.abs.ste);
        % result = errorbar(ax, det, data.abs.mean, data.abs.ste, ".", "MarkerSize",8);
        ylabel("Abs. signal (a.u.)");
        % title("abs spectrum");
        xlabel("detuning (MHz)");

    case 'fl1'
        result = errorbar(ax, det, data.fl.mean, data.fl.ste);
        % result = errorbar(ax, det, data.fl.mean, data.fl.ste, ".", "MarkerSize",8);
        ylabel("Fl. signal (a.u.)");
        % title("fl spectrum");
        xlabel("detuning (MHz)");

    case 'fl2'
        %         result = errorbar(ax, data.Det, data.FM2, data.FSte2, 'LineWidth',1);
        result = errorbar(ax, data.v, data.fl.mean2, data.fl.ste2);
        
        % title("fl spectrum-2");
        xlabel("velocity (m/s)",'fontsize',16);


    otherwise
        error('try abs or fl1 or fl2');
end


end