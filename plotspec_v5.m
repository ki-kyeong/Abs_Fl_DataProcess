function result = plotspec_v5(ax, data, kind)

% UV wavelength만 읽기 시기시작 24/04/30

set(groot, 'defaultAxesFontSize',18)
set(groot, 'defaultLineLineWidth',1)
set(groot, 'defaultAxesFontSize',18)
set(groot, 'defaultAxesLineWidth',1.5)
set(groot, 'defaultAxesBox','on')
pbaspect([1.5,1,1])

det = data.det.fit.mean;

switch kind
    case 'abs'
        % result = plot(ax, det, data.abs.mean);
        result = errorbar(ax, det.Q, data.abs.mean, data.abs.ste);
        % result = errorbar(ax, det, data.abs.mean, data.abs.ste, ".", "MarkerSize",8);
        ylabel("Abs. signal (a.u.)");
        % title("abs spectrum");
        xlabel("detuning from Q_{12}(1) (MHz)");
        

    case 'fl1'
        if data.beamaxis == 'z'
            result = errorbar(ax, det.Q, data.fl.mean, data.fl.ste);
            % result = errorbar(ax, det, data.fl.mean, data.fl.ste, ".", "MarkerSize",8);
            ylabel("Fl. signal (a.u.)");
            % title("fl spectrum");
            xlabel("detuning from Q_{12}(1) (MHz)");
                    


        elseif data.beamaxis == 'xy'
            result = errorbar(ax, det.P, data.fl.mean, data.fl.ste);
            % result = errorbar(ax, det, data.fl.mean, data.fl.ste, ".", "MarkerSize",8);
            ylabel("Fl. signal (a.u.)");
            % title("fl spectrum");
            xlabel("detuning from P_{1}(1) (MHz)");
        end

    case 'fl2'
        %         result = errorbar(ax, data.Det, data.FM2, data.FSte2, 'LineWidth',1);
        result = errorbar(ax, data.v, data.fl.mean2, data.fl.ste2);

        % title("fl spectrum-2");
        xlabel("velocity(P_1(1) ref) (m/s)",'fontsize',16);


    otherwise
        error('try abs or fl1 or fl2');
end


end