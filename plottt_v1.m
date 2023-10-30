function result = plottt_v1(gca, data, mode, freqidx )
set(groot, 'defaultLineLineWidth',1.5)
set(groot, 'defaultAxesFontSize',18)
set(groot, 'defaultAxesLineWidth',1)
pbaspect([1.5,1,1])

switch mode
    case 'abs'
        plot(data.t*1e-3, data.abs.tt(:,freqidx));
        xlim([0.08 data.t(end)*1e-3])
        ylabel("Abs. signal (a.u.)")
        

    case 'fl'
        plot(data.t*1e-3, data.fl.tt(:,freqidx));
        xlim([0.1 data.t(end)*1e-3])
        ylabel("Fl. signal (V)")
end

xlabel("time (ms)");

end