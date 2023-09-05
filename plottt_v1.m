function result = plottt_v1(gca, data, mode, freqidx )

switch mode
    case 'abs'
        plot(data.t*1e-3, data.abs.tt(:,freqidx));
        xlim([0.08 data.t(end)*1e-3])
        ylabel("absorption sigal (a.u.)")
        

    case 'fl'
        plot(data.t*1e-3, data.fl.tt(:,freqidx));
        xlim([0.1 data.t(end)*1e-3])
        ylabel("fluorescence sigal (a.u.)")
end

xlabel("time (ms)");

end