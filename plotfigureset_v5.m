function results = plotfigureset_v5(data,band)
clf;

if data.size.freq ~= 1

    figure('Name',data.savename+" abs spec")
    plotspec_v4(gca, data, 'abs', band);
    SaveFigToFile_v2(gcf,"K_results",data.savename+"_abs_spectrum")
    SaveFigToFile_v2(gcf,"K_results",data.savename+"_abs_spectrum")

    figure('Name',data.savename+" fl spec")
    plotspec_v4(gca, data, 'fl1', band);
    SaveFigToFile_v2(gcf,"K_results",data.savename+"_fl_spectrum");

    if data.beamaxis == 'xy'
        figure('Name',data.savename+" fl spec 2")
        plotspec_v4(gca, data, 'fl2', band);
        SaveFigToFile_v2(gcf,"K_results",data.savename+"_fl_spectrum-2");
    end

    % fig=figure('Name',data.savename+" absfl spec")
    % fig.Position = [120 130 2*560 1*420]
    % tiledlayout(1,2)
    %
    % plotspec_v4(nexttile, data, 'abs', band)
    %
    % plotspec_v4(nexttile, data, 'fl2', band)
    % saveas(gcf, './K_results/'+data.savename+'_absfl_spectrum.png');
    % saveas(gcf, './K_results/'+data.savename+'_absfl_spectrum.fig');

elseif data.size.freq == 1

    figure('Name',data.savename+" abs tt");
    plottt_v1(gca, data,'abs',1);

    xlim([0.08 18]);
    SaveFigToFile_v2(gcf, "K_results",data.savename+"_abs_tt");


    figure('Name',data.savename+" fl tt");
    plottt_v1(gca, data,'fl',1);
    xlim([0.08 18]);
    SaveFigToFile_v2(gcf, "K_results",data.savename+"_fl_tt");

end


%
figure('Name',data.savename+" abs 2D");
plot2Dtimedet_v4(gca, data, 'abs', 0.08, 18);
SaveFigToFile_v2(gcf, "K_results",data.savename+"_abs_2D");

%
%
figure('Name',data.savename+" fl 2D");
plot2Dtimedet_v4(gca, data, 'fl1', 0.08, 18);
SaveFigToFile_v2(gcf, "K_results",data.savename+"_fl_2D");


if data.beamaxis=='xy'
    figure('Name',data.savename+" fl 2D")
    plot2Dtimedet_v3(gca, data, 'fl2');
    saveas(gcf, './K_results/'+data.savename+'_fl_2D_2.png');
    saveas(gcf, './K_results/'+data.savename+'_fl_2D_2.fig');
end

end