function results = plotfigureset_v4(data, band)
clf;

% savename = data.savename;

if data.size.freq ~= 1

    figure('Name',data.savename+" abs spec")
    plotspec_v4(gca, data, 'abs', band);
    saveas(gcf, './K_results/'+data.savename+'_abs_spectrum.png');
    saveas(gcf, './K_results/'+data.savename+'_abs_spectrum.fig');


    figure('Name',data.savename+" fl spec")
    plotspec_v4(gca, data, 'fl1', band);
    saveas(gcf, './K_results/'+data.savename+'_fl_spectrum.png');
    saveas(gcf, './K_results/'+data.savename+'_fl_spectrum.fig');
    %
    %
    % figure('Name',data.savename+" fl spec 2")
    % plotspec_v4(gca, data, 'fl2', band);
    % saveas(gcf, './K_results/'+data.savename+'_fl_spectrum-2.png');
    % saveas(gcf, './K_results/'+data.savename+'_fl_spectrum-2.fig');


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
    plot(data.t*1e-3, data.abs.tt);
    xlim([0.08 10]);
    saveas(gcf, './K_results/'+data.savename+'_abs_tt.png');
    saveas(gcf, './K_results/'+data.savename+'_abs_tt.fig');

    figure('Name',data.savename+" fl tt");
    plot(data.t*1e-3, data.fl.tt);
    xlim([1 10]);
    saveas(gcf, './K_results/'+data.savename+'_fl_tt.png');
    saveas(gcf, './K_results/'+data.savename+'_fl_tt.fig');
end


% 
figure('Name',data.savename+" abs 2D");
plot2Dtimedet_v4(gca, data, 'abs', 0.08, 10);
saveas(gcf, './K_results/'+data.savename+'_abs_2D.png');
saveas(gcf, './K_results/'+data.savename+'_abs_2D.fig');
% 
%
figure('Name',data.savename+" fl 2D");
plot2Dtimedet_v4(gca, data, 'fl1', 0.5, 10);
saveas(gcf, './K_results/'+data.savename+'_fl_2D.png');
saveas(gcf, './K_results/'+data.savename+'_fl_2D.fig');

% fig = figure('Name',results.savename+"abs_fl 2D");
% fig.Position = [476,360,906,420];
% tiledlayout(1,2)
% nexttile
% plot2Dtimedet_v2(gca, results, 'abs');
% 
% nexttile
% plot2Dtimedet_v2(gca, results, 'fl1');
% 
% saveas(gcf, './K_results/'+savename+'_absfl_2D.png');
% saveas(gcf, './K_results/'+savename+'_absfl_2D.fig');
% 
% 
% figure('Name',data.savename+" fl 2D")
% plot2Dtimedet_v3(gca, data, 'fl2');
% saveas(gcf, './K_results/'+data.savename+'_fl_2D_2.png');
% saveas(gcf, './K_results/'+data.savename+'_fl_2D_2.fig');


end