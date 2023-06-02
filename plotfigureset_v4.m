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
    % figure('Name',results.savename+" fl spec 2")
    % plotspec_v3(gca, results, 'fl2', band);
    % saveas(gcf, './K_results/'+results.savename+'_fl_spectrum-2.png');
    % saveas(gcf, './K_results/'+results.savename+'_fl_spectrum-2.fig');
    %
    %
    % fig=figure('Name',results.savename+" absfl spec")
    % fig.Position = [120 130 2*560 1*420]
    % tiledlayout(1,2)
    % plotspec_v3(nexttile, results, 'abs', band)
    % % plotspec_v2(nexttile, results, 'fl1')
    % plotspec_v2(nexttile, results, 'fl2')
    % saveas(gcf, './K_results/'+results.savename+'_absfl_spectrum.png');
    % saveas(gcf, './K_results/'+results.savename+'_absfl_spectrum.fig');
    %
end
% 
figure('Name',data.savename+" abs 2D");
plot2Dtimedet_v3(gca, data, 'abs');
saveas(gcf, './K_results/'+data.savename+'_abs_2D.png');
saveas(gcf, './K_results/'+data.savename+'_abs_2D.fig');
% 
%
figure('Name',data.savename+" fl 2D");
plot2Dtimedet_v3(gca, data, 'fl1');
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
% % figure('Name',savename+" fl 2D")
% % plot2Dtimedet_v2(gca, results, 'fl2');
% % saveas(gcf, './K_results/'+savename+'_fl_2D_2.png');
% % saveas(gcf, './K_results/'+savename+'_fl_2D_2.fig');


end