function results = plotfigureset_v2(results, savename)
clf;

figure('Name',savename+" abs spec")
plotspec_v2(gca, results, 'abs');
saveas(gcf, './K_results/'+savename+'_abs_spectrum.png');
saveas(gcf, './K_results/'+savename+'_abs_spectrum.fig');


% figure('Name',savename+" fl spec")
% plotspec_v2(gca, results, 'fl1');
% saveas(gcf, './K_results/'+savename+'_fl_spectrum.png');
% saveas(gcf, './K_results/'+savename+'_fl_spectrum.fig');
% 
% 
% figure('Name',savename+" fl spec 2")
% plotspec_v2(gca, results, 'fl2');
% saveas(gcf, './K_results/'+savename+'_fl_spectrum-2.png');
% saveas(gcf, './K_results/'+savename+'_fl_spectrum-2.fig');
% 
% 
% fig=figure('Name',savename+" absfl spec")
% fig.Position = [120 130 2*560 1*420]
% tiledlayout(1,2)
% plotspec_v2(nexttile, results, 'abs')
% % plotspec_v2(nexttile, results, 'fl1')
% plotspec_v2(nexttile, results, 'fl2')
% saveas(gcf, './K_results/'+savename+'_absfl_spectrum.png');
% saveas(gcf, './K_results/'+savename+'_absfl_spectrum.fig');
% 
% 
% figure('Name',savename+" abs 2D")
% plot2Dtimedet_v2(gca, results, 'abs');
% saveas(gcf, './K_results/'+savename+'_abs_2D.png');
% saveas(gcf, './K_results/'+savename+'_abs_2D.fig');
% 
% %
% figure('Name',savename+" fl 2D")
% plot2Dtimedet_v2(gca, results, 'fl1');
% saveas(gcf, './K_results/'+savename+'_fl_2D.png');
% saveas(gcf, './K_results/'+savename+'_fl_2D.fig');
% 
% figure('Name',savename+" fl 2D")
% plot2Dtimedet_v2(gca, results, 'fl2');
% saveas(gcf, './K_results/'+savename+'_fl_2D_2.png');
% saveas(gcf, './K_results/'+savename+'_fl_2D_2.fig');


end