function results = plotfigureset_v1(results, savename)
clf;

figure('Name',savename+" abs spec")
plotspec_v1(gca, results, 'abs');
saveas(gcf, './K_results/'+savename+'_abs_spectrum.png');
saveas(gcf, './K_results/'+savename+'_abs_spectrum.fig');


figure('Name',savename+" fl spec")
plotspec_v1(gca, results, 'fl');
saveas(gcf, './K_results/'+savename+'_fl_spectrum.png');
saveas(gcf, './K_results/'+savename+'_fl_spectrum.fig');
% 
% 
fig=figure('Name',savename+" absfl spec")
fig.Position = [120 130 2*560 1*420]
tiledlayout(1,2)
plotspec_v1(nexttile, results, 'abs')
plotspec_v1(nexttile, results, 'fl')
saveas(gcf, './K_results/'+savename+'_absfl_spectrum.png');
saveas(gcf, './K_results/'+savename+'_absfl_spectrum.fig');
% 
% 
% figure('Name',savename+" abs 2D")
% plot2Dtimedet_v1(gca, results, 'abs');
% saveas(gcf, './K_results/'+savename+'_abs_2D.png');
% saveas(gcf, './K_results/'+savename+'_abs_2D.fig');
% 
% 
figure('Name',savename+" fl 2D")
plot2Dtimedet_v1(gca, results, 'fl');
saveas(gcf, './K_results/'+savename+'_fl_2D.png');
saveas(gcf, './K_results/'+savename+'_fl_2D.fig');


end