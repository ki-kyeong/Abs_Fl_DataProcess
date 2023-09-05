function result = SaveFigToFile_v1(fig, name)

saveas(fig, './K_results/'+name+'.fig');
saveas(fig, './K_results/'+name+'.png');

end