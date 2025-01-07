function result = SaveFigToFile_v2(fig,path, name)

saveas(fig, path+'/'+name+'.fig');
saveas(fig, path+'/'+name+'.png');


end