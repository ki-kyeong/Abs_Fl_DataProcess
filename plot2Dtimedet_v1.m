function result = plot2Dtimedet_v1(ax, Data, kind)

switch kind
    case 'abs'
        data = mean(Data.AN,2);
    case 'fl' 
        data = mean(Data.FN, 2);
    otherwise
        error('try abs or fl')
end
% 
% 
% result = image(ax, Data.t(Data.baselinerange:end), Data.v(1:89),...
%     reshape(data(Data.baselinerange:end,:,1:89), size(data(Data.baselinerange:end,:,1:89),1), size(data(Data.baselinerange:end,:,1:89),3)).',...
%     'CDataMapping','scaled'); hold on;



result = image(ax, Data.t(Data.baselinerange:end)*1e-3, Data.v,...
    reshape(data(Data.baselinerange:end,:,:),...
    size(data(Data.baselinerange:end,:,:),1), size(data(Data.baselinerange:end,:,:),3)).',...
    'CDataMapping','scaled'); hold on;

switch kind
    case 'fl'
        
plot(Data.t(Data.baselinerange:end)*1e-3,0.373./(Data.t(Data.baselinerange:end)*1e-6), '-.w', 'LineWidth',0.8); hold off;
end

set(ax, 'Ydir','normal');
xlim([0 6000]*1e-3);

ylabel("velocity (m/s)",'fontsize',16);
xlabel("time (ms)", 'FontSize',16);

colorbar

end