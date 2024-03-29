function result = plot2Dtimedet_v2(ax, Data, kind)

switch kind
    case 'abs'
        data = mean(Data.AN, [2 3]);
    case 'fl1'
        data = mean(Data.FN, [2 3]);
    case 'fl2'
        data = mean(Data.FN2, 2);
    otherwise
        error('try abs or fl')
end
%
%
% result = image(ax, Data.t(Data.baselinerange:end), Data.v(1:89),...
%     reshape(data(Data.baselinerange:end,:,1:89), size(data(Data.baselinerange:end,:,1:89),1), size(data(Data.baselinerange:end,:,1:89),3)).',...
%     'CDataMapping','scaled'); hold on;


% 
% result = image(ax, Data.t(Data.baselinerange:end)*1e-3, Data.v,...
%     reshape(data(Data.baselinerange:end,:,:),...
%     size(data(Data.baselinerange:end,:,:),1), size(data(Data.baselinerange:end,:,:),3)).',...
%     'CDataMapping','scaled'); hold on;

result = image(ax, Data.t(Data.baselinerange:end)*1e-3, Data.UVfittedDetM,...
    reshape(data(Data.baselinerange:end,:,:),...
    size(data(Data.baselinerange:end,:,:),1), size(data(Data.baselinerange:end,:,:),3)).',...
    'CDataMapping','scaled'); hold on;

set(ax, 'Ydir','normal');
% xlim([0 6000]*1e-3);

% ylabel("velocity (m/s)",'fontsize',16);
ylabel("det (MHz)",'fontsize',16);
xlabel("time (ms)", 'FontSize',16);

colorbar

% % 

% result = image(ax, Data.t(Data.baselinerange:end)*1e-3, Data.UVrealdet,...
%     reshape(data(Data.baselinerange:end,:,:),...
%     size(data(Data.baselinerange:end,:,:),1), size(data(Data.baselinerange:end,:,:),3)).',...
%     'CDataMapping','scaled'); hold on;

% switch kind
%     case 'abs'
% 
%     case 'fl1'
%         plot(Data.t(Data.baselinerange:end)*1e-3,Data.maxv, '-.w', 'LineWidth',0.8); hold off;
% 
%     case 'fl2'             
%         plot(Data.t(Data.baselinerange:end)*1e-3,Data.maxv, '-.w', 'LineWidth',0.8); hold off;
% end


end