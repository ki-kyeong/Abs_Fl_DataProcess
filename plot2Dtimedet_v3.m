function result = plot2Dtimedet_v3(ax, Data, kind)

if Data.size.freq == 1
    switch kind
        case 'abs'
            result = image(ax, Data.t*1e-3, 1:Data.size.rep, Data.abs.norm.', 'CDataMapping', 'Scaled');
        case 'fl1'
            result = image(ax, Data.t*1e-3, 1:Data.size.rep, Data.fl.norm.', 'CDataMapping', 'Scaled');
        otherwise
            error('try abs or fl')
    end

    set(ax, 'Ydir','normal');
    % xlim([0 6000]*1e-3);
    ylabel("rep #",'fontsize',16);
    xlabel("time (ms)", 'FontSize',16);

colorbar

else
    
    switch kind
        case 'abs'
            result = image(ax, Data.t*1e-3, Data.det.UV.wm.mean, Data.abs.tt.', 'CDataMapping', 'Scaled')
            ylabel("detuning (MHz)",'fontsize',16);
        case 'fl1'
            result = image(ax, Data.t*1e-3, Data.v, Data.fl.tt.', 'CDataMapping', 'Scaled')
            ylabel("velocity (m/s)",'fontsize',16);
        otherwise
            error('try abs or fl')
    end
    set(ax, 'Ydir','normal');
    % xlim([0 6000]*1e-3);
    xlabel("time (ms)", 'FontSize',16);
end



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