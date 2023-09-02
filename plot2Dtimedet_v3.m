function result = plot2Dtimedet_v3(ax, Data, kind)
set(groot, 'defaultAxesFontSize',18)
set(groot, 'defaultLineLineWidth',2)
if Data.size.freq == 1
    switch kind
        case 'abs'
            result = imagesc(ax, Data.t*1e-3, 1:Data.size.rep, Data.abs.norm.');
        case 'fl1'
            result = imagesc(ax, Data.t*1e-3, 1:Data.size.rep, Data.fl.norm.');
        case 'fl2'
            result = imagesc(ax, Data.t*1e-3, 1:Data.size.rep, Data.fl.norm2.');
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
            result = imagesc(ax, Data.t*1e-3, Data.det.UV.wm.mean, Data.abs.tt.')
            ylabel("detuning (MHz)",'fontsize',16);
        case 'fl1'
            % result = imagesc(ax, Data.t*1e-3, Data.v, Data.fl.tt.'); hold on;
            result = imagesc(ax, Data.t*1e-3, Data.det.UV.wm.mean, Data.fl.tt.'); 
            % hold on;
            % plot(Data.t(Data.baselineidx:end)*1e-3,Data.maxv+2*Data.vshift, '-.w', 'LineWidth',0.8); hold on;
            % plot(Data.t(Data.baselineidx:end)*1e-3,Data.maxv+Data.vshift, '-.w', 'LineWidth',0.8); hold on;
            % plot(Data.t(Data.baselineidx:end)*1e-3,Data.maxv, '-.w', 'LineWidth',0.8); hold on;
            % plot(Data.t(Data.baselineidx:end)*1e-3,Data.maxv-Data.vshift, '-.w', 'LineWidth',0.8); hold on;
            % plot(Data.t(Data.baselineidx:end)*1e-3,Data.maxv-2*Data.vshift, '-.w', 'LineWidth',0.8); hold off;
            % ylabel("velocity (m/s)",'fontsize',16);
            ylabel("detuning from Q_{12}(1) (MHz)",'fontsize',16);
        case 'fl2'
            result = imagesc(ax, Data.t*1e-3, Data.v, Data.fl.tt2.'); hold on;
            % plot(Data.t(Data.baselineidx:end)*1e-3,Data.maxv+2*Data.vshift, '-.w', 'LineWidth',0.8); hold on;
            % plot(Data.t(Data.baselineidx:end)*1e-3,Data.maxv+Data.vshift, '-.w', 'LineWidth',0.8); hold on;
            % plot(Data.t(Data.baselineidx:end)*1e-3,Data.maxv, '-.w', 'LineWidth',0.8); hold on;
            % plot(Data.t(Data.baselineidx:end)*1e-3,Data.maxv-Data.vshift, '-.w', 'LineWidth',0.8); hold on;
            % plot(Data.t(Data.baselineidx:end)*1e-3,Data.maxv-2*Data.vshift, '-.w', 'LineWidth',0.8); hold off;
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

% result = imagesc(ax, Data.t(Data.baselinerange:end)*1e-3, Data.UVrealdet,...
%     reshape(data(Data.baselinerange:end,:,:),...
%     size(data(Data.baselinerange:end,:,:),1), size(data(Data.baselinerange:end,:,:),3)).',...
%     'CDataMapping','scaled'); hold on;

% switch kind
%     case 'abs'
% 
%     case 'fl1'
%         plot(Data.t(Data.baselineidx:end)*1e-3,Data.maxv, '-.w', 'LineWidth',0.8); hold off;
% 
%     case 'fl2'             
%         plot(Data.t(Data.baselineidx:end)*1e-3,Data.maxv, '-.w', 'LineWidth',0.8); hold off;
% end


end