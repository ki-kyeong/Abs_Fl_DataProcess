function result = plot2Dtimedet_v5(ax, Data, kind, StartTime, EndTime)

% UV wavelength만 읽기 시기시작 24/04/30

set(groot, 'defaultAxesFontSize',18)
set(groot, 'defaultLineLineWidth',1)
set(groot, 'defaultAxesFontSize',18)
set(groot, 'defaultAxesLineWidth',1.5)
set(groot, 'defaultAxesBox','on')
pbaspect([1.5,1,1])

StartTimeIdx = Data.baselineidx+StartTime*1e3/Data.dt;
EndTimeIdx = Data.baselineidx+EndTime*1e3/Data.dt;

if Data.size.freq == 1
    switch kind
        case 'abs' 
            dd = reshape(Data.abs.norm,Data.size.time,Data.size.rep*Data.size.iter);
            result = imagesc(ax, Data.t(StartTimeIdx:EndTimeIdx)*1e-3, 1:Data.size.rep*Data.size.iter, dd(StartTimeIdx:EndTimeIdx,:).');
        case 'fl1'
            dd = reshape(Data.fl.norm,Data.size.time,Data.size.rep*Data.size.iter);
            result = imagesc(ax, Data.t(StartTimeIdx:EndTimeIdx)*1e-3, 1:Data.size.rep*Data.size.iter, dd(StartTimeIdx:EndTimeIdx,:).');
        case 'fl2'
            dd = reshape(Data.fl.norm2,Data.size.time,Data.size.rep*Data.size.iter);
            result = imagesc(ax, Data.t(StartTimeIdx:EndTimeIdx)*1e-3, 1:Data.size.rep*Data.size.iter, dd(StartTimeIdx:EndTimeIdx,:).');
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
            result = imagesc(ax, Data.t(StartTimeIdx:EndTimeIdx)*1e-3, Data.det.wm.mean.Q, Data.abs.tt(StartTimeIdx:EndTimeIdx,:).');
            ylabel("detuning from Q_{12}(1) (MHz)",'fontsize',16);
        case 'fl1'
            if Data.beamaxis =='xy'
                % result = imagesc(ax, Data.t(StartTimeIdx:EndTimeIdx)*1e-3, Data.v, Data.fl.tt(StartTimeIdx:EndTimeIdx,:).'); hold on;
                result = imagesc(ax, Data.t(StartTimeIdx:EndTimeIdx)*1e-3, Data.v, Data.fl.tt(StartTimeIdx:EndTimeIdx,:).'); hold on;
                % hold on;
                % plot(Data.t(Data.baselineidx:end)*1e-3,Data.maxv+2*Data.vshift, '-.w', 'LineWidth',0.8); hold on;
                % plot(Data.t(Data.baselineidx:end)*1e-3,Data.maxv+Data.vshift, '-.w', 'LineWidth',0.8); hold on;
                plot(Data.t(Data.baselineidx:end)*1e-3,Data.maxv, '-.w', 'LineWidth',0.8); hold off;
                % plot(Data.t(Data.baselineidx:end)*1e-3,Data.maxv-Data.vshift, '-.w', 'LineWidth',0.8); hold off;
                % plot(Data.t(Data.baselineidx:end)*1e-3,Data.maxv-2*Data.vshift, '-.w', 'LineWidth',0.8); hold off;
                % ylabel("velocity (m/s)",'fontsize',16);

                ylabel("velocity (m/s)",'fontsize',16);

            elseif Data.beamaxis == 'z'
                result = imagesc(ax, Data.t(StartTimeIdx:EndTimeIdx)*1e-3, Data.det.wm.mean.Q, Data.fl.tt(StartTimeIdx:EndTimeIdx,:).'); hold on;
                ylabel("detuning from Q_{12}(1) (MHz)",'fontsize',16);
            end
        case 'fl2'
            result = imagesc(ax, Data.t(StartTimeIdx:EndTimeIdx)*1e-3, Data.v, Data.fl.tt2(StartTimeIdx:EndTimeIdx,:).'); hold on;
            % plot(Data.t(Data.baselineidx:end)*1e-3,Data.maxv+2*Data.vshift, '-.w', 'LineWidth',0.8); hold on;
            % plot(Data.t(Data.baselineidx:end)*1e-3,Data.maxv+Data.vshift, '-.w', 'LineWidth',0.8); hold on;
            % plot(Data.t(Data.baselineidx:end)*1e-3,Data.maxv, '-.w', 'LineWidth',0.8); hold on;
            % plot(Data.t(Data.baselineidx:end)*1e-3,Data.maxv-Data.vshift, '-.w', 'LineWidth',0.8); hold on;
            % plot(Data.t(Data.baselineidx:end)*1e-3,Data.maxv-2*Data.vshift, '-.w', 'LineWidth',0.8); hold off;

            % ylabel("detuning from Q_{12}(1) (MHz)",'fontsize',16);
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