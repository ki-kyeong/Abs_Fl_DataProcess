classdef MgFSpecData_vtest
    properties
        name
        runnum
        plotmode
       
    end

    methods
        function self = MgFSpecData_vtest(name, plotmode)
            self.name = name;
            temp = split(name,'_');
            self.runnum = temp(end-1);
            self.plotmode = plotmode;

        end %MgFSpecData

        function result = RaedReadme(self)
            result = readlines(self.name+"readme.txt", Encoding="UTF-8");
        end %readreadme

        function self = ReadExpParams(self)

            opt = detectImportOptions(self.name+"parameter.csv");
            param = readmatrix(self.name+"parameter.csv",opt);
            for i = 1 : size(opt.VariableNames,2)
                self.params.(opt.VariableNames{i})= param(i);
            end

            opt2 = detectImportOptions(self.name+"info.csv");
            param2 = readmatrix(self.name+"info.csv",opt2);
            for i = 1 : size(opt2.VariableNames,2)
                self.params.(opt2.VariableNames{i})= param2(i);
            end

            self.size.iter = self.params.iteration;
            self.size.rep = self.params.repititionPerStep;
        end %readparams

        function self = ReadFreqData(self)

            for i = 1 : self.size.iter
                self.freq.wavemeterdata(:,:,i) = readmatrix(self.name+string(i-1)+"_wavemeter_data.csv");
            end

            self.size.freq = size(self.freq.wavemeterdata,1);

            % random array sorting
            self.freq.IR.set = self.freq.wavemeterdata(:,1,1);
            [self.freq.IR.set, randidx] = sort(self.freq.IR.set);

            % cavit UV freq cal
            self.freq.UV.cavity = 2*self.freq.IR.set;

            % wavemeter freq sorting and analysis
            self.freq.IR.wm.raw = self.freq.wavemeterdata(randidx,2:end,:);
            self.freq.IR.wm.mean = mean(self.freq.IR.wm.raw,[2 3]);
            self.freq.IR.wm.ste = std(self.freq.IR.wm.raw,0,[2 3])/sqrt(self.size.iter*self.size.rep);

            self.freq.UV.wm.raw = 2*self.freq.IR.wm.raw;
            self.freq.UV.wm.mean = 2*self.freq.IR.wm.mean;
            self.freq.UV.wm.ste = std(self.freq.UV.wm.raw,0,[2 3])/sqrt(self.size.iter*self.size.rep);

        end %readfreqdata

        function self = CalDetuning(self)
            %% Q12 detuning
            % freqQ11 = 417.170602; % Q1(1) F=1 → F'=1
            freqQ121 = 417.1472425;
            % IR cavity detuning cal
            self.det.IR.set.Q = (self.freq.IR.set-freqQ121)*1e6; % MHz, from Q12(1)
            % UV cavity detuning cal
            self.det.UV.set.Q = 2*self.det.IR.set.Q;

            %IR wavemeter detuning cal and analysis
            self.det.IR.wm.raw.Q = (self.freq.IR.wm.raw-freqQ121)*1e6; % MHz
            self.det.IR.wm.mean.Q = mean(self.det.IR.wm.raw.Q,[2 3]);
            self.det.IR.wm.ste.Q = std(self.det.IR.wm.raw.Q,0,[2 3])/sqrt(self.size.iter);

            %UV wavemeter detuning cal and analysis
            self.det.UV.wm.raw.Q = 2*self.det.IR.wm.raw.Q;
            self.det.UV.wm.mean.Q = 2*self.det.IR.wm.mean.Q;
            self.det.UV.wm.ste.Q = std(self.det.UV.wm.raw.Q,0,[2 3])/sqrt(self.size.iter);

            %% P1 detuning
            % IR cavity detuning cal
            freqP11 = 417.147178;
            self.det.IR.set.P = (self.freq.IR.set-freqP11)*1e6; % MHz, from P1(1), F=2
            % UV cavity detuning cal
            self.det.UV.set.P = 2*self.det.IR.set.P;

            %IR wavemeter detuning cal and analysis
            self.det.IR.wm.raw.P = (self.freq.IR.wm.raw-freqP11)*1e6; % MHz
            self.det.IR.wm.mean.P = mean(self.det.IR.wm.raw.P,[2 3]);
            self.det.IR.wm.ste.P = std(self.det.IR.wm.raw.P,0,[2 3])/sqrt(self.size.iter);

            %UV wavemeter detuning cal and analysis
            self.det.UV.wm.raw.P = 2*self.det.IR.wm.raw.P;
            self.det.UV.wm.mean.P = 2*self.det.IR.wm.mean.P;
            self.det.UV.wm.ste.P = std(self.det.UV.wm.raw.P,0,[2 3])/sqrt(self.size.iter);

        end %caldet

        function self = InitDataStruct(self)
            self.params.DataParams = detectImportOptions(self.name+'0_'+num2str(self.freq.IR.set(1),'%.6f')+".csv");
            self.params.DataParams.DataLines = [3 Inf];
            self.params.DataParams.VariableUnitsLine = 2;
            Data = readmatrix(self.name+'0_'+num2str(self.freq.IR.set(1),'%.6f')+".csv",self.params.DataParams);

            self.dt = (Data(2,1)-Data(1,1))*1e6; % µs unit

            self.size.time = size(Data,1);

            self.data.abs.raw = zeros(self.size.time, self.size.rep, self.size.iter, self.size.freq); % abs time signal row data - bacground
            self.data.abs.pfm = zeros(self.size.time, self.size.rep, self.size.iter, self.size.freq); % abs time signal row data - bacground
            self.data.fl.raw = zeros(self.size.time,self.size.rep, self.size.iter, self.size.freq); % fls time signal row data


            self.baselineidx = round(40/self.dt); % 40 us에서 ablation peak끝
            self.t = (Data(:,1)-Data(self.baselineidx,1))*1e6; % µs, 40 µs을 0초로 설정
        end %initdatastruct

        function self = ReadSpecData(self)
            wb = waitbar(0, ' Getting started');

            for j = 1 : self.size.freq
                f = self.freq.IR.set(j);

                waitbar(j/self.size.freq, wb, self.runnum+newline+...
                    num2str(f, '%.6f') + " THz"+newline+...
                    num2str(j/self.size.freq*100, '%.1f')+" % done",...
                    'WindowStyle','modal');

                for i = 1:self.size.iter
                    Data = readmatrix(self.name+string(i-1)+"_"+num2str(f,'%.6f')+".csv",self.params.DataParams);
                    self.data.abs.raw(:,1:self.size.rep,i,j) = Data(:,2:2+self.size.rep-1)-self.params.AbsBgVoltage_mV_*1e-3;
                    self.data.fl.raw(:,1:self.size.rep,i,j) = Data(:,2+self.size.rep:2+2*self.size.rep-1);
                    self.data.abs.pfm(:,1:self.size.rep,i,j) = Data(:,2+2*self.size.rep:end)-self.params.PFCBgVoltage_mV_*1e-3;
                    % TotalFlDatas(:,1:self.size.rep,i,j) = lowpass(Data(:,2*((self.size.rep+1)+1)+2:2:2*(2*(self.size.rep+1))),50e-4/(self.dt*1e-6), 1/(self.dt*1e-6)); % LPF 2 kHz
                end
            end

            close(wb)
        end %readspecdata

        function self = CalTTNorm(self, type, NormStartTime)

            switch type
                case 'abs'
                    for j = 1 : self.size.freq
                        for i = 1 : self.size.iter
                            for k = 1 : self.size.rep
                                % base(k, i, j) = mean(data.abs.raw(data.baselinerange,k,i,j));
                                self.data.abs.rawbase(k, i, j) = mean(self.data.abs.raw(NormStartTime*1e3/self.dt:19*1e3/self.dt,k,i,j));
                                self.data.abs.rawnorm(:,k,i,j) = 1-(self.data.abs.raw(:,k,i,j)/self.data.abs.rawbase(k,i,j));
                                self.data.abs.pfmbase(k, i, j) = mean(self.data.abs.pfm(NormStartTime*1e3/self.dt:19*1e3/self.dt,k,i,j));
                                self.data.abs.pfmnorm(:,k,i,j) = 1-(self.data.abs.pfm(:,k,i,j)/self.data.abs.pfmbase(k,i,j));
                            end
                        end
                    end
                    self.data.abs.norm = self.data.abs.rawnorm - self.data.abs.pfmnorm;

                case 'fl'
                    for j = 1 : self.size.freq
                        for i = 1 : self.size.iter
                            for k = 1 : self.size.rep
                                self.data.fl.base(k, i, j) = mean(self.data.fl.raw(NormStartTime*1e3/self.dt:end,k,i,j));
                                self.data.fl.norm(:,k, i,j) = (self.data.fl.raw(:,k, i,j)-self.data.fl.base(k,i,j)); % time trace 2D image를 그려보고 8 ms으로 정햇음...
                            end
                        end
                    end
                otherwise
                    error('try abs,or fl')
            end

        end %CalTTNorm

        function self = CalTTMean(self, type, IterArray)
            switch type
                case 'abs'
                    self.data.abs.tt = reshape(mean(self.data.abs.norm(:,:,IterArray,:), [2 3]),self.size.time, self.size.freq);

                case 'fl'
                    self.data.fl.tt = reshape(mean(self.data.fl.norm(:,:,IterArray,:), [2 3]),self.size.time, self.size.freq);
            end
        end %CalTTMean

        function self = CalTTSum(self, type, SumEndTime)
            switch type
                case 'abs'
                    for i = 1 : self.size.rep
                        for j = 1 : self.size.freq
                            for k = 1 : self.size.iter
                                self.data.abs.sum(j,i,k) = sum(self.data.abs.norm(self.baselineidx+1:self.baselineidx+round(SumEndTime*1e3/self.dt),i,k,j));
                            end
                        end
                    end
                case 'fl'
                    for i = 1 : self.size.rep
                        for j = 1 : self.size.freq
                            for k = 1 : self.size.iter

                                self.data.fl.sum(j,i,k) = sum(self.data.fl.norm(self.baselineidx+1:self.baselineidx+round(SumEndTime*1e3/self.dt),i,k,j));
                            end
                        end
                    end

            end
        end %CalTTSum

        function self = CalSpecMeanSte(self, type, IterArray)
            switch type
                case 'abs'
                    self.data.abs.mean = squeeze(mean(self.data.abs.sum(:,:,IterArray),[2 3])); % results.repititionPerStep, results.iteration 전부 평균
                    self.data.abs.ste = std(self.data.abs.sum(:,:,IterArray),0,[2 3])/sqrt(length(IterArray)*self.size.rep);

                    self.data.abs.mean_iter = squeeze(mean(self.data.abs.sum,2)); % results.repititionPerStep 평균
                    self.data.abs.ste_iter = std(self.data.abs.sum,0,2)/sqrt(self.size.rep);

                case 'fl'
                    self.data.fl.mean = squeeze(mean(self.data.fl.sum(:,:,IterArray),[2 3]));
                    self.data.fl.ste = std(self.data.fl.sum(:,:,IterArray),0,[2 3])/sqrt(length(IterArray)*self.size.rep);

                    self.data.fl.mean_iter = squeeze(mean(self.data.fl.sum,2));
                    self.data.fl.ste_iter = squeeze(std(self.data.fl.sum,0,2)/sqrt(self.size.rep));
            end

        end %CalSpecMeanSte

        function self = CalSpecData(self, type, IterArray, NormStartTime, SumEndTime)
            switch type
                case 'abs'
                    CalTTNorm(self,'abs',NormStartTime);
                    CalTTMean(self,'abs', IterArray);
                    CalTTSum(self,'abs',SumEndTime);
                    CalSpecMeanSte(self,'abs',IterArray);
                case 'fl'
                    CalTTNorm(self,'fl',NormStartTime);
                    CalTTMean(self,'fl',IterArray);
                    CalTTSum(self,'fl',SumEndTime);
                    CalSpecMeanSte(self,'fl',IterArray);
            end
        end %CalSpecData

        % Plot methods
        function result = PlotSpectrum(self, ax, kind, band)

            switch band
                case 'UV'
                    det = self.det.UV.wm.mean;

                case 'IR'
                    % det = data.IRrealdetM;
                    det = self.det.IR.wm.mean;

                otherwise
                    error('try UV or IR');
            end


            switch kind
                case 'abs'
                    result = errorbar(ax, det.Q, self.data.abs.mean, self.data.abs.ste);
                    ylabel("Abs. signal (a.u.)");
                    title("In-cell Abs Spectrum");
                    xlabel("detuning from Q_{12}(1) (MHz)");

                case 'fl1'
                    if self.beamaxis == 'x'
                        result = errorbar(ax, det.Q, self.data.fl.mean, self.data.fl.ste);
                        ylabel("LIF. signal (a.u.)");
                        title("LIF Spectrum(x)");
                        xlabel("detuning from Q_{12}(1) (MHz)");

                    elseif self.beamaxis == 'xz'
                        result = errorbar(ax, det.P, self.data.fl.mean, self.data.fl.ste);
                        ylabel("Fl. signal (a.u.)");
                        title("LIF spectrum(xz)");
                        xlabel("detuning from P_{1}(1) (MHz)");
                    end

                    % case 'fl2'
                    %     %         result = errorbar(ax, data.Det, data.FM2, data.FSte2, 'LineWidth',1);
                    %     result = errorbar(ax, data.v, data.fl.mean2, data.fl.ste2);
                    %
                    %     % title("fl spectrum-2");
                    %     xlabel("velocity(P_1(1) ref) (m/s)",'fontsize',16);
                    %

                otherwise
                    error('try abs or fl1 or fl2');
            end
            result.LineWidth=1;
        end %plotspec

        function result = PlotTimetrace(self, ax, kind, freqidx)
            switch kind
                case 'abs'
                    result = plot(ax, self.t*1e-3, self.data.abs.tt(:,freqidx));
                    % xlim([0.08 self.t(end)*1e-3])
                    ylabel("Abs. signal (a.u.)")


                case 'fl'
                    result = plot(ax,self.t*1e-3, self.data.fl.tt(:,freqidx));
                    % xlim([0.1 self.t(end)*1e-3])
                    ylabel("Fl. signal (V)")

                otherwise
                    error('try abs or fl');
            end

            xlabel("time (ms)");
            result.LineWidth=1;
            xlim([0.04 18])
        end %PlotTimetrace

        function result = Plot2Dmap(self, ax, kind)
            StartTime = 0.04;
            EndTime = 18;

            idx = self.baselineidx+StartTime*1e3/self.dt;
            jdx = self.baselineidx+EndTime*1e3/self.dt;

            if self.size.freq == 1
                switch kind
                    case 'abs'
                        dd = reshape(self.data.abs.norm,self.size.time,self.size.rep*self.size.iter);
                    case 'fl1'
                        dd = reshape(self.data.fl.norm,self.size.time,self.size.rep*self.size.iter);
                    % case 'fl2'
                    %     dd = reshape(self.data.fl.norm2,self.size.time,self.size.rep*self.size.iter);
                    otherwise
                        error('try abs or fl')
                end

                result = imagesc(ax, self.t(idx:jdx)*1e-3, 1:self.size.rep*self.size.iter, dd(idx:jdx,:).');
                set(ax, 'Ydir','normal');
                ylabel("rep #",'fontsize',16);
                xlabel("time (ms)", 'FontSize',16);

                colorbar

            else

                switch kind
                    case 'abs'
                        result = imagesc(ax, self.t(idx:jdx)*1e-3, self.det.UV.wm.mean.Q, self.data.abs.tt(idx:jdx,:).');
                        ylabel("detuning from Q_{12}(1) (MHz)",'fontsize',16);
                    case 'fl1'
                        if self.beamaxis =='xz'
                            % result = imagesc(ax, self.t*1e-3, data.v, self.data.fl.tt.'); hold on;
                            result = imagesc(ax, self.t(idx:jdx)*1e-3, self.det.UV.wm.mean.P, self.data.fl.tt(idx:jdx,:).'); hold off;
                            % hold on;
                            % plot(self.t(self.baselindidx:end)*1e-3,data.maxv+2*data.vshift, '-.w', 'LineWidth',0.8); hold on;
                            % plot(self.t(self.baselindidx:end)*1e-3,data.maxv+data.vshift, '-.w', 'LineWidth',0.8); hold on;
                            % plot(self.t*1e-3,data.maxv, '-.w', 'LineWidth',0.8); hold off;
                            % plot(self.t(self.baselindidx:end)*1e-3,data.maxv-data.vshift, '-.w', 'LineWidth',0.8); hold off;
                            % plot(self.t(self.baselindidx:end)*1e-3,data.maxv-2*data.vshift, '-.w', 'LineWidth',0.8); hold off;
                            % ylabel("velocity (m/s)",'fontsize',16);

                            ylabel("detuning from P_{1}(1) (MHz)",'fontsize',16);

                        elseif self.beamaxis == 'x'
                            result = imagesc(ax, self.t(idx:jdx)*1e-3, self.det.UV.wm.mean.Q, self.data.fl.tt(idx:jdx,:).'); hold off;
                            ylabel("detuning from Q_{12}(1) (MHz)",'fontsize',16);
                        end
                    case 'fl2'
                        % result = imagesc(ax, self.t*1e-3, data.v, self.data.fl.tt2.'); hold on;
                        % plot(self.t(self.baselindidx:end)*1e-3,data.maxv+2*data.vshift, '-.w', 'LineWidth',0.8); hold on;
                        % plot(self.t(self.baselindidx:end)*1e-3,data.maxv+data.vshift, '-.w', 'LineWidth',0.8); hold on;
                        % plot(self.t(self.baselindidx:end)*1e-3,data.maxv, '-.w', 'LineWidth',0.8); hold on;
                        % plot(self.t(self.baselindidx:end)*1e-3,data.maxv-data.vshift, '-.w', 'LineWidth',0.8); hold on;
                        % plot(self.t(self.baselindidx:end)*1e-3,data.maxv-2*data.vshift, '-.w', 'LineWidth',0.8); hold off;

                        % ylabel("detuning from Q_{12}(1) (MHz)",'fontsize',16);
                        ylabel("velocity (m/s)",'fontsize',16);
                    otherwise
                        error('try abs or fl')
                end
                set(ax, 'Ydir','normal');
                xlabel("time (ms)", 'FontSize',16);
            end



            colorbar
        end %Plot2Dmap

        function result = PlotFigureSet(self,band)
            clf;
            if self.size.freq ~= 1

                figure('Name',self.runnum+" abs spec")
                PlotSpectrum(self, gca, 'abs', band);
                SaveFigToFile_v2(gcf,"K_results",self.runnum+"_abs_spectrum")
                SaveFigToFile_v2(gcf,"K_results",self.runnum+"_abs_spectrum")

                figure('Name',self.runnum+" lif spec")
                PlotSpectrum(self, gca, 'fl1', band);
                SaveFigToFile_v2(gcf,"K_results",self.runnum+"_LIF_spectrum");

                % if self.beamaxis == 'xz'
                %     figure('Name',self.runnum+" lif spec xz")
                % PlotSpectrum(self, gca, 'fl2', band);
                %     SaveFigToFile_v2(gcf,"K_results",self.runnum+"_LIF_spectrum-xz");
                % end

            elseif self.size.freq == 1
                figure('Name',self.runnum+" abs tt");
                PlotTimetrace(self, gca,'abs',1);
                SaveFigToFile_v2(gcf, "K_results",self.runnum+"_abs_tt");

                figure('Name',self.runnum+" lif tt");
                PlotTimetrace(self, gca,'fl',1);
                SaveFigToFile_v2(gcf, "K_results",self.runnum+"_lif_tt");
            end

            figure('Name',self.runnum+" abs 2D");
            Plot2Dmap(self, gca, 'abs');
            SaveFigToFile_v2(gcf, "K_results",self.runnum+"_abs_2D");

            figure('Name',self.runnum+" fl 2D");
            Plot2Dmap(self, gca, 'fl1');
            SaveFigToFile_v2(gcf, "K_results",self.runnum+"_lif_2D");

            % if data.beamaxis=='xy'
            %     figure('Name',self.runnum+" fl 2D")
            %     plot2Dtimedet_v3(gca, data, 'fl2');
            %     saveas(gcf, './K_results/'+self.runnum+'_fl_2D_2.png');
            %     saveas(gcf, './K_results/'+self.runnum+'_fl_2D_2.fig');
            % end

        end %PlotFigureSet

    end %method
end %classdef

