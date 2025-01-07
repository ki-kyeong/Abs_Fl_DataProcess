classdef MgFSpecData_v_laserslowing_2
    % LabView Sequence Program을 위한 대대적인 변경! 24/7/26
    % slowing data processing 추가, det/v calculation 수정 24/9/11
    % freq ref shift 기능 추가 24/10/2
    % AOM freq 기능 추가 24/10/언젠가..
    % data merge 기능 추가 24/10/18
    % 87Rb MTS freq data 추가 24/12/31
    % iteration random sorting 추가 25/01/07
    
    properties
        ExpType
        name
        runnum
        plotmode
        readme
        params = struct
        datasize = struct
        freq = struct
        det = struct
        dt
        data = struct
        baselineidx
        t
        beamaxis
        v
        maxv
        mint
        FreqRefArray = [417.147178, 417.1472425, 417.147300;
            546.742645, 0, 0;
            814.044490, 0, 0];
        LaserArray = ["Tet", "Puyo", "Zelda", 'Yeager'];
        TransitionArray = ["X(0,1)-A(0,0), P_{1}(1) F=2,1^{+}","X(0,1)-A(0,0), Q_{12}(1) F=0","X(0,1)-A(0,0), P_{1}(1) F=1^{-}";
            "X(2,1)-B(0,0), pP_{12}+pQ_{12}", "empty", "empty"];
        TransitionArray2 = ["P_{1}(1) F=2,1^{+}","Q_{12}(1) F=0","P_{1}(1) F=1^{-}";
            "pP_{12}+pQ_{12}", "empty", "empty"];
        frefshift
        scanlaser
        randidx=1;
        AOMshift

    end

    methods
        function self = MgFSpecData_v_laserslowing_2(name, plotmode, varargin)
            if nargin > 0
                p = inputParser;
                addParameter(p, 'ExpType', "MgF24")
                addParameter(p, 'beamaxis', "xz")
                addParameter(p, 'frefshift', 0)
                addParameter(p, 'scanlaser', 1)
                addParameter(p, 'AOMshift', 80)

                parse(p , varargin{:});
                self.ExpType = p.Results.ExpType;
                self.beamaxis = p.Results.beamaxis;
                self.frefshift = p.Results.frefshift;
                self.scanlaser = p.Results.scanlaser;
                self.AOMshift = p.Results.AOMshift;
                self = Main_Procedure(self, name, plotmode);
            end
        end %MgFSpecData

        function self = RaedReadme(self)
            self.readme = readlines(self.name+"readme.txt", Encoding="UTF-8");
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
            self.datasize.iter = self.params.iteration;
            self.datasize.rep = self.params.repitition;
        end %readparams

        function self = ReadFreqData(self) % 그냥 전부다 읽고나서 randidx찾는 것으로 바꿔야겠음..
            % freq.{IR, UV}.raw[freqscan,rep,iter,{1:Tet, 2:Puyo, 3:Zelda, 4:Yeager, 5: MTS}]
            for i = 1 : self.datasize.iter
                self.freq.IR.raw(:,:,i,1) = readmatrix(self.name+string(i-1)+"_wavemeter_data_TET IR.csv",Range=[2,3]);

                self.freq.IR.raw(:,:,i,2) = readmatrix(self.name+string(i-1)+"_wavemeter_data_PUYO IR.csv",Range=[2,3]);
                self.freq.UV.raw(:,:,i,2) = readmatrix(self.name+string(i-1)+"_wavemeter_data_PUYO UV.csv",Range=[2,3]);

                self.freq.IR.raw(:,:,i,3) = readmatrix(self.name+string(i-1)+"_wavemeter_data_Zelda VIS.csv",Range=[2,3]);

                self.freq.UV.raw(:,:,i,4) = readmatrix(self.name+string(i-1)+"_wavemeter_data_Yeager.csv",Range=[2,3]); 

                self.freq.IR.raw(:,:,i,5) = readmatrix(self.name+string(i-1)+"_wavemeter_data_MTS.csv",Range=[2,3]);               
            end
            for i = [1 3]
                self.freq.UV.raw(:,:,:,i) = 2*self.freq.IR.raw(:,:,:,i);
            end
            self.datasize.freq = size(self.freq.IR.raw,1);
            setfreqrange = [2 2 self.datasize.freq+1 2];

            % freq.{IR, UV}.set[freqscan,{1:Tet, 2:Puyo, 3:Zelda, 4:Yeager}]
            for i = 1 : self.datasize.iter
                self.freq.IR.set(:,i,1) = readmatrix(self.name+string(i-1)+"_wavemeter_data_TET IR.csv",Range=setfreqrange);
                self.freq.IR.set(:,i,2) = readmatrix(self.name+string(i-1)+"_wavemeter_data_PUYO IR.csv",Range=setfreqrange);
                self.freq.IR.set(:,i,3) = readmatrix(self.name+string(i-1)+"_wavemeter_data_Zelda VIS.csv",Range=setfreqrange);
                self.freq.UV.set(:,i,4) = readmatrix(self.name+string(i-1)+"_wavemeter_data_Yeager.csv",Range=setfreqrange);
                self.freq.IR.set(:,i,5) = readmatrix(self.name+string(i-1)+"_wavemeter_data_MTS.csv",Range=setfreqrange);
            end

            for i = [1 2 3]
                self.freq.UV.set(:,:,i) = 2*self.freq.IR.set(:,:,i);
            end

            self.freq.MTSshift = self.freq.IR.raw(:,:,:,5) - mean(self.freq.IR.set(:,:,5),'all');

            for i = 1 : self.datasize.iter
                for j = 1:4
                    self.freq.IR.raw(:,:,i,j) = self.freq.IR.raw(:,:,i,j)-self.frefshift*1e-6-self.freq.MTSshift(:,:,i);
                    self.freq.UV.raw(:,:,i,j) = self.freq.UV.raw(:,:,i,j)-self.frefshift*1e-6-self.freq.MTSshift(:,:,i)*2;
                end
            end

            self.freq.IR.set = self.freq.IR.set-2*self.frefshift*1e-6;
            self.freq.UV.set = self.freq.UV.set-2*self.frefshift*1e-6;

            self.freq.IR.raw(:,:,:,2) = self.freq.IR.raw(:,:,:,2)+0.5*self.AOMshift*1e-6;
            self.freq.IR.set(:,2) = self.freq.IR.set(:,2)+0.5*self.AOMshift*1e-6;
            self.freq.UV.raw(:,:,:,2) = self.freq.UV.raw(:,:,:,2)+self.AOMshift*1e-6;
            self.freq.UV.set(:,2) = self.freq.UV.set(:,2)+self.AOMshift*1e-6;
                      

            % random array sorting
            [~, self.randidx] = sort(self.freq.UV.set(:,:,self.scanlaser),1);
            % [~, self.randidx] = sort(self.freq.IR.set(:,self.scanlaser));
            
            for i = 1 : self.datasize.iter
                self.freq.IR.set(:,i,:) = self.freq.IR.set(self.randidx(:,i),i,:);
                self.freq.IR.raw(:,:,i,:) = self.freq.IR.raw(self.randidx(:,i),:,i,:);

                self.freq.UV.set(:,i,:) = self.freq.UV.set(self.randidx(:,i),i,:);
                self.freq.UV.raw(:,:,i,:) = self.freq.UV.raw(self.randidx(:,i),:,i,:);
            end

            % freq cal for Tet and Zelda
            for i = [1 2 3 4]
                if i ~=4
                self.freq.IR.mean(:,i) = mean(self.freq.IR.raw(:,:,:,i),[2 3]);
                self.freq.IR.ste(:,i) = std(self.freq.IR.raw(:,:,:,i),0,[2 3])/sqrt(self.datasize.iter*self.datasize.rep);
                end
                self.freq.UV.mean(:,i) = mean(self.freq.UV.raw(:,:,:,i),[2 3]);
                self.freq.UV.ste(:,i) = std(self.freq.UV.raw(:,:,:,i),0,[2 3])/sqrt(self.datasize.iter*self.datasize.rep);
            end
        end %readfreqdata

        function self = CalDetuning(self)
            % IR cavity detuning cal, det.{IR, UV}.set(freqscan,transition, {Tet, Puyo}]
            for i = 1:2
                for j = 1:3
                    self.det.IR.set(:,j,i) = (self.freq.IR.set(:,i)-self.FreqRefArray(1,j))*1e6; % MHz
                    self.det.UV.set(:,j,i) = 2*self.det.IR.set(:,j,i);

                    self.det.IR.raw(:,:,:,j,i) = (self.freq.IR.raw(:,:,:,i)-self.FreqRefArray(1,j))*1e6; % MHz
                    self.det.IR.mean(:,j,i) = mean(self.det.IR.raw(:,:,:,j,i),[2 3]);
                    self.det.IR.ste(:,j,i) = std(self.det.IR.raw(:,:,:,j,i),0,[2 3])/sqrt(self.datasize.iter*self.datasize.rep);
                end
            end
            % Detuning Cal for Zelda (X(1,1)-B(0,0) band)

            i=3;
            j=1;
            self.det.IR.set(:,j,i) = (self.freq.IR.set(:,i)-self.FreqRefArray(2,j))*1e6; % MHz
            self.det.UV.set(:,j,i) = 2*self.det.IR.set(:,j,i);

            self.det.IR.raw(:,:,:,j,i) = (self.freq.IR.raw(:,:,:,i)-self.FreqRefArray(2,j))*1e6; % MHz
            self.det.IR.mean(:,j,i) = mean(self.det.IR.raw(:,:,:,j,i),[2 3]);
            self.det.IR.ste(:,j,i) = std(self.det.IR.raw(:,:,:,j,i),0,[2 3])/sqrt(self.datasize.iter*self.datasize.rep);


            self.det.UV.raw = 2*self.det.IR.raw;
            self.det.UV.mean = 2*self.det.IR.mean;
            self.det.UV.ste = 2*self.det.IR.ste;

            % Detuning Cal for Yeager
            i=4;
            j=1;
            self.det.UV.set(:,j,i) = (self.freq.UV.set(:,i)-self.FreqRefArray(3,j))*1e6; % MHz
            
            self.det.UV.raw(:,:,:,j,i) = (self.freq.UV.raw(:,:,:,i)-self.FreqRefArray(3,j))*1e6; % MHz
            self.det.UV.mean(:,j,i) = mean(self.det.UV.raw(:,:,:,j,i),[2 3]);
            self.det.UV.ste(:,j,i) = std(self.det.UV.raw(:,:,:,j,i),0,[2 3])/sqrt(self.datasize.iter*self.datasize.rep);
        end %caldet

        function self = InitDataStruct(self)
            self.params.DataParams = detectImportOptions(self.name+"0_0.csv");
            self.params.DataParams.DataLines = [2 Inf];
            Data = readmatrix(self.name+"0_0.csv",self.params.DataParams);
            self.datasize.time = size(Data,1);
            self.data.abs.raw = zeros(self.datasize.time, self.datasize.rep, self.datasize.iter, self.datasize.freq); % abs time signal row data - bacground
            self.data.abs.pfm = zeros(self.datasize.time, self.datasize.rep, self.datasize.iter, self.datasize.freq); % abs time signal row data - bacground
            self.data.fl.raw = zeros(self.datasize.time,self.datasize.rep, self.datasize.iter, self.datasize.freq); % fls time signal row data
            % self.baselineidx = round(40/self.params.TimeDiv_ms_*1e-6); % 40 us에서 ablation peak끝
            % self.t = (Data(:,1)-Data(self.baselineidx,1)); % ms, 40 µs을 0초로 설정
            self.baselineidx = 1;
            self.t = Data(:,1)*1e3; % ms
            self.dt = self.params.TimeDiv_ms_; % ms
        end %initdatastruct

        function self = ReadSpecData(self)
            wb = waitbar(0, ' Getting started');
            for j = 1 : self.datasize.freq
                waitbar(j/self.datasize.freq, wb, self.runnum+newline+...
                    num2str(j/self.datasize.freq*100, '%.1f')+" % done",...
                    'WindowStyle','modal');
                for i = 1:self.datasize.iter
                    jdx = self.randidx(j,i);
                    Data = readmatrix(self.name+string(i-1)+"_"+string(jdx-1)+".csv",self.params.DataParams);
                    self.data.abs.raw(:,1:self.datasize.rep,i,j) = Data(:,2:2+self.datasize.rep-1)-self.params.AbsBgVoltage_mV_*1e-3;
                    self.data.fl.raw(:,1:self.datasize.rep,i,j) = Data(:,2+self.datasize.rep:2+2*self.datasize.rep-1);
                    self.data.abs.pfm(:,1:self.datasize.rep,i,j) = Data(:,2+2*self.datasize.rep:end)-self.params.PFCBgVoltage_mV_*1e-3;
                    % TotalFlDatas(:,1:self.size.rep,i,j) = lowpass(Data(:,2*((self.size.rep+1)+1)+2:2:2*(2*(self.size.rep+1))),50e-4/(self.dt*1e-6), 1/(self.dt*1e-6)); % LPF 2 kHz
                end
            end
            close(wb)
        end %readspecdata

        function self = CalTTNorm(self, type, varargin)

            p = inputParser;
            p.CaseSensitive=false;
            addParameter(p, 'absNormStartTime', 18)
            addParameter(p, 'flNormStartTime', 18)
            
            parse(p , varargin{:});

            
            absNormStartTime = p.Results.absNormStartTime;
            flNormStartTime = p.Results.flNormStartTime;

            switch type
                case 'abs'
                    for j = 1 : self.datasize.freq
                        for i = 1 : self.datasize.iter
                            for k = 1 : self.datasize.rep
                                % base(k, i, j) = mean(data.abs.raw(data.baselinerange,k,i,j));
                                self.data.abs.rawbase(k, i, j) = mean(self.data.abs.raw(absNormStartTime/self.dt:19/self.dt,k,i,j));
                                self.data.abs.rawnorm(:,k,i,j) = 1-(self.data.abs.raw(:,k,i,j)/self.data.abs.rawbase(k,i,j));
                                self.data.abs.pfmbase(k, i, j) = mean(self.data.abs.pfm(absNormStartTime/self.dt:19/self.dt,k,i,j));
                                self.data.abs.pfmnorm(:,k,i,j) = 1-(self.data.abs.pfm(:,k,i,j)/self.data.abs.pfmbase(k,i,j));
                            end
                        end
                    end
                    self.data.abs.norm = self.data.abs.rawnorm - self.data.abs.pfmnorm;

                case 'fl'
                    for j = 1 : self.datasize.freq
                        for i = 1 : self.datasize.iter
                            for k = 1 : self.datasize.rep
                                self.data.fl.base(k, i, j) = mean(self.data.fl.raw(flNormStartTime/self.dt:end,k,i,j));
                                self.data.fl.norm(:,k, i,j) = (self.data.fl.raw(:,k, i,j)-self.data.fl.base(k,i,j)); % time trace 2D image를 그려보고 8 ms으로 정햇음...
                            end
                        end
                    end
                otherwise
                    error('try abs,or fl')
            end
        end

        function self = CalTTMean(self, type, varargin)
            p = inputParser;
            p.CaseSensitive=false;
            addParameter(p, 'IterArray', 1 : self.datasize.iter)

            parse(p , varargin{:});
            IterArray = p.Results.IterArray;

            switch type
                case 'abs'
                    self.data.abs.tt = reshape(mean(self.data.abs.norm(:,:,IterArray,:), [2 3]),self.datasize.time, self.datasize.freq);

                case 'fl'
                    self.data.fl.tt = reshape(mean(self.data.fl.norm(:,:,IterArray,:), [2 3]),self.datasize.time, self.datasize.freq);
            end
        end %CalTTMean

        function self = CalTTSum(self, type,varargin)
            p = inputParser;
            p.CaseSensitive=false;
           
            addParameter(p, 'absSumStartTime', 0.2)
            addParameter(p, 'flSumStartTime', 0.2)
            addParameter(p, 'SumEndTime', 18)

            parse(p , varargin{:});

            SumEndTime = p.Results.SumEndTime;
            absSumStartTime = p.Results.absSumStartTime;
            flSumStartTime = p.Results.flSumStartTime;
            
            switch type
                case 'abs'
                    for i = 1 : self.datasize.rep
                        for j = 1 : self.datasize.freq
                            for k = 1 : self.datasize.iter
                                self.data.abs.sum(j,i,k) = sum(self.data.abs.norm(self.baselineidx+round(absSumStartTime/self.dt):self.baselineidx+round(SumEndTime/self.dt),i,k,j));
                            end
                        end
                    end
                case 'fl'
                    for i = 1 : self.datasize.rep
                        for j = 1 : self.datasize.freq
                            for k = 1 : self.datasize.iter
                                self.data.fl.sum(j,i,k) = sum(self.data.fl.norm(self.baselineidx+round(flSumStartTime/self.dt):self.baselineidx+round(SumEndTime/self.dt),i,k,j)); % slowing laser PMT saturation 이후부터 더해야함..
                            end
                        end
                    end

            end
        end %CalTTSum

        function self = CalSpecMeanSte(self, type, varargin)
            p = inputParser;
            p.CaseSensitive=false;
            addParameter(p, 'IterArray', 1 : self.datasize.iter)

            parse(p , varargin{:});
            IterArray = p.Results.IterArray;

            switch type
                case 'abs'
                    self.data.abs.mean = squeeze(mean(self.data.abs.sum(:,:,IterArray),[2 3])); % results.repititionPerStep, results.iteration 전부 평균
                    self.data.abs.ste = std(self.data.abs.sum(:,:,IterArray),0,[2 3])/sqrt(length(IterArray)*self.datasize.rep);

                    self.data.abs.mean_iter = squeeze(mean(self.data.abs.sum,2)); % results.repititionPerStep 평균
                    self.data.abs.ste_iter = std(self.data.abs.sum,0,2)/sqrt(self.datasize.rep);

                case 'fl'
                    self.data.fl.mean = squeeze(mean(self.data.fl.sum(:,:,IterArray),[2 3]));
                    self.data.fl.ste = std(self.data.fl.sum(:,:,IterArray),0,[2 3])/sqrt(length(IterArray)*self.datasize.rep);

                    self.data.fl.mean_iter = squeeze(mean(self.data.fl.sum,2));
                    self.data.fl.ste_iter = squeeze(std(self.data.fl.sum,0,2)/sqrt(self.datasize.rep));
            end

        end %CalSpecMeanSte

        function self = CalSpecData_total(self, type, varargin)
            p = inputParser;
            p.CaseSensitive=false;
            addParameter(p, 'IterArray', 1 : self.datasize.iter)
            addParameter(p, 'absNormStartTime', 18)
            addParameter(p, 'flNormStartTime', 18)
            addParameter(p, 'absSumStartTime', 0.2)
            addParameter(p, 'flSumStartTime', 0.2)
            addParameter(p, 'SumEndTime', 18)

            parse(p , varargin{:});

            IterArray = p.Results.IterArray;
            absNormStartTime = p.Results.absNormStartTime;
            flNormStartTime = p.Results.flNormStartTime;
            SumEndTime = p.Results.SumEndTime;
            absSumStartTime = p.Results.absSumStartTime;
            flSumStartTime = p.Results.flSumStartTime;

            switch type
                case 'abs'
                    self = CalTTNorm(self,'abs',absNormStartTime=absNormStartTime);
                    self = CalTTMean(self,'abs', IterArray=IterArray);
                    self = CalTTSum(self,'abs',absSumStartTime=absSumStartTime,SumEndTime=SumEndTime);
                    self = CalSpecMeanSte(self,'abs',IterArray=IterArray);
                case 'fl'
                    self = CalTTNorm(self,'fl',flNormStartTime=flNormStartTime);
                    self = CalTTMean(self,'fl',IterArray=IterArray);
                    self = CalTTSum(self,'fl',flSumStartTime=flSumStartTime,SumEndTime=SumEndTime);
                    self = CalSpecMeanSte(self,'fl',IterArray=IterArray);
            end
        end %CalSpecData

        function self = CalSpecData_merge(self, type, varargin)
            p = inputParser;
            addParameter(p, 'IterArray', 1 : self.datasize.iter)
            addParameter(p, 'absSumStartTime', 0.2)
            addParameter(p, 'flSumStartTime', 3)
            addParameter(p, 'SumEndTime', 18)

            parse(p , varargin{:});

            IterArray = p.Results.IterArray;
            
            SumEndTime = p.Results.SumEndTime;
            absSumStartTime = p.Results.absSumStartTime;
            flSumStartTime = p.Results.flSumStartTime;

            switch type
                case 'abs'
                    self = CalTTMean(self,'abs', IterArray=IterArray);
                    self = CalTTSum(self,'abs',absSumStartTime=absSumStartTime,SumEndTime=SumEndTime);
                    self = CalSpecMeanSte(self,'abs',IterArray=IterArray);
                case 'fl'
                    self = CalTTMean(self,'fl',IterArray=IterArray);
                    self = CalTTSum(self,'fl',flSumStartTime=flSumStartTime,SumEndTime=SumEndTime);
                    self = CalSpecMeanSte(self,'fl',IterArray=IterArray);
            end
        end %CalSpecData

        function self = CalVelocity(self, varargin)
            p = inputParser;
            addParameter(p, 'BeamDistance', 0.70);
            addParameter(p, 'BeamAngle', 45*pi/180);

            parse(p, varargin{:});

            BeamDistance = p.Results.BeamDistance;
            BeamAngle = p.Results.BeamAngle;


            for FreqRef = 1:3
                self.v(:,FreqRef) = -self.det.IR.mean(:,FreqRef,1)/(self.FreqRefArray(1,FreqRef)*1e6)*299792458/cos(BeamAngle);
                self.maxv(:,FreqRef) = BeamDistance./(self.t(self.baselineidx:end)*1e-3); % possible maximum velocity
                self.mint(:,FreqRef) = BeamDistance./self.v(:,FreqRef)*1e3; % possible minimum arrival time, ms; % possible minimum arrival time, µs
            end

        end

        function self = CalTimeBinSpectra(self, varargin)
            p = inputParser;
            addParameter(p, 'StartTime', 3);
            addParameter(p, 'TimeBin', 0.5);
            addParameter(p, 'EndTime', 10);

            addParameter(p, 'a1', [1 1 1 1 1 1 1 1 1 1 1 1]);
            addParameter(p, 'c1', 400);

            parse(p, varargin{:});

            StartTime = p.Results.StartTime;
            TimeBin = p.Results.TimeBin;
            EndTime = p.Results.EndTime;

            a1 = p.Results.a1;
            c1 = p.Results.c1;

            StartTimeIdx = (StartTime+0.2)/self.dt;
            TimeBinIdx = TimeBin/self.dt;
            EndTimeIdx = EndTime/self.dt;
            TimeBinArrayIdx = (StartTimeIdx:TimeBinIdx/2:EndTimeIdx-TimeBinIdx/2);

            % for i = 3:3
            for i = 1: length(TimeBinArrayIdx)
                self.data.TimeBinSpectra(:,i) = sum(self.data.fl.tt(TimeBinArrayIdx(i):TimeBinArrayIdx(i)+TimeBinIdx-1,:),1);
                self.data.TimeBinSpectraSum(i) = sum(self.data.TimeBinSpectra(:,i));
                switch self.beamaxis
                    case 'xz'
                        if i == 1
                            bshiftlim = -1000;
                            binitguess = -400;
                        else
                            % bshiftlim = self.data.TimeBinSpectrafitresult{i-1}.b1;
                            % binitguess = -0.7/self.t(TimeBinArrayIdx(i))*1e-3/299792458*2*417147242*cos(pi/4);
                            binitguess = self.data.TimeBinSpectrafitresult{i-1}.b1;
                        end
                        self.data.TimeBinSpectrafitresult{i} = HFEOM_Lorentzian(self, self.det.UV.mean(:,2,1), self.data.TimeBinSpectra(:,i),...
                            [a1 binitguess c1 0], bshiftlim);
                    case 'x'
                        self.data.TimeBinSpectrafitresult{i} = HFEOM_Lorentzian(self, self.det.UV.mean(:,2,1), self.data.TimeBinSpectra(:,i),...
                            [a1 0 c1 0], 0);

                end
                title(string(self.t(TimeBinArrayIdx(i)))+" ms spectrum");
                SaveFigToFile_v2(gcf, "K_results",self.runnum+"fit result_"+string(i));
            end

            self.data.TimeBinSpectraIdx = TimeBinArrayIdx;

            for i = 1 : length(self.data.TimeBinSpectrafitresult)
                self.data.TimeBinVelocity(i)=-self.data.TimeBinSpectrafitresult{i}.b1/(2*417147242.5)*299792458/cos(pi/4);
            end

            % figure(12031230)
            % plot(self.data.TimeBinVelocity,self.t(TimeBinArrayIdx(1:(end-1))))
            %
            % self.data.TimeBinVelocityDist = sum(self.data.TimeBinSpectra,1);
            % figure(12301203)
            % plot(self.data.TimeBinVelocity, self.data.TimeBinVelocityDist);


        end

        function merge_Data = MergeRunData(merge_Data, Data ,merge_runs,varargin) % 만드는 중...
            p = inputParser;
            addParameter(p, 'absSumStartTime', 0.2)
            addParameter(p, 'flSumStartTime', 3)
            addParameter(p, 'SumEndTime', 18)

            parse(p , varargin{:});

            SumEndTime = p.Results.SumEndTime;
            absSumStartTime = p.Results.absSumStartTime;
            flSumStartTime = p.Results.flSumStartTime;
            
            merge_Data = Data(merge_runs(1));
            merge_length = length(merge_runs);
            
            for i = 2 : merge_length
                
            end

            for i = 2 : merge_length
                d = Data(merge_runs(i));
                diter = d.datasize.iter;
                % norm data들 합치기
                merge_Data.datasize.iter = merge_Data.datasize.iter + diter;
                merge_Data.data.abs.norm(:,:,end+1:end+diter,:)= d.data.abs.norm;
                merge_Data.data.fl.norm(:,:,end+1:end+diter,:)= d.data.fl.norm;

                % raw freq, det들 합치기
                merge_Data.freq.IR.raw(:,:,end+1:end+diter,:) = d.freq.IR.raw;
                merge_Data.freq.UV.raw(:,:,end+1:end+diter,:) = d.freq.UV.raw;
                merge_Data.det.IR.raw(:,:,end+1:end+diter,:,:) = d.det.IR.raw;
                merge_Data.det.UV.raw(:,:,end+1:end+diter,:,:) = d.det.UV.raw;
            end

            % 데이터 평균과 표준오차 구하기
            merge_Data = CalSpecData_merge(merge_Data,'abs',absSumStartTime=absSumStartTime, SumEndTime=SumEndTime);
            merge_Data = CalSpecData_merge(merge_Data,'fl',flSumStartTime=flSumStartTime, SumEndTime=SumEndTime);

            % 주파수도 다시 평균과 표준오차 구하기
            merge_Data.freq.IR.mean = squeeze(mean(merge_Data.freq.IR.raw, [2 3]));
            merge_Data.freq.IR.ste = squeeze(std(merge_Data.freq.IR.raw,0, [2 3])/sqrt(merge_Data.datasize.iter*merge_Data.datasize.rep));
            merge_Data.freq.UV.mean = squeeze(mean(merge_Data.freq.UV.raw, [2 3]));
            merge_Data.freq.UV.ste = squeeze(std(merge_Data.freq.UV.raw,0, [2 3])/sqrt(merge_Data.datasize.iter*merge_Data.datasize.rep));
            merge_Data.det.IR.mean = squeeze(mean(merge_Data.det.IR.raw, [2 3]));
            merge_Data.det.IR.ste = squeeze(std(merge_Data.det.IR.raw,0, [2 3])/sqrt(merge_Data.datasize.iter*merge_Data.datasize.rep));
            merge_Data.det.UV.mean = squeeze(mean(merge_Data.det.UV.raw, [2 3]));
            merge_Data.det.UV.ste = squeeze(std(merge_Data.det.UV.raw,0, [2 3])/sqrt(merge_Data.datasize.iter*merge_Data.datasize.rep));

        end

        function self = Main_Procedure(self, name, plotmode)
            clf;
            self.name = name;
            temp = split(name,'_');
            self.runnum = temp(end-1);
            self.plotmode = plotmode;

            self = RaedReadme(self);
            self = ReadExpParams(self);

            self = ReadFreqData(self);
            self = CalDetuning(self);

            self = InitDataStruct(self);
            self = ReadSpecData(self);
            
            self = CalSpecData_total(self,'abs');
            self = CalSpecData_total(self,'fl');

            if self.beamaxis == 'xz'
                self = CalVelocity(self);
            end

            % if self.plotmode == 1
            %     PlotFigureSet(self)
            % end

        end %MgF24_Procedure


        function [Mean, Ste] = Data_massage(~, data1, data2, type)
            switch type
                case "mul"
                    Mean = data1.mean.*data2.mean;
                    Ste = abs(Mean).*sqrt((data1.ste./data1.mean).^2+(data2.ste./data2.mean).^2);
                case "div"
                    Mean = data1.mean./data2.mean;
                    Ste = abs(Mean).*sqrt((data1.ste./data1.mean).^2+(data2.ste./data2.mean).^2);
                case "sum"
                    Mean = data1.mean+data2.mean;
                    Ste = sqrt((data1.ste).^2+(data2.ste).^2);
                case "sub"
                    Mean = data1.mean-data2.mean;
                    Ste = sqrt((data1.ste).^2+(data2.ste).^2);
            end
        end%Data_massage

        % fitting methods
        function result = SingleFreqX_Lorentzian(self, varargin)

            p = inputParser;
            addParameter(p, 'fitparams', [1 1 1 0 400 0]);
            % addParameter(p, 'fitparams', [1 1 1 1 0 400 0]);
            
            parse(p, varargin{:});

            fitparams = p.Results.fitparams;

            % hfs = [129.5 0 -109.8];
            % hfsst = ["+129.5" "+0" "-109.8"]; % hyperfine strunture detuning from F=0, [2 1+ 0 1-]
            hfs = [126 0 -109.8];
            hfsst = ["+126" "+0" "-109.8"]; % hyperfine strunture detuning from F=0, [2 1+ 0 1-]
            % hfs = [129.5 120.3 0 -109.7];
            % hfsst = ["+129.5" "+120.3" "+0" "-109.7"]; % hyperfine strunture detuning from F=0, [2 1+ 0 1-]
            gamma = 22;
            gamma_factor = 22^2/4;
            AAA = string(gamma_factor);
            % amparray = [5 3 1 3];
            % Amparray = string(amparray);

            [xData, yData] = prepareCurveData( self.det.UV.mean(:,2,1), self.data.fl.mean);

            % Set up fittype and options.
            FitFunction = "";
            for i = 1:3
                fitfunction_temp = "a"+string(i)+"*2*"+AAA+"/((x-b1"+hfsst(i)+")^2+("+AAA+")*(1+c1))+";
                % fitfunction_temp = Amparray(i)+"*a1*2*"+AAA+"/((x-b1"+hfsst(i)+")^2+("+AAA+")*(1+c1))+";
                FitFunction = FitFunction+fitfunction_temp;
            end
            FitFunction = FitFunction+"d1";
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.Lower = [0 0 0 -30  0 -Inf ];
            % opts.Lower = [0 0 0 0 -30  0 -Inf ];
            opts.StartPoint = fitparams;

            % Fit model to data.
            [result, gof] = fit( xData, yData, FitFunction, opts);

            figure1=figure( 'Name', 'L fit' );
            errorbar(gca, xData, yData, self.data.fl.ste, LineWidth=1.5); hold on;
            h = plot( result); hold on;
            h.LineWidth=2;
            for i = 1:3
                plot(xData, eval("result.a"+string(i))*2*(gamma^2/4)./((xData-result.b1+hfs(i)).^2+gamma^2/4*(1+result.c1))); hold on;
            end
            hold off
            xlabel("detuning from F=0 (MHz)")


            legend("Data", "\Delta = "+string(result.b1)+" MHz")
            SaveFigToFile_v2(gcf, "K_results",self.runnum+"_LF");

        end

        function result = HFEOM_Lorentzian(self, xData, yData, fitparams, bshiftlim)
            % hfs = [124.9 0 -109.7];
            % hfsst = ["+124.9" "+0" "-109.7"]; % hyperfine strunture detuning from F=0, [2 1+ 0 1-]
            % hfs = [129.5 120.3 0 -109.7];
            % hfsst = ["+129.5" "+120.3" "+0" "-109.7"]; % hyperfine strunture detuning from F=0, [2 1+ 0 1-]
            % ModFreq = [-115 0 115];
            % ModFreqst = ["+115","+0","-115"];

            hfs = [240 120.3 0 -115 -225];
            hfsst = ["+240" "+120.3" "+0" "-115" "-225"];
            gamma = 22;
            gamma_factor = gamma^2/4;
            AAA = string(gamma_factor);
            ModAmpst= [1 1 1];


            % Set up fittype and options.
            FitFunction = "";
            for i = 1:5
                
                    fitfunction_temp = "a"+string(i)+"*2*"+AAA+"/((x-b1"+hfsst(i)+")^2+("+AAA+")*(1+c"+string(i)+"))+";
                    FitFunction = FitFunction+fitfunction_temp;
                
            end
            % for i = 1:4
            %     for j = 1:3
            %         fitfunction_temp = "a"+string((i-1)*3+j)+"*2*"+AAA+"/((x-b1"+hfsst(i)+ModFreqst(j)+")^2+("+AAA+")*(1+"+ModAmpst(j)+"*c1))+";
            %         FitFunction = FitFunction+fitfunction_temp;
            %     end
            % end
            FitFunction = FitFunction+"d1";
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.Lower = [0 0 0 0 0 bshiftlim 0 0 0 0 0 0];
            % opts.Lower = [0 0 0 0 0 0 0 0 0 0 0 0 bshiftlim -Inf ];
            opts.StartPoint = fitparams;

            % Fit model to data.
            [result, gof] = fit( xData, yData, FitFunction, opts);

            figure( 'Name', 'L fit' );
            plot(gca, xData, yData); hold on;
            h = plot(result); hold on;
            for i = 1:5
                    plot(xData, eval("result.a"+string(i))*2*(gamma_factor)./((xData-result.b1+hfs(i)).^2+gamma_factor*(1+eval("result.c"+string(i))))); hold on;
            end
            % for i = 1:3
            %     for j = 1:3
            %         plot(xData, eval("result.a"+string(3*(i-1)+j))*2*(gamma_factor)./((xData-result.b1+hfs(i)+ModFreq(j)).^2+gamma_factor*(1+result.c1))); hold on;
            %     end
            % end
            hold off
            legend("Data", "v = "+string(-result.b1/(2*417147242)*299792458/cos(pi/4))+" m/s",...
                "-2", "-1","0","+1","+2");
            % SaveFigToFile_v2(gcf, "K_results",self.runnum+"_LF");

        end



        % Plot methods

        function result = PlotSpectrum(self, kind, varargin)
            p = inputParser;
            addParameter(p, 'ax', gca);
            addParameter(p, 'FreqRef', 2);
            addParameter(p, 'Laser', 1);
            addParameter(p, 'lw', 1);
            addParameter(p, 'fs', 18);

            parse(p, varargin{:});

            ax = p.Results.ax;
            FreqRef = p.Results.FreqRef;
            Laser = p.Results.Laser;
            lw = p.Results.lw;
            fs = p.Results.fs;

            detval = self.det.UV.mean(:,FreqRef, Laser);

            switch kind
                case 'abs'
                    result = errorbar(ax, detval, self.data.abs.mean, self.data.abs.ste);
                    ylabel("Abs. signal (a.u.)");
                    title("In-cell Abs Spectrum, "+self.LaserArray(2));
                    xlabel(self.LaserArray(Laser)+" detuning from "+self.TransitionArray(1, FreqRef)+"(MHz)");

                case 'fl'
                    if self.beamaxis == 'x'
                        result = errorbar(ax, detval, self.data.fl.mean, self.data.fl.ste);
                        ylabel("LIF. signal (a.u.)");
                        title("LIF Spectrum(x), "+self.LaserArray(Laser));
                        xlabel(self.LaserArray(Laser)+" detuning from "+self.TransitionArray(1, FreqRef)+"(MHz)");

                    elseif self.beamaxis == 'xz'
                        result = errorbar(ax, detval, self.data.fl.mean, self.data.fl.ste);
                        ylabel("Fl. signal (a.u.)");
                        title("LIF spectrum(xz), "+self.LaserArray(Laser));
                        xlabel(self.LaserArray(Laser)+" detuning from "+self.TransitionArray(1, FreqRef)+"(MHz)");
                    end



                otherwise
                    error('try abs or fl');
            end
            result.LineWidth=lw;
            result.Parent.FontSize=fs;


        end %plotspec

        function result = PlotTimetrace(self, kind, varargin)
            p = inputParser;
            addParameter(p, 'ax', gca);
            addParameter(p, 'freqidx', 1);
            addParameter(p, 'FreqRef', 2);
            addParameter(p, 'Laser', 1);
            addParameter(p, 'tlim', [0 15]);
            addParameter(p, 'LW', 1);
            addParameter(p, 'fs', 18);

            parse(p, varargin{:});

            ax = p.Results.ax;
            FreqRef = p.Results.FreqRef;
            Laser = p.Results.Laser;
            freqidx = p.Results.freqidx;
            tlim = p.Results.tlim;
            LW = p.Results.LW;
            fs = p.Results.fs;

            switch kind
                case 'abs'
                    result = plot(ax, self.t, self.data.abs.tt(:,freqidx));
                    xlim(tlim)
                    ylabel("Abs. signal (a.u.)")
                    title("In-cell Abs timetrace, "+self.LaserArray(Laser));


                case 'fl'
                    result = plot(ax,self.t, self.data.fl.tt(:,freqidx));
                    xlim(tlim)
                    ylabel("Fl. signal (V)")

                otherwise
                    error('try abs or fl');
            end

            xlabel("time (ms)");
            result.LineWidth=LW;
            result.Parent.FontSize=fs;

        end %PlotTimetrace

        function result = Plot2Dmap(self, varargin)
            p = inputParser;
            addParameter(p, 'StartTime', 0);
            addParameter(p, 'EndTime', 10);
            addParameter(p, 'SlowingEndTime', 3);
            addParameter(p, 'FreqRef', 2);
            addParameter(p, 'Laser', 1);
            addParameter(p, 'kind', 'fl');

            parse(p, varargin{:});

            StartTime = p.Results.StartTime;
            EndTime = p.Results.EndTime;
            SlowingEndTime = p.Results.SlowingEndTime;
            FreqRef = p.Results.FreqRef;
            Laser = p.Results.Laser;
            kind = p.Results.kind;
            
            if self.datasize.freq == 1
                Plot2Dmap_single(ax, self, SlowingEndTime=SlowingEndTime, kind=kind);
            else
                switch self.beamaxis
                    case 'x'
                        result = Plot2Dmap_x(self, StartTime=StartTime, EndTime=EndTime, SlowingEndTime=SlowingEndTime, ...
                            FreqRef=FreqRef, Laser=Laser, kind=kind);
                    case 'xz'
                        result = Plot2Dmap_xz(self, StartTime=StartTime, EndTime=EndTime, SlowingEndTime=SlowingEndTime, ...
                            FreqRef=FreqRef, Laser=Laser, kind=kind);
                end
            end

        end %Plot2Dmap

        function result = Plot2Dmap_xz(self, varargin)
            p = inputParser;
            addParameter(p, 'StartTime', 0);
            addParameter(p, 'EndTime', 10);
            addParameter(p, 'SlowingEndTime', 3);
            addParameter(p, 'FreqRef', 2);
            addParameter(p, 'Laser', 1);
            addParameter(p, 'kind', 'fl');
            parse(p, varargin{:});
            StartTime = p.Results.StartTime;
            EndTime = p.Results.EndTime;
            SlowingEndTime = p.Results.SlowingEndTime;
            FreqRef = p.Results.FreqRef;
            Laser = p.Results.Laser;
            kind = p.Results.kind;
            ax = gca;
            idx = self.baselineidx+StartTime/self.dt;
            jdx = self.baselineidx+EndTime/self.dt;
            switch kind
                case 'abs'
                    result = imagesc(ax, self.t(idx:jdx), self.det.UV.mean(:,FreqRef, Laser), self.data.abs.tt(idx:jdx,:).');
                    ylabel(self.LaserArray(Laser)+"detuning from "+newline+self.TransitionArray(1, FreqRef)+" (MHz)",'fontsize',12);
                case 'fl'
                    result = imagesc(ax, self.t(idx:jdx), self.v(:,FreqRef), self.data.fl.tt(idx:jdx,:).'); hold on;
                    ylabel("velocity (m/s)"+newline+self.LaserArray(Laser)+', '+self.TransitionArray2(1, FreqRef)+" (MHz))",'fontsize',16);
                    plot(self.t, self.maxv(:,FreqRef), LineStyle="--", LineWidth=2, Color="w"); hold off;
                    set(ax, 'Ydir','normal');
                otherwise
                    error('try abs or fl')
            end
            xlabel("time (ms)", 'FontSize',16);
            clim([min(self.data.fl.tt((SlowingEndTime+0.1)/self.dt:end,:),[],'all') max(self.data.fl.tt((SlowingEndTime+0.1)/self.dt:end,:),[],'all')])
            colorbar()
        end %Plot2Dmap_xz

        function result = Plot2Dmap_x(self, varargin)
            p = inputParser;
            addParameter(p, 'StartTime', 0);
            addParameter(p, 'EndTime', 10);
            addParameter(p, 'SlowingEndTime', 3);
            addParameter(p, 'FreqRef', 2);
            addParameter(p, 'Laser', 1);
            addParameter(p, 'kind', 'fl');
            parse(p, varargin{:});
            StartTime = p.Results.StartTime;
            EndTime = p.Results.EndTime;
            SlowingEndTime = p.Results.SlowingEndTime;
            FreqRef = p.Results.FreqRef;
            Laser = p.Results.Laser;
            kind = p.Results.kind;
            ax = gca;
            idx = self.baselineidx+StartTime/self.dt;
            jdx = self.baselineidx+EndTime/self.dt;
            switch kind
                case 'abs'
                    result = imagesc(ax, self.t(idx:jdx), self.det.UV.mean(:,FreqRef, Laser), self.data.abs.tt(idx:jdx,:).');
                    ylabel(self.LaserArray(Laser)+"detuning from "+newline+self.TransitionArray(1, FreqRef)+" (MHz)",'fontsize',12);
                case 'fl'
                    result = imagesc(ax, self.t(idx:jdx), -self.det.UV.mean(:,FreqRef, Laser)/2/417147242*299792458, self.data.fl.tt(idx:jdx,:).'); hold on;
                    ylabel("target center velocity (m/s)"+newline+"detuning from "+self.TransitionArray(1, FreqRef)+" (MHz)",'fontsize',16);
                    % result = imagesc(ax, self.t(idx:jdx), self.det.UV.mean(:,FreqRef, Laser), self.data.fl.tt(idx:jdx,:).'); hold on;
                    % ylabel(self.LaserArray(Laser)+" detuning from "+newline+self.TransitionArray(1, FreqRef)+" (MHz)",'fontsize',16);
                otherwise
                    error('try abs or fl')
            end
            xlabel("time (ms)", 'FontSize',16);
            % clim([min(self.data.fl.tt((SlowingEndTime+0.1)/self.dt:end,:),[],'all') max(self.data.fl.tt((SlowingEndTime+0.1)/self.dt:end,:),[],'all')]);
            colorbar()
            set(ax, 'Ydir','normal');
        end %Plot2Dmap_x

        function result = Plot2Dmap_single(self, varargin)
            p = inputParser;
            addParameter(p, 'SlowingEndTime', 3);
            addParameter(p, 'kind', 'fl');
            parse(p, varargin{:});
            SlowingEndTime = p.Results.SlowingEndTime;
            kind = p.Results.kind;
            ax = gca;
            switch kind
                case 'abs'
                    dd = reshape(self.data.abs.norm,self.datasize.time,self.datasize.rep*self.datasize.iter);
                case 'fl'
                    dd = reshape(self.data.fl.norm,self.datasize.time,self.datasize.rep*self.datasize.iter);
                    % case 'fl2'
                    %     dd = reshape(self.data.fl.norm2,self.size.time,self.size.rep*self.size.iter);
                otherwise
                    error('try abs or fl')
            end
            result = imagesc(ax, self.t(idx:jdx), 1:self.datasize.rep*self.datasize.iter, dd(idx:jdx,:).');
            set(ax, 'Ydir','normal');
            ylabel("rep #",'fontsize',16);
            xlabel("time (ms)", 'FontSize',16);
            clim([min(self.data.fl.tt((SlowingEndTime+0.1)/self.dt:end,:),[],'all') max(self.data.fl.tt((SlowingEndTime+0.1)/self.dt:end,:),[],'all')]);
            colorbar()
        end %Plot2Dmap_single



        function result = Plot2DmapSub(Data1, Data2, varargin)
            p = inputParser;
            addParameter(p, 'StartTime', 0);
            addParameter(p, 'EndTime', 10);
            addParameter(p, 'SlowingEndTime', 3);
            addParameter(p, 'FreqRef', 2);
            addParameter(p, 'Laser', 1);
            addParameter(p, 'kind', 'fl');

            parse(p, varargin{:});

            StartTime = p.Results.StartTime;
            EndTime = p.Results.EndTime;
            SlowingEndTime = p.Results.SlowingEndTime;
            FreqRef = p.Results.FreqRef;
            Laser = p.Results.Laser;
            kind = p.Results.kind;

            ax = gca;

            idx = Data1.baselineidx+StartTime/Data1.dt;
            jdx = Data1.baselineidx+EndTime/Data1.dt;

            plotdata = Data1.data.fl.tt - Data2.data.fl.tt;

            result = imagesc(ax, Data1.t(idx:jdx), -Data1.det.UV.mean(1:end-1,FreqRef, Laser)/2/417147242.5*299792458, plotdata(idx:jdx,1:end-1).'); hold on;
            % result = imagesc(ax, self.t(idx:jdx), -self.det.UV.mean(1:end-1,FreqRef, Laser)/2/546742645*299792458, self.data.fl.tt(idx:jdx,1:end-1).'); hold on;
            ylabel(Data1.LaserArray(Laser)+" detuning from "+Data1.TransitionArray(1, FreqRef)+"(MHz)",'fontsize',16);

            % result = imagesc(ax, Data1.t(idx:jdx), Data1.v(:,FreqRef), plotdata(idx:jdx,:).'); hold on;
            % ylabel("velocity (m/s, "+Data1.LaserArray(Laser)+", "+Data1.TransitionArray2(1, FreqRef)+"(MHz))",'fontsize',16);
            % result = imagesc(ax, Data1.t(idx:jdx), Data1.det.UV.mean(:,FreqRef, Laser), plotdata(idx:jdx,:).'); hold on;
            % ylabel(Data1.LaserArray(Laser)+" detuning from "+Data1.TransitionArray(1, FreqRef)+"(MHz)",'fontsize',16);
            xlabel("time (ms)", 'FontSize',16);
            colorbar()
            % plot(Data1.t, Data1.maxv(:,FreqRef), LineStyle="--", LineWidth=2, Color="w"); hold off;

            set(ax, 'Ydir','normal');

            clim([min(plotdata((SlowingEndTime+0.1)/Data1.dt:end,:),[],'all') max(plotdata((SlowingEndTime+0.1)/Data1.dt:end,:),[],'all')])

        end

        function result = Plot2DmapTotal(Data1, Data2, varargin)
            p = inputParser;
            addParameter(p, 'StartTime', 0);
            addParameter(p, 'EndTime', 10);
            addParameter(p, 'SlowingEndTime', 3);
            addParameter(p, 'FreqRef', 2);
            addParameter(p, 'Laser', 1);
            addParameter(p, 'kind', 'fl');

            parse(p, varargin{:});

            StartTime = p.Results.StartTime;
            EndTime = p.Results.EndTime;
            SlowingEndTime = p.Results.SlowingEndTime;
            FreqRef = p.Results.FreqRef;
            Laser = p.Results.Laser;
            kind = p.Results.kind;

            ax = gca;

            ax(1).Parent.Position(3) = 560*3;

            tiledlayout(1,3)
            nexttile
            p=Plot2Dmap(Data1, StartTime = StartTime, EndTime=EndTime,...
                SlowingEndTime = SlowingEndTime, FreqRef = FreqRef, Laser=Laser, kind = kind);
            title("no slowing")
            nexttile
            Plot2Dmap(Data2, StartTime = StartTime, EndTime=EndTime,...
                SlowingEndTime = SlowingEndTime, FreqRef = FreqRef, Laser=Laser, kind = kind);
            title("slowing")
            clim([p.Parent.CLim])
            nexttile
            Plot2DmapSub(Data2, Data1, StartTime = StartTime, EndTime=EndTime,...
                SlowingEndTime = SlowingEndTime, FreqRef = FreqRef, Laser=Laser, kind = kind );
            title("slowing-no slowing")
        end



        function result = PlotFigureSet(self,band)

            clf;
            if nargin==1 || isempty(band)
                band = 'UV';
            end

            if self.datasize.freq ~= 1

                % figure('Name',self.runnum+" abs spec")
                % PlotSpectrum(self, 'abs', band);
                % SaveFigToFile_v2(gcf,"K_results",self.runnum+"_abs_spectrum")
                % SaveFigToFile_v2(gcf,"K_results",self.runnum+"_abs_spectrum")

                % figure('Name',self.runnum+" lif spec")
                % PlotSpectrum(self, 'fl', 2, );
                % SaveFigToFile_v2(gcf,"K_results",self.runnum+"_LIF_spectrum");

                % if self.beamaxis == 'xz'
                %     figure('Name',self.runnum+" lif spec xz")
                % PlotSpectrum(self, gca, 'fl2', band);
                %     SaveFigToFile_v2(gcf,"K_results",self.runnum+"_LIF_spectrum-xz");
                % end

                % elseif self.datasize.freq == 1
                %     figure('Name',self.runnum+" abs tt");
                %     PlotTimetrace(self, 'abs',1);
                %     SaveFigToFile_v2(gcf, "K_results",self.runnum+"_abs_tt");
                %
                %     figure('Name',self.runnum+" lif tt");
                %     PlotTimetrace(self, 'fl',1);
                %     SaveFigToFile_v2(gcf, "K_results",self.runnum+"_lif_tt");
                % end

                % figure('Name',self.runnum+" abs 2D");
                % Plot2Dmap(self, 'abs');
                % SaveFigToFile_v2(gcf, "K_results",self.runnum+"_abs_2D");

                figure('Name',self.runnum+" fl 2D");
                Plot2Dmap(self, 'fl', 2, 1);
                SaveFigToFile_v2(gcf, "K_results",self.runnum+"_lif_2D");

                % if data.beamaxis=='xy'
                %     figure('Name',self.runnum+" fl 2D")
                %     plot2Dtimedet_v3(gca, data, 'fl2');
                %     saveas(gcf, './K_results/'+self.runnum+'_fl_2D_2.png');
                %     saveas(gcf, './K_results/'+self.runnum+'_fl_2D_2.fig');
                % end
                result = gcf;
            end

        end %PlotFigureSet


    end %method
end %classdef

