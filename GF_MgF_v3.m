function FitResults = GF_MgF_v3(data)

kB = 1.380649*1e-23; % Boltzmann constant, J/K
c = 299792458; % speed of light, m/s
amu = 1.66053906660*1e-27; % atomic mass unit, kg
hfs = [-129.5 -120.3 0 109.7]; % hyperfine strunture detuning from F=0, [2 1+ 0 1-]

[xData, yData] = prepareCurveData( data.det.UV.mean(:,2,2), data.data.abs.mean);
% [xData, yData] = prepareCurveData( data.det.wm.mean.Q, data.abs.mean);

% Set up fittype and options.
% ft = fittype( 'a1/c1/sqrt(pi)*exp(-((x-b1+129.5)/c1)^2)+a2/c1/sqrt(pi)*exp(-((x-b1+120.3)/c1)^2)+a3/c1/sqrt(pi)*exp(-((x-b1)/c1)^2)+a4/c1/sqrt(pi)*exp(-((x-b1-109.7)/c1)^2)+d1', 'independent', 'x', 'dependent', 'y' );
ft = fittype( 'a1/c1/sqrt(pi)*exp(-((x-b1+129.5)/c1)^2)+a3/c1/sqrt(pi)*exp(-((x-b1)/c1)^2)+a4/c1/sqrt(pi)*exp(-((x-b1-109.7)/c1)^2)+d1', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0 -150 0 0];
opts.StartPoint = [1 1 1 -100 200 0];
% opts.StartPoint = [100 20 100 0 100 200 300 0];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

f0 = 834.294356*1e6; % MHz

u = fitresult.c1/f0*c;

T = u^2*43*amu/(2*kB); % K

Dwidth = 1.7*u/c*f0; % MHz

% FitResults = struct('a1', fitresult.a1, 'a2', fitresult.a2,'a3', fitresult.a3, 'a4', fitresult.a4,...
%     'b1', fitresult.b1,'b2', fitresult.b2, 'b3', fitresult.b3,'b4', fitresult.b4,...
%     'c1', fitresult.c1,'d1', fitresult.d1, 'f0', f0,...
%     'u', u, 'T', T, 'width', Dwidth, 'result', fitresult);

% FitResults = struct('a1', fitresult.a1, 'a2', fitresult.a2,'a3', fitresult.a3, ...
%     'b1', fitresult.b1,'b2', fitresult.b2, 'b3', fitresult.b3,...
%     'c1', fitresult.c1,'d1', fitresult.d1, 'f0', f0,...
%     'u', u, 'T', T, 'width', Dwidth, 'result', fitresult);

FitResults = struct('a1', fitresult.a1, 'a3', fitresult.a3,'a4', fitresult.a4, ...
    'b1', fitresult.b1+hfs(1), 'b3', fitresult.b1+hfs(3),'b4', fitresult.b1+hfs(4),...
    'c1', fitresult.c1,'d1', fitresult.d1, 'f0', f0,...
    'u', u, 'T', T, 'width', Dwidth, 'result', fitresult);

set(groot, 'defaultAxesFontSize',18)
set(groot, 'defaultLineLineWidth',1)
set(groot, 'defaultAxesFontSize',18)
set(groot, 'defaultAxesLineWidth',1.5)
set(groot, 'defaultAxesBox','on')
set(groot, 'DefaultAxesPlotBoxAspectRatio', [1.5 1 1]);

figure1=figure( 'Name', 'G fit' );
h = plot( fitresult, xData, yData ); hold on;

plot(xData, fitresult.a1/fitresult.c1/sqrt(pi)*exp(-((xData-fitresult.b1-hfs(1))/fitresult.c1).^2)); hold on;
% plot(xData, fitresult.a2/fitresult.c1/sqrt(pi)*exp(-((xData-fitresult.b1-hfs(2))/fitresult.c1).^2)); hold on;
plot(xData, fitresult.a3/fitresult.c1/sqrt(pi)*exp(-((xData-fitresult.b1-hfs(3))/fitresult.c1).^2)); hold on;
plot(xData, fitresult.a4/fitresult.c1/sqrt(pi)*exp(-((xData-fitresult.b1-hfs(4))/fitresult.c1).^2)); hold off;
% plot(xData, fitresult.a4/fitresult.c1/sqrt(pi)*exp(-((xData-fitresult.b4)/fitresult.c1).^2)); hold off;
% annotation(figure1,'textbox',...
%     [0.529571428571428 0.490476190476191 0.329357142857143 0.197619047619048],...
%     'String',{"\Delta_{F=2}=  "+num2str(hfs(1)+fitresult.b1)+" MHz, F=0 = "+num2str(fitresult.b1)+" MHz, F=1- = "+num2str(hfs(4)+fitresult.b1)+" MHz"},...
%     'FitBoxToText','on');
legend( 'data', "G fit : T="+num2str(T, 2)+" K",'a1/a3 = '+string(fitresult.a1/fitresult.a3),'$\Delta$ = '+string(fitresult.b1)+' MHz','a4/a3 = '+string(fitresult.a4/fitresult.a3),...
    'Location', 'best', 'Interpreter', 'Latex' );

SaveFigToFile_v2(gcf, "K_results", data.runnum+"_abs_GF")
% saveas(gcf, "./K_results/"+data.savename+"_GF.fig");
% saveas(gcf, "./K_results/"+data.savename+"_GF.png");

end