function FitResults = LF_MgF_v1(data)

kB = 1.380649*1e-23; % Boltzmann constant, J/K
c = 299792458; % speed of light, m/s
amu = 1.66053906660*1e-27; % atomic mass unit, kg
hfs = [-129.5 -120.3 0 109.7]; % hyperfine strunture detuning from F=0, [2 1+ 0 1-]
gamma = 22;
A = gamma^2/4;
C = 0.1;
[xData, yData] = prepareCurveData( data.det.UV.wm.mean, data.fl.mean);

% Set up fittype and options.
% ft = fittype( 'a1/c1/sqrt(pi)*exp(-((x-b1+129.5)/c1)^2)+a2/c1/sqrt(pi)*exp(-((x-b1+120.3)/c1)^2)+a3/c1/sqrt(pi)*exp(-((x-b1)/c1)^2)+a4/c1/sqrt(pi)*exp(-((x-b1-109.7)/c1)^2)+d1', 'independent', 'x', 'dependent', 'y' );
% Lorentz = @(a, b, c, x, f) a*A./((x-b-f).^2+A*(1+c));

% ft = fittype( 'a1*(22^2/4)/((x-b1-249)^2+(22^2/4)*(1+c1))+a2*(22^2/4)/((x-b1-129)^2+(22^2/4)*(1+c1))+a3*(22^2/4)/((x-b1)^2+(22^2/4)*(1+c1))+a4*(22^2/4)/((x-b1+110)^2+(22^2/4)*(1+c1))+a5*(22^2/4)/((x-b1+230)^2+(22^2/4)*(1+c1))+d1',...
%     'independent', 'x', 'dependent', 'y' );
ft = fittype( 'a1*(5*(22^2/4)/((x-b1+129.5)^2+(22^2/4)*(1+c1))+3*(22^2/4)/((x-b1+120.3)^2+(22^2/4)*(1+c1))+1*(22^2/4)/((x-b1)^2+(22^2/4)*(1+c1))+3*(22^2/4)/((x-b1-109.7)^2+(22^2/4)*(1+c1)))+d1',...
    'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 -50  -100 0 ];
opts.StartPoint = [1e4  0 100 0];
% opts.StartPoint = [100 20 100 0 100 200 300 0];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

f0 = 834.294356*1e6; % MHz

% u = fitresult.c1/f0*c;

% T = u^2*43*amu/(2*kB); % K

% Dwidth = 1.7*u/c*f0; % MHz

% FitResults = struct('a1', fitresult.a1, 'a2', fitresult.a2,'a3', fitresult.a3, 'a4', fitresult.a4,...
%     'b1', fitresult.b1,'b2', fitresult.b2, 'b3', fitresult.b3,'b4', fitresult.b4,...
%     'c1', fitresult.c1,'d1', fitresult.d1, 'f0', f0,...
%     'u', u, 'T', T, 'width', Dwidth, 'result', fitresult);

% FitResults = struct('a1', fitresult.a1, 'a2', fitresult.a2,'a3', fitresult.a3, ...
%     'b1', fitresult.b1,'b2', fitresult.b2, 'b3', fitresult.b3,...
%     'c1', fitresult.c1,'d1', fitresult.d1, 'f0', f0,...
%     'u', u, 'T', T, 'width', Dwidth, 'result', fitresult);

% FitResults = struct('a1', fitresult.a1, 'a3', fitresult.a3,'a4', fitresult.a4, ...
%     'b1', fitresult.b1+hfs(1), 'b3', fitresult.b1+hfs(3),'b4', fitresult.b1+hfs(4),...
%     'c1', fitresult.c1,'d1', fitresult.d1, 'f0', f0,...
%     'u', u, 'T', T, 'width', Dwidth, 'result', fitresult);

figure1=figure( 'Name', 'L fit' );
h = plot( fitresult, xData, yData ); hold on;

FitResults = fitresult;

% plot(xData, fitresult.a1/fitresult.c1/sqrt(pi)*exp(-((xData-fitresult.b1+hfs(1))/fitresult.c1).^2)); hold on;
% plot(xData, fitresult.a2/fitresult.c1/sqrt(pi)*exp(-((xData-fitresult.b1+hfs(2))/fitresult.c1).^2)); hold on;
% plot(xData, fitresult.a3/fitresult.c1/sqrt(pi)*exp(-((xData-fitresult.b1+hfs(3))/fitresult.c1).^2)); hold on;
% plot(xData, fitresult.a4/fitresult.c1/sqrt(pi)*exp(-((xData-fitresult.b1+hfs(4))/fitresult.c1).^2)); hold off;
% plot(xData, fitresult.a4/fitresult.c1/sqrt(pi)*exp(-((xData-fitresult.b4)/fitresult.c1).^2)); hold off;

plot(xData, fitresult.a1*5*(22^2/4)./((xData-fitresult.b1-hfs(1)).^2+(22^2/4)*(1+fitresult.c1))); hold on;
plot(xData, fitresult.a1*3*(22^2/4)./((xData-fitresult.b1-hfs(2)).^2+(22^2/4)*(1+fitresult.c1))); hold on;
plot(xData, fitresult.a1*(22^2/4)./((xData-fitresult.b1-hfs(3)).^2+(22^2/4)*(1+fitresult.c1))); hold on;
plot(xData, fitresult.a1*3*(22^2/4)./((xData-fitresult.b1-hfs(4)).^2+(22^2/4)*(1+fitresult.c1))); hold off;

% annotation(figure1,'textbox',...
%     [0.529571428571428 0.490476190476191 0.329357142857143 0.197619047619048],...
%     'String',{"\Delta_{F=2}=  "+num2str(hfs(1)+fitresult.b1)+" MHz, F=0 = "+num2str(fitresult.b1)+" MHz, F=1- = "+num2str(hfs(4)+fitresult.b1)+" MHz"},...
%     'FitBoxToText','on');
legend( 'data', "fit result",'F=2','F=1+','F=0','F=1-',...
    'Location', 'best', 'Interpreter', 'Latex' );
% 
saveas(gcf, "./K_results/"+data.savename+"_GF.fig");
saveas(gcf, "./K_results/"+data.savename+"_GF.png");

end