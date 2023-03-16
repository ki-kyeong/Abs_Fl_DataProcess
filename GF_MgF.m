function FitResults = GF_MgF(x, y)

kB = 1.380649*1e-23; % Boltzmann constant, J/K
c = 299792458; % speed of light, m/s
amu = 1.66053906660*1e-27; % atomic mass unit, kg

[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
% ft = fittype( 'a1/c1/sqrt(pi)*exp(-((x-b1)/c1)^2)+a2/c1/sqrt(pi)*exp(-((x-b2)/c1)^2)+a3/c1/sqrt(pi)*exp(-((x-b3)/c1)^2)+a4/c1/sqrt(pi)*exp(-((x-b4)/c1)^2)+d1', 'independent', 'x', 'dependent', 'y' );
ft = fittype( 'a1/c1/sqrt(pi)*exp(-((x-b1)/c1)^2)+a2/c1/sqrt(pi)*exp(-((x-b2)/c1)^2)+a3/c1/sqrt(pi)*exp(-((x-b3)/c1)^2)+d1', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% opts.Lower = [0 0 -Inf 0 0 0];
% opts.StartPoint = [100 100 20 100 0 0 100 200 200 0];
opts.StartPoint = [100 20 100 0 100 200 300 0];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

f0 = 834.294356*1e6+fitresult.b1; % MHz

u = fitresult.c1/f0*c;

T = u^2*7*amu/(2*kB); % K

Dwidth = 1.7*u/c*f0; % MHz

% FitResults = struct('a1', fitresult.a1, 'a2', fitresult.a2,'a3', fitresult.a3, 'a4', fitresult.a4,...
%     'b1', fitresult.b1,'b2', fitresult.b2, 'b3', fitresult.b3,'b4', fitresult.b4,...
%     'c1', fitresult.c1,'d1', fitresult.d1, 'f0', f0,...
%     'u', u, 'T', T, 'width', Dwidth, 'result', fitresult);

FitResults = struct('a1', fitresult.a1, 'a2', fitresult.a2,'a3', fitresult.a3, ...
    'b1', fitresult.b1,'b2', fitresult.b2, 'b3', fitresult.b3,...
    'c1', fitresult.c1,'d1', fitresult.d1, 'f0', f0,...
    'u', u, 'T', T, 'width', Dwidth, 'result', fitresult);

figure1=figure( 'Name', 'G fit' );
h = plot( fitresult, xData, yData ); hold on;
legend( h, 'yData vs xData', 'G fit : T='+string(T), 'Location', 'best', 'Interpreter', 'none' );
plot(xData, fitresult.a1/fitresult.c1/sqrt(pi)*exp(-((xData-fitresult.b1)/fitresult.c1).^2)); hold on;
plot(xData, fitresult.a2/fitresult.c1/sqrt(pi)*exp(-((xData-fitresult.b2)/fitresult.c1).^2)); hold on;
plot(xData, fitresult.a3/fitresult.c1/sqrt(pi)*exp(-((xData-fitresult.b3)/fitresult.c1).^2)); hold off;
% plot(xData, fitresult.a4/fitresult.c1/sqrt(pi)*exp(-((xData-fitresult.b4)/fitresult.c1).^2)); hold off;
% annotation(figure1,'textbox',...
%     [0.529571428571428 0.490476190476191 0.329357142857143 0.197619047619048],...
%     'String',{'T = '+string(T)+' K','FWHM = '+string(Dwidth)+' MHz','det1 = '+string(b1)+' MHz','det2 = '+string(b2)+' MHz'},...
%     'FitBoxToText','on');

end