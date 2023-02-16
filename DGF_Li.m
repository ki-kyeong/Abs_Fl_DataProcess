function FitResults = DGF_Li(x, y)

kB = 1.380649*1e-23; % Boltzmann constant, J/K
c = 299792458; % speed of light, m/s
amu = 1.66053906660*1e-27; % atomic mass unit, kg

[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'a1/c1/sqrt(pi)*exp(-((x-b1)/c1)^2)+a2/c1/sqrt(pi)*exp(-((x-b2)/c1)^2)+d1', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% opts.Lower = [0 0 -Inf 0 0 0];
opts.StartPoint = [200 100 0 770 200 0];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

a1 = fitresult.a1;
b1 = fitresult.b1;
c1 = fitresult.c1;
d1 = fitresult.d1;
a2 = fitresult.a2;
b2 = fitresult.b2;

f0 = 446.80987*1e6+b1; % MHz

u = c1/f0*c;

T = u^2*7*amu/(2*kB); % K

Dwidth = 1.7*u/c*f0; % MHz

FitResults = struct('a1', fitresult.a1, 'a2', fitresult.a2, 'b1', fitresult.b1,...
    'b2', fitresult.b2, 'c1', fitresult.c1,'d1', fitresult.d1, 'f0', f0,...
    'u', u, 'T', T, 'width', Dwidth, 'result', fitresult);

figure1=figure( 'Name', 'd-G fit' );
h = plot( fitresult, xData, yData );
legend( h, 'yData vs xData', 'd-G fit : T='+string(T), 'Location', 'best', 'Interpreter', 'none' );

% annotation(figure1,'textbox',...
%     [0.529571428571428 0.490476190476191 0.329357142857143 0.197619047619048],...
%     'String',{'T = '+string(T)+' K','FWHM = '+string(Dwidth)+' MHz','det1 = '+string(b1)+' MHz','det2 = '+string(b2)+' MHz'},...
%     'FitBoxToText','on');

end