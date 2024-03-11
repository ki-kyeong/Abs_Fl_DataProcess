function results = EOMfit_v1(x,y,fm,n,sps)

% 2023.12.28. 생성
% parameter 설명 추가
% x,y    : n*1 double each
% x      : time or frequency
% y      : PD voltage
% fm     : modulation frequency
% sps    : {'amp','Gamma','x0','beta','offset'}
%        : Fitting value start point에요.

%    amp    : maximum peak size
%    Gamma  : linewidth of lorentzian signal (typically, 7 MHz)
%    x0     : x-axis(time or frequency) offset
%    beta   : modulation index
%    offset : PD voltage offset
% ========== % ========== % ========== % ========== % ========== % ========== %
% Fitting이 잘 안되면 이렇게 해보세요.(주로 r^2값이 0.7 이하일 때 쓰면 유용함)
% 1. lower bound / upper bound 살짝 조정해보기
% 2. sps를 직접 설정해서 집어넣기
%    -> 이 때는 함수에서 sps를 없애고 아래에서 직접 설정해야 합니다.

my_EOM_fit_function = "((amp*Gamma^2/((x-x0)^2+Gamma^2))*(besselj(0,beta)^2)) + offset";

for i = 1:n
    fitf_i = "((amp*Gamma^2)/((x-x0-"+num2str(i*fm)+")^2+Gamma^2)+(amp*Gamma^2)/((x-x0+"+num2str(i*fm)+")^2+(Gamma^2)))*(besselj("+num2str(i)+",beta)^2)";
    my_EOM_fit_function = my_EOM_fit_function + " + " + fitf_i;
end

results.model=my_EOM_fit_function;

myfittype = fittype(my_EOM_fit_function,...
        'dependent','y','independent',{'x'},...
        'coefficients',{'amp','Gamma','x0','beta','offset'});


[f,gof] = fit(x,y,myfittype,...
    'startpoint',sps,...[amp, Gamma, x0, beta, DCoffset],...
    'lower',[0 0 -Inf 0 -Inf],'upper',[Inf 10 Inf Inf Inf]);

eom_fit = figure(1)
plot(f,x,y)

title("EOM fitting results, \beta = " + num2str(f.beta) + " rad")
xlabel("Frequency (MHz)")
ylabel("Cavity transmission signal (V)")
grid on

results.amp = f.amp;
results.Gamma = f.Gamma;
results.x0 = f.x0;
results.beta = f.beta;
% results.fit = eom_fit
results.DC_offset = f.offset;

results.r2 = gof.rsquare;
% 
% plotlegend = legend("Data points","fitting curve");
% % plotlegend.Title.String = "Fit data type : " + num2str(type);
% 
% lp = plotlegend.Position; % LegendPosition -> lp
% 
% annotation('textbox',[lp(1) 0.85*lp(2) lp(3) lp(4)],...
%     'String',["\beta = "+num2str(results.beta,'%.3f')+" rad" "r^{2} = "+num2str(results.r2,'%.3f')],...
%     'VerticalAlignment','middle',...
%     'BackgroundColor','w',...
%     'FitBoxToText','on')