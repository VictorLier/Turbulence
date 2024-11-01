clear all
close all
clc

%% AI1
load("Exercise3.mat")
u = Jet(12).u;
t = Jet(12).t;
ubar = mean(u);
uprime = u-ubar;
%plot(t,uprime)
%xlim([6,6.25])
writematrix([t uprime], 'AI1.txt', Delimiter=' ')

%% AI2

ubar;
sigma2u = var(u);
skew = skewness(u);
kurt = kurtosis(u);
sigmau = sqrt(sigma2u);
turb_int = sigmau/ubar;


%% AI3

p = pdf('Normal',uprime,mean(uprime),sigmau);
p_calc = (1/(sigmau.*sqrt(2.*pi))).*exp(-(uprime.^2)./(2.*sigma2u));

%figure()
%scatter(uprime,p, Color='r')
%hold on
%scatter(uprime,p_calc, Color='b')
%legend('pdf','calculated')
%hold off

%writematrix([uprime p], 'AI3_pdf.txt', Delimiter=' ')
%writematrix([uprime p_calc], 'AI3_pdfcalc.txt', Delimiter=' ')

%% AI4




%% AI5


%% AI6


%% AI7


%% AI8


%% AI9



%% AI10

%% AI11