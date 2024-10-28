clear all; close all; clc;
load("MatRANS/TurbulenceBook/Exercise4/Case1/out_MatRANS.mat");

u = MatRANS.u;
y = MatRANS.y;
U_0m = MatRANS.U0m;
a = 0.355806790776241;
n_t = MatRANS.n_t;
t = MatRANS.t
T = MatRANS.T
nu = MatRANS.nu

omega = 2 * pi / T

% Stokes length
sigma_1 = (2 * nu / omega)^(1/2) % (5.13)

% Theoretical velocity (5.12)
u_theory1 = U_0m .* sin(0*pi/180) - U_0m .* exp(-y./sigma_1) .* sin((0*pi/180) - y./sigma_1);
u_theory2 = U_0m .* sin(45*pi/180) - U_0m .* exp(-y./sigma_1) .* sin((45*pi/180) - y./sigma_1);
u_theory3 = U_0m .* sin(90*pi/180) - U_0m .* exp(-y./sigma_1) .* sin((90*pi/180) - y./sigma_1);
u_theory4 = U_0m .* sin(135*pi/180) - U_0m .* exp(-y./sigma_1) .* sin((135*pi/180) - y./sigma_1);

period_data_length = (n_t-1)/5;

omegat(1) = period_data_length*4+1;
omegat(2) = omegat(1) + period_data_length/360*45+1;
omegat(3) = omegat(2) + period_data_length/360*45+1;
omegat(4) = omegat(3) + period_data_length/360*45+1;


A0 = u(omegat(1),:)/U_0m;
A45 = u(omegat(2),:)/U_0m;
A90 = u(omegat(3),:)/U_0m;
A135 = u(omegat(4),:)/U_0m;
A0_t = u_theory1/U_0m;
A45_t = u_theory2/U_0m;
A90_t = u_theory3/U_0m;
A135_t = u_theory4/U_0m;
B = y/a;
data0 = [A0; B]';
data45 = [A45; B]';
data90 = [A90; B]';
data135 = [A135; B]';
data0_t = [A0_t; B]';
data45_t = [A45_t; B]';
data90_t = [A90_t; B]';
data135_t = [A135_t; B]';
writematrix(data0, 'Case1Data/0.txt', Delimiter=' ');
writematrix(data45, 'Case1Data/45.txt', Delimiter=' ');
writematrix(data90, 'Case1Data/90.txt', Delimiter=' ');
writematrix(data135, 'Case1Data/135.txt', Delimiter=' ');
writematrix(data0_t, 'Case1Data/0_t.txt', Delimiter=' ');
writematrix(data45_t, 'Case1Data/45_t.txt', Delimiter=' ');
writematrix(data90_t, 'Case1Data/90_t.txt', Delimiter=' ');
writematrix(data135_t, 'Case1Data/135_t.txt', Delimiter=' ');




% Create a 1x4 grid of subplots
figure(1)
subplot(1, 4, 1); % First subplot
plot(u(omegat(1),:)/U_0m,y/a)
plot(u_theory1/U_0m,y/a)
title('Plot 1');
xlim([0,1.5])
ylim([0,0.02])

subplot(1, 4, 2); % Second subplot
plot(u(omegat(2),:)/U_0m,y/a)
plot(u_theory2/U_0m,y/a)
title('Plot 2');
xlim([0,1.5])
ylim([0,0.02])

subplot(1, 4, 3); % Third subplot
plot(u(omegat(3),:)/U_0m,y/a)
plot(u_theory3/U_0m,y/a)
title('Plot 3');
xlim([0,1.5])
ylim([0,0.02])

subplot(1, 4, 4); % Fourth subplot
plot(u(omegat(4),:)/U_0m,y/a)
plot(u_theory4/U_0m,y/a)
title('Plot 4')
xlim([0,1.5])
ylim([0,0.02])



% bed shear stress time series
tau = MatRANS.tau0(omegat(1)+1:end)
rho = MatRANS.rho

mu = nu*rho

angle = linspace(0,2*pi, period_data_length)

tau_t = sqrt(2) * U_0m * mu/sigma_1 * sin(angle + pi/4) %(5.20)

bss_t = tau_t / (rho * U_0m^2)

%angle = linspace(0,360,period_data_length);




bss = tau / (rho * U_0m^2)

figure(2)
hold on
plot(angle*180/pi, bss, Color='r')
plot(angle*180/pi,bss_t, Color='b')
hold off

x = angle*180/pi

bss_data = [x; bss']'

writematrix(bss_data, 'Case1Data/tau.txt', Delimiter=' ');
writematrix([angle*180/pi; bss_t]', 'Case1Data/tau_t.txt', Delimiter=' ');


