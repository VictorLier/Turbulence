clear all; close all; clc;
load("MatRANS/TurbulenceBook/Exercise4/Case2/out_MatRANS.mat");
load("Exercise4.mat")

%% Bed shear stress
u = MatRANS.u;
y = MatRANS.y;
U_0m = MatRANS.U0m;
n_t = MatRANS.n_t;
t = MatRANS.t;
T = MatRANS.T;
nu = MatRANS.nu;

omega = 2 * pi / T;

a = U_0m/omega

% Stokes length
sigma_1 = (2 * nu / omega)^(1/2); % (5.13)

period_data_length = (n_t-1)/5;
omegat(1) = period_data_length*4; % 0
omegat(2) = omegat(1) + period_data_length/360*45; % 45
omegat(3) = omegat(2) + period_data_length/360*45; % 90
omegat(4) = omegat(3) + period_data_length/360*45; % 135


% bed shear stress time series
tau = MatRANS.tau0(omegat(1)+2:end);
rho = MatRANS.rho;

mu = nu*rho;

angle = linspace(0,2*pi, period_data_length);

tau_t = sqrt(2) * U_0m * mu/sigma_1 * sin(angle + pi/4); %(5.20)

bss_t = tau_t / (rho * U_0m^2);

%angle = linspace(0,360,period_data_length);


bss = tau / (rho * U_0m^2);

tau_d = WBL(2).tau0;
bss_d = tau_d / (rho * U_0m^2);
angle_d = WBL(2).omegat_tau0

figure(2)
hold on
plot(angle*180/pi, bss, Color='r')
plot(angle*180/pi,bss_t, Color='b')
plot(angle_d,bss_d, Color='g')
hold off

x = angle*180/pi;

bss_data = [x; bss']';

writematrix(bss_data, 'Case2Data/tau.txt', Delimiter=' ');
writematrix([angle*180/pi; bss_t]', 'Case2Data/tau_t.txt', Delimiter=' ');
writematrix([angle_d; bss_d]', 'Case2Data/tau_d.txt', Delimiter=' ')


%% Friction coefficient
Fw_model = (2 .* tau')./rho./(U_0m^2 .* sin(angle + pi/4)); %(5.22)

Fw_data = (2 .* tau_d)./rho./(U_0m^2 .* sin((angle_d*pi/180) + pi/4)); %(5.22)

Re = a * U_0m/nu

f_w = 2/sqrt(Re) %(5.59)

figure(3)
hold on
plot(angle*180/pi,Fw_model, Color='r')
plot(angle_d,Fw_data,Color='b')
yline(f_w, Color='g')
hold off

writematrix([angle*180/pi;Fw_model]', 'Case2Data/Fw_model.txt', Delimiter=' ');
writematrix([angle_d;Fw_data]', 'Case2Data/Fw_data.txt', Delimiter=' ');
writematrix(f_w, 'Case2Data/Fw_laminar.txt', Delimiter=' ')


