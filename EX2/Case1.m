clear all; close all; clc;
load("MatRANS/TurbulenceBook/Exercise4/Case1/out_MatRANS.mat");

u = MatRANS.u;
y = MatRANS.y;
U_0m = MatRANS.U0m;
a = 0.355806790776241;
n_t = MatRANS.n_t;
t = MatRANS.t

period_data_length = (n_t-1)/5

omegat(1) = period_data_length*4
omegat(2) = omegat(1) + period_data_length/360*45
omegat(3) = omegat(2) + period_data_length/360*45
omegat(4) = omegat(3) + period_data_length/360*45







% Create a 1x4 grid of subplots
subplot(1, 4, 1); % First subplot
plot(u(omegat(1),:)/U_0m,y/a)
title('Plot 1');
xlim([0,1.5])
ylim([0,0.02])

subplot(1, 4, 2); % Second subplot
plot(u(omegat(2),:)/U_0m,y/a)
title('Plot 2');
xlim([0,1.5])
ylim([0,0.02])

subplot(1, 4, 3); % Third subplot
plot(u(omegat(3),:)/U_0m,y/a)
title('Plot 3');
xlim([0,1.5])
ylim([0,0.02])

subplot(1, 4, 4); % Fourth subplot
plot(u(omegat(4),:)/U_0m,y/a)
title('Plot 4')
xlim([0,1.5])
ylim([0,0.02])