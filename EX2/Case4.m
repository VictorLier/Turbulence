clear all; close all; clc;
load("MatRANS/TurbulenceBook/Exercise4/Case4/out_MatRANS.mat");
load("Exercise4.mat")

%% Data velocity
U_0m = WBL(4).U0m;

y_data = WBL(4).y_u;

omegat_data(1) = WBL(4).omegat(1);
omegat_data(2) = WBL(4).omegat(4);
omegat_data(3) = WBL(4).omegat(7);
omegat_data(4) = WBL(4).omegat(10);

u_data_0 = WBL(4).u(:,1);
u_data_45 = WBL(4).u(:,4);
u_data_90 = WBL(4).u(:,7);
u_data_135 = WBL(4).u(:,10);


omega = 2 * pi / WBL(4).T;

a = U_0m/omega;

writematrix([u_data_0/U_0m, y_data/a], 'Case4Data/d_vel_0.txt',Delimiter=' ');
writematrix([u_data_45/U_0m, y_data/a], 'Case4Data/d_vel_45.txt',Delimiter=' ');
writematrix([u_data_90/U_0m, y_data/a], 'Case4Data/d_vel_90.txt',Delimiter=' ');
writematrix([u_data_135/U_0m, y_data/a], 'Case4Data/d_vel_135.txt',Delimiter=' ');

%% Model Velocity
n_t = MatRANS.n_t - 1;
one_phase = n_t/5;
fortyfive = one_phase/8;
omegat_0 = n_t - one_phase+1;
omegat_45 = omegat_0 + fortyfive;
omegat_90 = omegat_45 + fortyfive;
omegat_135 = omegat_90 + fortyfive;

y_model = MatRANS.y';
u_model_0 = MatRANS.u(omegat_0,:)';
u_model_45 = MatRANS.u(omegat_45,:)';
u_model_90 = MatRANS.u(omegat_90,:)';
u_model_135 = MatRANS.u(omegat_135,:)';

writematrix([u_model_0/U_0m, y_model/a], 'Case4Data/m_vel_0.txt',Delimiter=' ');
writematrix([u_model_45/U_0m, y_model/a], 'Case4Data/m_vel_45.txt',Delimiter=' ');
writematrix([u_model_90/U_0m, y_model/a], 'Case4Data/m_vel_90.txt',Delimiter=' ');
writematrix([u_model_135/U_0m, y_model/a], 'Case4Data/m_vel_135.txt',Delimiter=' ');

%figure(1)
%hold on
%plot(u_data_45/U_0m, y_data/a)
%plot(u_model_45/U_0m, y_model/a)
%hold off


%% Data turbulent kinetic energy
yumb_data = WBL(4).y_uuvv;
umb_data_0 = WBL(4).uu(:,1);
umb_data_45 = WBL(4).uu(:,4);
umb_data_90 = WBL(4).uu(:,7);
umb_data_135 = WBL(4).uu(:,10);

vmb_data_0 = WBL(4).vv(:,1);
vmb_data_45 = WBL(4).vv(:,4);
vmb_data_90 = WBL(4).vv(:,7);
vmb_data_135 = WBL(4).vv(:,10);

k_data_0 = 0.65 * (umb_data_0 + vmb_data_0);
k_data_45 = 0.65 * (umb_data_45 + vmb_data_45);
k_data_90 = 0.65 * (umb_data_90 + vmb_data_90);
k_data_135 = 0.65 * (umb_data_135 + vmb_data_135);

%plot(k_0/U_0m^2, yumb_data/a)

writematrix([k_data_0/U_0m^2, yumb_data/a], 'Case4Data/d_k_0.txt',Delimiter=' ');
writematrix([k_data_45/U_0m^2, yumb_data/a], 'Case4Data/d_k_45.txt',Delimiter=' ');
writematrix([k_data_90/U_0m^2, yumb_data/a], 'Case4Data/d_k_90.txt',Delimiter=' ');
writematrix([k_data_135/U_0m^2, yumb_data/a], 'Case4Data/d_k_135.txt',Delimiter=' ');


%% Model turbulent kinetic energy
yumb_model = MatRANS.y';

k_model_0 = MatRANS.k(omegat_0,:)';
k_model_45 = MatRANS.k(omegat_45,:)';
k_model_90 = MatRANS.k(omegat_90,:)';
k_model_135 = MatRANS.k(omegat_135,:)';

%figure(2)
%hold on
%plot(k_data_45/U_0m, yumb_data/a)
%plot(k_model_45/U_0m, yumb_model/a)
%hold off

writematrix([k_model_0/U_0m^2, yumb_model/a], 'Case4Data/m_k_0.txt',Delimiter=' ');
writematrix([k_model_45/U_0m^2, yumb_model/a], 'Case4Data/m_k_45.txt',Delimiter=' ');
writematrix([k_model_90/U_0m^2, yumb_model/a], 'Case4Data/m_k_90.txt',Delimiter=' ');
writematrix([k_model_135/U_0m^2, yumb_model/a], 'Case4Data/m_k_135.txt',Delimiter=' ');


%% Data Reynolds stress
yuv_data = WBL(4).y_uv;
uv_data_0 = WBL(4).uv(:,1);
uv_data_45 = WBL(4).uv(:,4);
uv_data_90 = WBL(4).uv(:,7);
uv_data_135 = WBL(4).uv(:,10);

%plot(uv_data_0,yuv_data)

writematrix([-uv_data_0/U_0m^2, yuv_data/a], 'Case4Data/d_r_0.txt',Delimiter=' ');
writematrix([-uv_data_45/U_0m^2, yuv_data/a], 'Case4Data/d_r_45.txt',Delimiter=' ');
writematrix([-uv_data_90/U_0m^2, yuv_data/a], 'Case4Data/d_r_90.txt',Delimiter=' ');
writematrix([-uv_data_135/U_0m^2, yuv_data/a], 'Case4Data/d_r_135.txt',Delimiter=' ');


%% Model Reyonlds stress
nut_0 = MatRANS.nu_t(omegat_0,:)';
nut_45 = MatRANS.nu_t(omegat_45,:)';
nut_90 = MatRANS.nu_t(omegat_90,:)';
nut_135 = MatRANS.nu_t(omegat_135,:)';



uv_model_0 =  nut_0 .* gradient(u_model_0, y_model);
uv_model_45 = nut_45 .* gradient(u_model_45, y_model);
uv_model_90 = nut_90 .* gradient(u_model_90, y_model);
uv_model_135 = nut_135 .* gradient(u_model_135, y_model);

writematrix([uv_model_0/U_0m^2, y_model/a], 'Case4Data/m_r_0.txt',Delimiter=' ');
writematrix([uv_model_45/U_0m^2, y_model/a], 'Case4Data/m_r_45.txt',Delimiter=' ');
writematrix([uv_model_90/U_0m^2, y_model/a], 'Case4Data/m_r_90.txt',Delimiter=' ');
writematrix([uv_model_135/U_0m^2, y_model/a], 'Case4Data/m_r_135.txt',Delimiter=' ');


%figure(3)
%hold on
%plot(-uv_data_135/(U_0m^2), yuv_data/a)
%plot(uv_model_135/(U_0m^2), y_model/a)
%hold off


%% Data bed shear stress time series
rho = MatRANS.rho

omega_data = WBL(4).omegat_tau0;
tau_data = WBL(4).tau0;

writematrix([omega_data, tau_data/(rho * U_0m^2)], 'Case4Data/data_shear.txt',Delimiter=' ');




%plot(omega_data, tau_data)


%% Model bed shear stress time series
omega_model = linspace(0, 360, one_phase);
tau_model = MatRANS.tau0(omegat_0+1:end);

writematrix([omega_model', tau_model/(rho * U_0m^2)], 'Case4Data/model_shear.txt',Delimiter=' ');

%
%figure(4)
%hold on
%plot(omega_data, tau_data/(rho * U_0m^2))
%plot(omega_model, tau_model/(rho * U_0m^2))
%hold off

