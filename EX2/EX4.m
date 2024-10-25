clear all
load("Exercise4.mat")

for i = 1:4
U_0m(i) = WBL(i).U0m;
end
T = WBL(1).T;
nu = WBL(1).nu;

%% 1
% Find reynolds number
omega = 2 * pi / T;     %(P. 301)

a = U_0m ./ omega;      %(P. 238)
Re = U_0m .* a / nu;    %(5.1)


% friction factor laminar
F_w(1) = 2 / sqrt(Re(1)); % (5.59)


% Friction factor transitional
F_w(2) = 0.005;         % (5.61) - Det giver den største friction velocity


% Friction factor turbulent
F_w(3) = 0.035/Re(3)^(0.16);  % (5.60)


k_s(4) = WBL(4).k_s;
F_w(4) = exp(5.5 * (a(4)/k_s(4))^(-0.16) - 6.7);  % (5.69)


% Maximum friction velocity

U_fm = sqrt(F_w/2) .* U_0m;    %(5.57) 


%% 2
Deltay = nu./U_fm; % (9.49)

k_s(1:3) = nu ./ U_fm(1:3) .* 0.1; %(712)



% Expected k+_s
kp_s = k_s .* U_fm ./ nu %(712)