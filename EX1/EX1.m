clear all
clc
close all
load("Exercise1.mat")

%% AI 1
n  = length(Channel);

for j=1:n;

N = length(Channel(j).tt);
a1 = 0;
b1 = 0;
for i=1:N;
    a1 = a1 + Channel(j).u(i) * Channel(j).tt(i);
    b1 = b1 + Channel(j).tt(i);
end

ubar(j) = a1/b1;

y(j) = Channel(j).y;

end



ubar = [0, ubar, 0.3];

h = Channel(1).h;
y = [0, y, h];

% plot(y,ubar)
% title("Mean velocity as a function of y")
% ylabel("Ubar")
% xlabel("y")

%% AI 2
V = 1/h * trapz(y,ubar)

%% AI 3

b = 0.3;
nu = Channel(1).nu;

A = h*b;
P = 2*h + b;

r_h = A/P;

Re = r_h*V/nu;

f = 0.0557/Re^0.25

U_f = sqrt(f/2)*V

%% AI 4

Re_tau = h * U_f/nu;

yplus = y*U_f/nu;
ubarlim = 0;
ylim = 0;
for i=1:length(ubar);
    if yplus(i) <= 0.1*Re_tau && yplus(i) >= 30;
        ubarlim = [ubarlim, ubar(i)];
        ylim = [ylim, y(i)];
    end
end
ubarlim = ubarlim(2:end);
ylim = ylim(2:end);

semilogy(ubarlim,ylim)


xyzb2 = polyfit(log(ylim), ubarlim, 1);
slope = xyzb2(1);
offset = xyzb2(2);

newA = 2.5;

newU_f = slope / newA;

%% AI 5
dim1 = ubar/newU_f;
dim2 = y/h;

figure()
semilogy(dim1,yplus)
yline(5)
yline(30)
yline(0.1*Re_tau)
xlabel("ubar/U_f")
ylabel("y+")
title("ubar/newU_f")

figure()
semilogy(dim1, dim2)
xlabel("ubar/U_f")
ylabel("y/h")
title("y/h")


%% AI 6

kappa = 0.4;
A_d = 25;


vDriest = 2.*newU_f.*cumtrapz(1./(1+(1+4.*kappa.^2*yplus.^2.*(1-exp(-yplus./A_d)).^2).^(1/2)));

figure
semilogy(vDriest/newU_f,yplus)
hold on
semilogy(dim1,yplus)
yline(5)
yline(30)
yline(0.1*Re_tau)
legend("vDriest","Ubar/U_f")


%% AI 7


for jj=1:n

  udiff = diff(Channel(jj).u);
  turb_u(jj) = sqrt(mean(udiff.^2));

  vdiff = diff(Channel(jj).v);
  turb_v(jj) = sqrt(mean(vdiff.^2));

  turb_uv(jj) = sqrt(-mean(udiff.*vdiff));
end

figure()
plot(yplus(2:end-1),turb_u/newU_f)
hold on
plot(yplus(2:end-1),turb_v/newU_f)
plot(yplus(2:end-1),turb_uv/newU_f)
xlim([0 100])
legend("turb_u","turb_v","turb_uv")

%% AI 8

figure()
plot(y(2:end-1)/h,turb_u/newU_f)
hold on
plot(y(2:end-1)/h,turb_v/newU_f)
plot(y(2:end-1)/h,turb_uv/newU_f)
xlim([0 1])
legend("turb_u","turb_v","turb_uv")


%% AI 9

turb_w = 1.8*turb_v.^2;
k = 1/2 * (turb_u.^2+turb_v.^2+turb_w.^2);

figure()
plot(y(2:end-1)/h, k/newU_f^2)