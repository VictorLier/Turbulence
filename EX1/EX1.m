clear all
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

semilogx(ylim, ubarlim)


xyzb2 = polyfit(log(ylim), ubarlim, 1);
slope = xyzb2(1);
offset = xyzb2(2);

newA = 2.5;

newU_f = slope / newA

%% AI 5
dim1 = ubar/newU_f
dim2 = y/h

figure()
plot(yplus, dim1)

figure()
plot(dim2, dim1)





