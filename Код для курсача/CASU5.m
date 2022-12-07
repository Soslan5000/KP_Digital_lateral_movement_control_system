clc
clear

syms s z K11 K12 K14
% Зададим постоянные времени
Tn = 0.11;
Te = 0.11;
% Зададим время дискретизации
T0 = 0.18;

% Зададим матрицы системы уравнений
A = [-0.6, -5.71, -0.04, -2.9, 0.55;...
     1, -0.26, 0.065, -0.04, 0;...
     -0.7, -24, -2.5, -3.3, -19;...
     0, 0, 0, -1/Tn, 0;...
     0, 0, 0, 0, -1/Te];

B = [0, 0;...
     0, 0;...
     0, 0;...
     1/Tn, 0;...
     0, 1/Te];

C = [1, 0, 0, 0, 0;...
     0, 1, 0, 0, 0;...
     0, 0, 1, 0, 0;...
     0, 0, 0, 1, 0;...
     0, 0, 0, 0, 1];

D = 0;

I = eye(5, 5);

% Вычислим матрицу перехода G
G = 0;
n = 20;
for i = 0:n
    G = G + A^i * T0^i / factorial(i);
end

% Вычислим матрицу перехода Dd
Dd = 0;
for i = 0:n
    Dd = Dd + A^(i) * T0^(i+1) / factorial(i+1);
end
Dd = Dd * B;

K11 = 1.6127;
K12 = 2.1617;
K14 = -0.4477;
K23 = 0.05;
K = [K11, K12, 0, K14, 0;
     0, 0, K23, 0, 0];

G_zam = G + Dd*K;

Kx = -0.2799;

W_zam = collect((inv(z*I-G_zam))*Dd);

Wx_zam = W_zam(3, 2);
tf_Wx_zam = convert_to_tf(Wx_zam, true);

figure;
opt = stepDataOptions;
opt.StepAmplitude = Kx;
step(d2c(tf_Wx_zam), opt)
grid on;
title('wx(t) при ступенчатом воздействии Ue');
ylabel('wx(t)', Rotation=0);