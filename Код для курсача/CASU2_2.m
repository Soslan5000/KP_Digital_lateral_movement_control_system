clc
clear

syms s z L1 L2 L4 lambd
% Зададим постоянную времени
Tn = 0.11;
% Зададим время дискретизации
T0 = 0.18;

% Зададим матрицы системы уравнений
A = [-0.6, -5.71, -2.9;...
     1, -0.26, -0.04;...
     0, 0, -1/Tn];

B = [0;...
     0;...
     1/Tn];

C = [0, 1, 0];

D = 0;

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

I = eye(3, 3);

K = [1.6127, 2.1617, -0.4477];

L = [L1;...
     L2;...
     L4];

% Вычислим характеристический многочлен замкнутого фильтра Калмана
Dz = vpa(collect(det(lambd*I - (G + Dd*K - L*C)), lambd));
coef = coeffs(Dz, lambd);
a0 = coef(1);
a1 = coef(2);
a2 = coef(3);
a3 = coef(4);

% Получим желаемый характеристический многочлен
tg = 0.3;
tn = 6.3;
t = tg / tn;
Dg = vpa(collect((lambd + 1/t)^3));
TF_DG = convert_to_tf(1/Dg, false);
TFd_DG = c2d(TF_DG, T0);
b0 = TFd_DG.denominator{1}(4);
b1 = TFd_DG.denominator{1}(3);
b2 = TFd_DG.denominator{1}(2);
b3 = TFd_DG.denominator{1}(1);

eq0 = a0 == b0;
eq1 = a1 == b1;
eq2 = a2 == b2;
root = struct2array(solve(eq0, eq1, eq2, L1, L2, L4));

L(1) = root(1);
L(2) = root(2);
L(3) = root(3);
L = double(L);

Q = [G + Dd*K, zeros(3, 3);...
     L*C, G+Dd*K-L*C];
R = [1, 0, 0, 0, 0, 0;
     0, 1, 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0;
     0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 1, 0;
     0, 0, 0, 0, 0, 1];
P = [Dd;
     Dd];

W = inv(z*eye(6) - Q) * P;

W_1 = convert_to_tf(vpa(collect(W(1))), true)
W_2 = convert_to_tf(vpa(collect(W(2))), true)
W_3 = convert_to_tf(vpa(collect(W(3))), true)
W_4 = convert_to_tf(vpa(collect(W(4))), true)
W_5 = convert_to_tf(vpa(collect(W(5))), true)
W_6 = convert_to_tf(vpa(collect(W(6))), true)

figure;
step(W_1, W_4)
figure;
step(W_2, W_5)
figure;
step(W_3, W_6)
