clc
clear

syms s z K11 K12 K14
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

C = [1, 0, 0;...
     0, 1, 0;...
     0, 0, 1];

D = 0;

I = eye(3, 3);

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

% Вычислим матрицу передаточных функций разомкнутой цифровой системы
Wz = collect((inv(z*I-G))*Dd);

% Выделим передаточные функции по всем каналам
W_omy_Un = collect(Wz(1));
W_b_Un = collect(Wz(2));
W_sigmn_Un = collect(Wz(3));
% Преобразуем ПФ к типу tf
tf_W_omy_Un = convert_to_tf(Wz(1), true);
tf_W_b_Un = convert_to_tf(Wz(2), true);
tf_W_sigmn_Un = convert_to_tf(Wz(3), true);

% Вычислим передаточные функции замкнутого цифрового контура и приведём их
% к нормальному виду
W_omy_Un_z = collect(W_omy_Un / (1 - K11*W_omy_Un - K12*W_b_Un - K14*W_sigmn_Un));
[a, b] = numden(W_omy_Un_z);
c = coeffs(b, z);
c = c(end);
a_omy_Un_z = vpa(collect(a / c));
b_omy_Un_z = vpa(collect(b / c));

W_b_Un_z = collect(W_b_Un / (1 - K11*W_omy_Un - K12*W_b_Un - K14*W_sigmn_Un));
[a, b] = numden(W_b_Un_z);
c = coeffs(b, z);
c = c(end);
a_b_Un_z = vpa(collect(a / c));
b_b_Un_z = vpa(collect(b / c));

W_sigmn_Un_z = collect(W_sigmn_Un / (1 - K11*W_omy_Un - K12*W_b_Un - K14*W_sigmn_Un));
[a, b] = numden(W_sigmn_Un_z);
c = coeffs(b, z);
c = c(end);
a_sigmn_Un_z = vpa(collect(a / c));
b_sigmn_Un_z = vpa(collect(b / c));

% Вычислим номинальную передаточную функцию, исходя из требований
T = sqrt(1 / ( 16 + (log(20))^2 ));
psi = log(20) * T;
Wn = (1 / (T^2*Tn*s^3 + (2*Tn*psi*T + T^2)*s^2 + (Tn + 2*T*psi)*s + 1));
Wn = convert_to_tf(Wn, false);
Wn_z = c2d(Wn, T0);

% Найдём коэффициенты обратной связи, проводя сравнение коэффициентов с
% номинальным характеристическим многочленом
A_K = [0.24877, 0.02045, -0.80531;...
       -0.10019, 0.04396, 1.35556;...
       -0.13607, 0.00514, -0.68982];
B_K = [0.80596;...
       -0.67349;...
       0.10054];
K = inv(A_K)*B_K;
K = K';

% Подставим полученные значения коэффициентов в передаточные функции
% замкнутого цифрового контура
b_omy_Un_z = subs(b_omy_Un_z, [K11, K12, K14], [K(1), K(2), K(3)]);
b_b_Un_z = subs(b_b_Un_z, [K11, K12, K14], [K(1), K(2), K(3)]);
b_sigmn_Un_z = subs(b_sigmn_Un_z, [K11, K12, K14], [K(1), K(2), K(3)]);

% Приведём всё к типу данным tf
a_omy_Un_z = sym2poly(a_omy_Un_z);
b_omy_Un_z = sym2poly(b_omy_Un_z);
W_omy_Un_z = tf(a_omy_Un_z, b_omy_Un_z, T0);

a_b_Un_z = sym2poly(a_b_Un_z);
b_b_Un_z = sym2poly(b_b_Un_z);
W_b_Un_z = tf(a_b_Un_z, b_b_Un_z, T0);

a_sigmn_Un_z = sym2poly(a_sigmn_Un_z);
b_sigmn_Un_z = sym2poly(b_sigmn_Un_z);
W_sigmn_Un_z = tf(a_sigmn_Un_z, b_sigmn_Un_z, T0);

% Построим графики переходных функций
d_step(W_b_Un_z)
grid on;
title('wy(t) при ступенчатом воздействии Un');
ylabel('wy(t)', Rotation=0);

d_step(W_omy_Un_z)
grid on;
title('b(t) при ступенчатом воздействии Un');
ylabel('b(t)', Rotation=0);

d_step(W_sigmn_Un_z)
grid on;
title('dn(t) при ступенчатом воздействии Un');
ylabel('dn(t)', Rotation=0);
