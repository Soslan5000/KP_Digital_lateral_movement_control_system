clc
clear

syms s z K11 K12 K14 w
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

% Запишем коэффициенты обратной связи
K12 = 2.1617;
K = [K11, 2.1617, K14];
K = 2*K;

% Характеристический многочлен системы
% Вычислим матрицу передаточных функций разомкнутой цифровой системы
W = (inv((z*I - (G + Dd*K))))*Dd;

% Выделим передаточные функции по всем каналам
W_omy_Un_z = collect(W(1));
W_b_Un_z = collect(W(2));
W_sigmn_Un_z = collect(W(3));

% Заменим во всех передаточных функциях z на (1+w)/(1-w)
W_omy_Un_w = subs(W_omy_Un_z, z, (1+w)/(1-w));
[a, b] = numden(W_omy_Un_w);
c = max(coeffs(b));
a_omy_Un_w = vpa(collect(a / c));
b_omy_Un_w = vpa(collect(b / c));

W_b_Un_w = subs(W_b_Un_z, z, (1+w)/(1-w));
[a, b] = numden(W_b_Un_w);
c = max(coeffs(b));
a_b_Un_w = vpa(collect(a / c));
b_b_Un_w = vpa(collect(b / c));

W_sigmn_Un_w = subs(W_sigmn_Un_z, z, (1+w)/(1-w));
[a, b] = numden(W_sigmn_Un_w);
c = max(coeffs(b));
a_sigmn_Un_w = vpa(collect(a / c));
b_sigmn_Un_w = vpa(collect(b / c));

% Вычислим характеристический многочлен и выделим коэффициенты
raus_znam = b_b_Un_w;
coef_raus = coeffs(raus_znam, w);
a0 = coef_raus(1);
a1 = coef_raus(2);
a2 = coef_raus(3);
a3 = coef_raus(4);

% Запишем неравенства ждя критерия Рауса
eq1 = a0 > 0;
eq2 = a1 > 0;
eq3 = a2 > 0;
eq4 = a3 > 0;
eq5 = simplify(a1*a2 - a0*a3) > 0;

% Преобразуем неравенства из символьного типа данных в функции
eq1 = matlabFunction(eq1);
eq2 = matlabFunction(eq2);
eq3 = matlabFunction(eq3);
eq4 = matlabFunction(eq4);
eq5 = matlabFunction(eq5);

% Задание сетки
l11=[-2:0.05:6];
l14=[-2:0.005:2];
[k11, k14] = meshgrid(l11,l14);

% Проверка устойчивости узлов сетки
c1 = eq1(k11,k14);
c2 = eq2(k11,k14);
c3 = eq3(k11,k14);
c4 = eq4(k11,k14);
c5 = eq5(k11,k14);

Kv11 = 1.6127;
Kv14 = -0.4477;
figure;
contourf(l11, l14, c1 & c2 & c3 & c4 & c5, [1 1 1 1 1])
colormap lines
hold on
pnt = scatter(Kv11, Kv14,'r','filled');
hold off
legend([pnt],"Выбранное значение коэффициентов")

xlabel("K11")
ylabel("K14")