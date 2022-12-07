clc
clear

syms s z L1 L2 L4 lambd

Tn = 0.11;

% Зададим матрицы системы уравнений
A = [-0.26, 0.065, -0.04;...
     -24, -2.5, -3.3;...
     0, 0, -1/Tn];

B = [0;...
     0;...
     1/Tn];

C = [1, 0, 0];

D = 0;

I = eye(3, 3);

K = [10.2403, 0.6995, -0.3973];

L = [L1;...
     L2;...
     L4];

% Вычислим характеристический многочлен замкнутого фильтра Калмана
Dz = vpa(collect(det(lambd*I - (A + B*K + L*C)), lambd));
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
b0 = TF_DG.denominator{1}(4);
b1 = TF_DG.denominator{1}(3);
b2 = TF_DG.denominator{1}(2);
b3 = TF_DG.denominator{1}(1);


eq0 = a0 == b0;
eq1 = a1 == b1;
eq2 = a2 == b2;
root = solve(eq0, eq1, eq2, L1, L2, L4);

L(1) = root.L1;
L(2) = root.L2;
L(3) = root.L4;

Q = [A + B*K, [0, 0, 0; 0, 0, 0; 0, 0, 0];...
     A + B*K + L*C, -L*C];
R = [1, 0, 0, 1, 0, 0];
P = [B;
     0;
     0;
     0];
W = inv(s*eye(6) - Q) * P;
W_1 = convert_to_tf(vpa(collect(W(1))), false)
W_2 = convert_to_tf(vpa(collect(W(2))), false)
W_3 = convert_to_tf(vpa(collect(W(3))), false)
W_4 = convert_to_tf(vpa(collect(W(4))), false)

W_5 = vpa(collect(W(5))
[a, b] = numden(vpa(collect(W(6))));
a = sym2poly(a);
b = sym2poly(b);
a = a / b(1);
b = b / b(1);
W_5 = convert_to_tf(W_5, false)

W_6 = convert_to_tf(vpa(collect(W(6))), false)
[a, b] = numden(vpa(collect(W(6))));
a = sym2poly(a);
b = sym2poly(b);
a = a / b(1);
b = b / b(1);


