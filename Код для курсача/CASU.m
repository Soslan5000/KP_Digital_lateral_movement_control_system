clc
clear

syms s

% Зададим матрицы системы уравнений
A = [-0.6, -5.71, -0.04;...
     1, -0.26, 0.065;...
     -0.7, -24, -2.5];

B = [-2.9, 0.55;...
     -0.04, 0;...
     -3.3, -19];

C = [1, 0, 0;...
     0, 1, 0;...
     0, 0, 1];

D = 0;

I = eye(3, 3);

% Получаем матрицу передаточных функций от всех управляющих воздействий ко
% всем наблюдаемым величинам
W_matrix = simplify(C * inv((s*I - A)) * B);

% Выводим передаточную функцию W_wx_dn в стандартном виде отношения полиномов
W_wy_dn = W_matrix(1, 1);
W_wy_dn = convert_to_tf(W_wy_dn, false)

% Выводим передаточную функцию W_b_dn в стандартном виде отношения полиномов
W_b_dn = W_matrix(2, 1);
W_b_dn = convert_to_tf(W_b_dn, false)

% Выводим передаточную функцию W_wy_dn в стандартном виде отношения полиномов
W_wx_dn = W_matrix(3, 1);
W_wx_dn = convert_to_tf(W_wx_dn, false)

% Выводим передаточную функцию W_wx_de в стандартном виде отношения полиномов
W_wy_de = W_matrix(1, 2);
W_wy_de = convert_to_tf(W_wy_de, false)

% Выводим передаточную функцию W_b_de в стандартном виде отношения полиномов
W_b_de = W_matrix(2, 2);
W_b_de = convert_to_tf(W_b_de, false)

% Достаём передаточную функцию W_wy_de в стандартном виде отношения полиномов
W_wx_de = W_matrix(3, 2);
W_wx_de = convert_to_tf(W_wx_de, false)

% Строим графики при помощи созданной нами функции tf_and_step(W)
figure;
step(W_wy_dn);
grid on;
title('wy(t) при ступенчатом воздействии dn');
ylabel('wy(t)', Rotation=0);

figure;
step(W_b_dn);
grid on;
title('b(t) при ступенчатом воздействии dn');
ylabel('b(t)', Rotation=0);

figure;
step(W_wx_dn);
grid on;
title('wx(t) при ступенчатом воздействии dn');
ylabel('wx(t)', Rotation=0);

figure;
step(W_wy_de);
grid on;
title('wy(t) при ступенчатом воздействии de');
ylabel('wy(t)', Rotation=0);

figure;
step(W_b_de);
grid on;
title('b(t) при ступенчатом воздействии de');
ylabel('b(t)', Rotation=0);

figure;
step(W_wx_de);
grid on;
title('wx(t) при ступенчатом воздействии de');
ylabel('wx(t)', Rotation=0);

% Получаем полюса передаточных функций при помощи созданной нами функции
% POLUSES(W)
disp([newline, 'Полюса (корни характеристического многочлена)'])
POLUSES_W_wx_dn = roots(sym2poly(det(s*I - A)))