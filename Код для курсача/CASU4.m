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
K = [K11, K12, 0, K14, 0;
     0, 0, 0, 0, 0];

G_raz = G + Dd*K;

W_raz = collect((inv(z*I-G_raz))*Dd);

Wx_raz = W_raz(3, 2);
tf_Wx_raz = convert_to_tf(Wx_raz, true);

figure;
Wx_raz_w = (LAFCH(Wx_raz, T0));

K23 = 0.05;
K = [K11, K12, 0, K14, 0;
     0, 0, K23, 0, 0];
G_zam = G + Dd*K;

W_zam = collect((inv(z*I-G_zam))*Dd);

Wx_zam = W_zam(3, 2);
tf_Wx_zam = convert_to_tf(Wx_zam, true);

Wx_zam_w = (LAFCH(Wx_zam, T0));

figure;
step(d2c(tf_Wx_zam))
title('wx(t) при ступенчатом воздействии Ue');
grid on;
ylabel('wx(t)', Rotation=0);

Wy_zam_un = W_zam(1, 1)
Wy_zam_ue = W_zam(1, 2)
Wb_zam_un = W_zam(2, 1)
Wb_zam_ue = W_zam(2, 2)
Wx_zam_un = W_zam(3, 1)
Wx_zam_ue = W_zam(3, 2)
Wsn_zam_un = W_zam(4, 1)
Wsn_zam_ue = W_zam(4, 2)
Wse_zam_un = W_zam(5, 1)
Wse_zam_ue = W_zam(5, 2)

tf_Wx_zam_un = convert_to_tf(Wx_zam_un, true)
tf_Wx_zam_ue = convert_to_tf(Wx_zam_ue, true)
tf_Wb_zam_un = convert_to_tf(Wb_zam_un, true)
tf_Wb_zam_ue = convert_to_tf(Wb_zam_ue, true)
tf_Wy_zam_un = convert_to_tf(Wy_zam_un, true)
tf_Wy_zam_ue = convert_to_tf(Wy_zam_ue, true)
tf_Wsn_zam_un = convert_to_tf(Wsn_zam_un, true)
tf_Wsn_zam_ue = convert_to_tf(Wsn_zam_ue, true)
tf_Wse_zam_un = convert_to_tf(Wse_zam_un, true)
tf_Wse_zam_ue = convert_to_tf(Wse_zam_ue, true)

figure;
step(d2c(tf_Wy_zam_un))
grid on;
title('wy(t) при ступенчатом воздействии Un');
ylabel('wy(t)', Rotation=0);

figure;
step(d2c(tf_Wy_zam_ue))
grid on;
title('wy(t) при ступенчатом воздействии Ue');
ylabel('wy(t)', Rotation=0);

figure;
step(d2c(tf_Wb_zam_un))
grid on;
title('b(t) при ступенчатом воздействии Un');
ylabel('b(t)', Rotation=0);

figure;
step(d2c(tf_Wb_zam_ue))
grid on;
title('b(t) при ступенчатом воздействии Ue');
ylabel('b(t)', Rotation=0);

figure;
step(d2c(tf_Wx_zam_un))
grid on;
title('wx(t) при ступенчатом воздействии Un');
ylabel('wx(t)', Rotation=0);

figure;
step(d2c(tf_Wx_zam_ue))
grid on;
title('wx(t) при ступенчатом воздействии Ue');
ylabel('wx(t)', Rotation=0);

figure;
step(d2c(tf_Wsn_zam_un))
grid on;
title('sn(t) при ступенчатом воздействии Un');
ylabel('sn(t)', Rotation=0);

figure;
step(d2c(tf_Wsn_zam_ue))
grid on;
title('sn(t) при ступенчатом воздействии Ue');
ylabel('sn(t)', Rotation=0);

figure;
step(d2c(tf_Wse_zam_un))
grid on;
title('se(t) при ступенчатом воздействии Un');
ylabel('se(t)', Rotation=0);

figure;
step(d2c(tf_Wse_zam_ue))
grid on;
title('se(t) при ступенчатом воздействии Ue');
ylabel('se(t)', Rotation=0);

K2 = 2*K;
G_zam = G + Dd*K2;

W_zam = collect((inv(z*I-G_zam))*Dd);

Wx_zam = W_zam(1, 2);
tf_Wx_zam = convert_to_tf(Wx_zam, true);

Wx_zam_w = (LAFCH(Wx_zam, T0));

figure;
step(d2c(tf_Wx_zam))

Wy_zam_un = W_zam(1, 1)
Wy_zam_ue = W_zam(1, 2)
Wb_zam_un = W_zam(2, 1)
Wb_zam_ue = W_zam(2, 2)
Wx_zam_un = W_zam(3, 1)
Wx_zam_ue = W_zam(3, 2)
Wsn_zam_un = W_zam(4, 1)
Wsn_zam_ue = W_zam(4, 2)
Wse_zam_un = W_zam(5, 1)
Wse_zam_ue = W_zam(5, 2)

tf_Wx_zam_un = convert_to_tf(Wx_zam_un, true)
tf_Wx_zam_ue = convert_to_tf(Wx_zam_ue, true)
tf_Wb_zam_un = convert_to_tf(Wb_zam_un, true)
tf_Wb_zam_ue = convert_to_tf(Wb_zam_ue, true)
tf_Wy_zam_un = convert_to_tf(Wy_zam_un, true)
tf_Wy_zam_ue = convert_to_tf(Wy_zam_ue, true)
tf_Wsn_zam_un = convert_to_tf(Wsn_zam_un, true)
tf_Wsn_zam_ue = convert_to_tf(Wsn_zam_ue, true)
tf_Wse_zam_un = convert_to_tf(Wse_zam_un, true)
tf_Wse_zam_ue = convert_to_tf(Wse_zam_ue, true)

figure;
step(d2c(tf_Wy_zam_un))
grid on;
title('wy(t) при ступенчатом воздействии Un');
ylabel('wy(t)', Rotation=0);

figure;
step(d2c(tf_Wy_zam_ue))
grid on;
title('wy(t) при ступенчатом воздействии Ue');
ylabel('wy(t)', Rotation=0);

figure;
step(d2c(tf_Wb_zam_un))
grid on;
title('b(t) при ступенчатом воздействии Un');
ylabel('b(t)', Rotation=0);

figure;
step(d2c(tf_Wb_zam_ue))
grid on;
title('b(t) при ступенчатом воздействии Ue');
ylabel('b(t)', Rotation=0);

figure;
step(d2c(tf_Wx_zam_un))
grid on;
title('wx(t) при ступенчатом воздействии Un');
ylabel('wx(t)', Rotation=0);

figure;
step(d2c(tf_Wx_zam_ue))
grid on;
title('wx(t) при ступенчатом воздействии Ue');
ylabel('wx(t)', Rotation=0);

figure;
step(d2c(tf_Wsn_zam_un))
grid on;
title('sn(t) при ступенчатом воздействии Un');
ylabel('sn(t)', Rotation=0);

figure;
step(d2c(tf_Wsn_zam_ue))
grid on;
title('sn(t) при ступенчатом воздействии Ue');
ylabel('sn(t)', Rotation=0);

figure;
step(d2c(tf_Wse_zam_un))
grid on;
title('se(t) при ступенчатом воздействии Un');
ylabel('se(t)', Rotation=0);

figure;
step(d2c(tf_Wse_zam_ue))
grid on;
title('se(t) при ступенчатом воздействии Ue');
ylabel('se(t)', Rotation=0);

