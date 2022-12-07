function outputArg = tf_and_step(W)

[a, b] = numden(W); % Разбивает дробь на числитель (a) и знаменатель (b)
a = sym2poly(a); % Конвертация числителя в вектор с коэффициентами полинома
b = sym2poly(b); % Конвертация знаменателя в вектор с коэффициентами полинома
transfer_func = tf(a, b); % получение передаточной функции в необходимом формате

% Строим графики переходных процессов по передаточной функции при помощи
% встроенное в матлаб функции step(W)
figure;
step(transfer_func);
grid on;
xlabel('t, c');
ylabel('h(t)', Rotation=0);

outputArg = transfer_func;

end