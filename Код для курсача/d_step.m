function [outputArg] = d_step(W)
% Построение графика переходного процесса по заданной передаточной функции
[y,t] = step(W);
step(W)

plot(t, y)
grid on;
xlabel('t, c');
ylabel('h(t)', Rotation=0);

outputArg = W;

end