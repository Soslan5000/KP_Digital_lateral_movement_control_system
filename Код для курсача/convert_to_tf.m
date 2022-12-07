function [outputArg] = convert_to_tf(W_in, is_z_func)
% Преобразование дроби к типу данных tf
T0 = 0.18;
[a, b] = numden(W_in);
a = sym2poly(a);
b = sym2poly(b);
a = a / b(end);
b = b / b(end);
if is_z_func
    W_out = tf(a, b, T0);
else
    W_out = tf(a, b);
end

outputArg = W_out;

end