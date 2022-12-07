function [outputArg] = tf2sym(TF)
syms x

den = vpa(poly2sym(TF.denominator{1}));
num = vpa(poly2sym(TF.numerator{1}));

W = num / den;

outputArg = W;
end