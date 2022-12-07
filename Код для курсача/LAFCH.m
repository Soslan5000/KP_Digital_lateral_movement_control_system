function [outputArg] = LAFCH(W, T0)

syms w z lambda

z2w = (1 + w*T0/2)/(1 - w*T0/2);

H = subs(W, z, z2w);
[a, b] = numden(H);
c = coeffs(a, w);
c = max(c);
a_H = vpa(collect(a / c));
b_H = vpa(collect(b / c));

a_H_l = subs(a_H, w, 1i*lambda);
b_H_l = subs(b_H, w, 1i*lambda);

H = a_H_l / b_H_l;
F = matlabFunction(H);

x = logspace(-1,2,1000);
y = F(x);
L = 20*log10(abs(y));

phi = [];
for i = 1:length(x)
    if (real(y(i)) > 0 && imag(y(i)) > 0)
        phi(i) = rad2deg(atan(imag(y(i))/real(y(i))));
    end

    if (real(y(i)) > 0 && imag(y(i)) < 0)
        phi(i) = rad2deg(atan(imag(y(i))/real(y(i))));
    end

    if (real(y(i)) < 0 && imag(y(i)) > 0)
        phi(i) = rad2deg(atan(imag(y(i))/real(y(i)))) + 180;
    end

    if (real(y(i)) < 0 && imag(y(i)) < 0)
        phi(i) = rad2deg(atan(imag(y(i))/real(y(i)))) - 180;
    end
end

figure;
subplot(2, 1, 1)
semilogx(x, L)
xlabel("lg(\lambda)")
ylabel("20 lg |W|, Дб", 'FontSize', 15)
grid on
subplot(2, 1, 2)
semilogx(x, phi)
xlabel("lg(\lambda)")
ylabel("\phi, градусы",'FontSize', 15)
ylim([-240 240])
grid on

outputArg = H;

end