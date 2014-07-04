function y = zeta_func(t, q, gn2)
%define the integration fuction part, and gn2=(Gamma*n)^2
y = exp(t*q^2-pi^2*gn2./t).*(pi./t).^(1.5);
