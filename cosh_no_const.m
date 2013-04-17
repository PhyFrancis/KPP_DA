function y = cosh_no_const(p,x,T)
	% p[1] is mass
	% p[2] is normalization factor
	% x is time argument
	% T is the total time distance
	y = p(2) * (exp(-p(1)*x) + exp(-p(1)*(T-x)));
end
