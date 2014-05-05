function y = exp_with_no_const(p,x,e)
	% p[1] is norm
	% x is time argument
	y = p(1) * exp(x*e);
end
