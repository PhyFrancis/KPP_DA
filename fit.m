function y = fit(data, cinv, fit_func, guess, T) 
	% guess[1] = mass
	% guess[2] = a (normalization factor)
	% guess[3] = const if exist
	p = ...
		fminsearch(@(p) chi_sq(p, cinv, data,fit_func, T), ...
		guess, ...
		optimset('MaxFunEvals', 100000000,'MaxIter', 100000,'TolFun',1e-12,'TolX',1e-12));
	y = p;
end

function chi2 = chi_sq(p, cov_inv , data, fit_func, T)
	% p[1] = mass
	% p[2] = a (normalization factor)
	% data stors the "time \t corr \n"
	if size(data,1) ~= size(cov_inv,2)
		disp('The dimension of data and cov_inv do not match!');
		exit;
	end

	var = [];
	length = size(data,1);
	for i = 1:length
		var = [var, fit_func(p,data(i,1),T) - data(i,2)];
	end
	chi2 = var * cov_inv * var';
end
