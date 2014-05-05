function eff = eff_mass(jackknifed_data, t_size, const) 
  % FIXME how to use the const parameter? Now just treat it as 0.
	eff = [];
	eff_all = [];
	if size(jackknifed_data{1,1},1) ~= t_size
		fprintf('Should use full T-range for effective mass!');
	end
	n_conf = size(jackknifed_data,1);

	for i = 1:n_conf
		eff_tmp = [];
		for t = 1:t_size-1
			obj = (jackknifed_data{i,1}(t+1,2)-const(i))...
				/ (jackknifed_data{i,1}(t,2)-const(i));
			ratio = @(x) (exp(-x*t)+exp(-x*(t_size-t)))...
				/ (exp(-x*(t-1))+exp(-x*(t_size-(t-1)))) - obj;
			eff_tmp = [eff_tmp; fsolve(ratio,0.3)];
		end
		eff_all = [eff_all, eff_tmp];
	end

	for t = 1:t_size-1
		eff = [eff; [t, mean(eff_all(t,:)), std(eff_all(t,:))*sqrt(n_conf-1)]];
	end

end
