function jackknifed_data = jackknife(all_data,fit_range,frozen,correlated)
	% frozen is true when doing a frozen jackknife, otherwise false
	% correlated is true when doing correlated fitting, otherwise false

	count = size(all_data,2);
	jackknifed_data = cell(count,2); % first column stores the correlator, second colume is the cov_inv
	if correlated
		frozen_cinv = inv(cov(all_data(fit_range,:)'));
	else
		frozen_cinv = inv(diag(std(all_data(fit_range,:),0,2).^2));
	end

	for i = 1:count

		tmp = [];
		for j = 1:count
			if j==i
				continue; % omitting i's piece
			end
			tmp = [tmp, all_data(fit_range, j)];
		end
		jackknifed_data{i,1} = [(fit_range-1)',mean(tmp,2)]; % -1 because this should be the real t value

		if frozen
			jackknifed_data{i,2} = frozen_cinv;
		else
			if correlated
				jackknifed_data{i,2} = inv(cov(tmp'));
			else
				jackknifed_data{i,2} = inv(diag(mean(tmp,2)));
			end
		end

	end

clear count;
clear frozen_cinv;
clear tmp;
end
