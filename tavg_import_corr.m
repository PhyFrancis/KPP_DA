function corr = tavg_import_corr(place,name,conf,t_range,field,T)
	% place: dir of the data file
	% name: what kind of correlator
	% conf: all trajactories
	% t_range: the source positions to do time average
	% field: which field in data file to import

	% check whether t_range is within T
	for i = 1:numel(t_range)
		if ~any(t_range(i)==0:T-1)
			fprintf('t_range is not valid')
		end
	end
	% end of check

	corr_tmp = [];
	corr = [];
	for traj = conf
		trajs=int2str(traj);
		if ~exist([place,'traj_',trajs,'_',name],'file')
			fprintf('File %s doesnt exist!\n', ['traj_',trajs,'_',name]);
			continue;
		end
		corr_tmp1 = importdata([place,'traj_', trajs,'_',name]);
		
	  corr_tmp2 = zeros(T,size(field,2));
		nTslice = 0;
		for i = 0:T-1
			if ~any(i==t_range)
				continue;
			end
			if sum(sum(abs(corr_tmp1(i*T+1:(i+1)*T,field)))) == 0  % to rule out those un-computed tslice in Type1-4 contractions
				continue;
			end
			corr_tmp2 = corr_tmp2 + corr_tmp1(i*T+1:(i+1)*T,field);
			nTslice = nTslice + 1;
		end
		corr_tmp2 = corr_tmp2 / nTslice;
		corr = [corr, corr_tmp2];
	end

	clear corr_tmp1;
	clear corr_tmp2;
	clear nTslice;
	clear trajs;
	clear traj;
end
