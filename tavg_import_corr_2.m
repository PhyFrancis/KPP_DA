function corr = tavg_import_corr_2(place,name,conf,t_range,field,T)
	% place: dir of the data file
	% name: what kind of correlator
	% conf: all trajactories
	% t_range: the source positions to do time average ( actually a sum)
	% field: which field in data file to import

	% check whether t_range is within T
	for i = 1:numel(t_range)
		if ~any(t_range(i)==0:T-1)
			fprintf('t_range is not valid')
			exit;
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
		corr_tmp2 = [];
		LN = size(corr_tmp1,1);
		for i = 0:T-1
		  if ~any(i==t_range)
				continue;
			end
			corr_tmp2 = [corr_tmp2, corr_tmp1(i*T+1:(i+1)*T,field)];
		end
		corr = [corr , sum(corr_tmp2, 2)];
	end

	clear corr_tmp1;
	clear corr_tmp2;
	clear trajs;
	clear traj;
end
