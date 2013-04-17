function corr = import_corr(place,name,conf,field)
	% place: dir of the data file
	% name: what kind of correlator
	% conf: all trajactories
	% T: length in T direction
	% field: which field in data file to import

	corr_tmp = [];
	corr = [];
	conf_count=0;
	for traj = conf
		trajs=int2str(traj);
		if ~exist([place,'traj_',trajs,'_',name],'file')
			continue;
		end
		conf_count = conf_count + 1;
		corr_tmp = importdata([place,'traj_', trajs,'_',name]);
		corr = [corr , corr_tmp(:,field)];
	end

	clear corr_tmp;
	clear conf_count;
	clear trajs;
	clear traj;
end
