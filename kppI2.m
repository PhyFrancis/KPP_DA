% This file does K->PiPi(I=0) analysis

% initialize parameters
run('conf.m');

% analysis meson mass
run('pipi_I2_sep.m');

%data analysis for K -> PiPi (I=0)
%1 initialize fitting parameters:
frozen = 0;% 1 if doing frozen jackknife
correlated = 0; % 1 if doing correlated fitting
cal_effMass = 0;
t_trans = 0:(t_size-1); 
sep = 4;
DeltaT = 12:2:12;

Qi_result_glb = [];
for deltat = DeltaT

	fit_range = 5:(deltat-3);
	% Type1_name = ['type1','_deltat_',int2str(deltat),'_sep_',int2str(sep)];
	Type1_name = ['type1','_deltat_',int2str(deltat),'_sep_',int2str(sep)];

	%2 I2 kpipi calculation

	%2.1 import data
	Type1 = tavg_import_corr_bin(place,Type1_name,all_traj,t_trans,2*[1:48]+1,t_size,bin_size);
	nType1 = size(Type1,2)/48; fprintf('Number of type1 file:\t%d\n',nType1);

	%2.2 combine data
	Qi = combine_type1_I2(Type1);

	%2.4 jackknife data and subtract mixed graph
	jackknifed_Qi = cell(10,1);
	for Q = [1,7,8]
		jackknifed_Qi{Q,1} = jackknife(Qi{Q,1},fit_range,frozen,correlated);
	end

	%2.5 fit Qi
	Qi_result = cell(10,1);
	for Q = [1,7,8]
		for conf = 1:nType1
			deltaE = I2_result(conf,1) - kaon_result(conf,1);
			zPiPi = I2_result(conf,2);
			zK = kaon_result(conf,2);
			zTotal = (zPiPi * zK)^0.5 * exp(- I2_result(conf,1) * deltat) / sqrt(2); % sqrt(2) because only 1/2 of the kaon is contracting
			guess = jackknifed_Qi{Q,1}{conf,1}(1,2);
			Qi_result{Q,1} = [Qi_result{Q,1}; fit(jackknifed_Qi{Q,1}{conf,1},jackknifed_Qi{Q,1}{conf,2},@(x1,x2,x3) exp_with_no_const(x1,x2,x3),guess,deltaE) / zTotal];
		end
	end

	%2.5 fit each contraction
	% jackknifed_Ci = cell(48,1);
	% for C = 1:48
	% 	jackknifed_Ci{C,1} = jackknife(Type1(:,C+[0:nType1-1]*48),fit_range,frozen,correlated);
	% end
	% Ci_result = cell(48,1);
	% for C = 1:48
	% 	for conf = 1:nType1
	% 		deltaE = I2_result(conf,1) - kaon_result(conf,1);
	% 		zPiPi = I2_result(conf,2);
	% 		zK = kaon_result(conf,2);
	% 		zTotal = (zPiPi * zK)^0.5 * exp(- I2_result(conf,1) * deltat) / sqrt(2); % sqrt(2) because only 1/2 of the kaon is contracting
	% 		guess = jackknifed_Ci{C,1}{conf,1}(1,2);
	% 		Ci_result{C,1} = [Ci_result{C,1}; fit(jackknifed_Ci{C,1}{conf,1},jackknifed_Ci{C,1}{conf,2},@(x1,x2,x3) exp_with_no_const(x1,x2,x3),guess,deltaE) / zTotal];
	% 	end
	% end
	% for C = 1:48
	% 	fprintf('C%d(I=2):\t%.6e\t std:\t %.6e\n',C,mean(Ci_result{C,1}),std(Ci_result{C,1})*(nType1-1)^0.5);
	% end

	%2.6 display result
	for Q = [1,7,8]
		fprintf('Q%d(I=2):\t%.6e\t std:\t %.6e\n',Q,mean(Qi_result{Q,1}),std(Qi_result{Q,1})*(nType1-1)^0.5);
	end

	fn = ['Qi_result/Qi_I2_result_deltat_',int2str(deltat),'_sep_',int2str(sep)];
	save(fn,'Qi_result');

	% calculate each operator
	Qi_avg = zeros(size(jackknifed_Qi{1,1}{1,1},1),20);
	for Q = [1,7,8]
		Qi_avg_tmp = [];
		for conf = 1:nType1
			deltaE = I2_result(conf,1) - kaon_result(conf,1);
			zPiPi = I2_result(conf,2);
			zK = kaon_result(conf,2);
			zTotal = (zPiPi * zK)^0.5 * exp(- I2_result(conf,1) * deltat) / sqrt(2); % sqrt(2) because only 1/2 of the kaon is contracting
			Qi_avg_tmp = [Qi_avg_tmp, jackknifed_Qi{Q,1}{conf,1}(:,2) / zTotal];
		end
		Qi_avg(:,2*Q-1) = mean(Qi_avg_tmp,2);
		Qi_avg(:,2*Q)   = std(Qi_avg_tmp,0,2)*(nType1-1)^0.5; 
	end
	fn = ['Qi_result/Qi_I2_avg_delta_',int2str(deltat),'_sep_',int2str(sep)];
	csvwrite(fn, [(fit_range-1)',Qi_avg]);
	system(['./change_csv.sh ',fn] ,'-echo');
	
end
