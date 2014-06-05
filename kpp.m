% This file does K->PiPi(I=0) analysis

place = '/home/daiqian/BGQ/32nt_kpp/L900_11_S0_11_Mu_0.000_Ms_0.045/';
all_traj = [302:4:494];
% place = '/home/daiqian/BGQ/32nt_kpp/L500_11_S0_11_Mu_0.000_Ms_0.045/';
% all_traj = [328,336,342:8:398];
% place = '/home/daiqian/BGQ/32nt_kpp/L500_L900_avg/';
% all_traj = [342:8:398];

%data analysis for meson mass
run('pipi_sep.m');

%data analysis for K -> PiPi (I=0)
%1 initialize fitting parameters:
frozen = 0;% 1 if doing frozen jackknife
correlated = 0; % 1 if doing correlated fitting
cal_effMass = 0;
t_size = 64;
t_trans = 0:63; 
sep = 4;
DeltaT = 12:2:18;
FigureVdis_name = ['FigureVdis_sep',int2str(sep)];
bin_size = 1;

Qi_result_glb = [];
for deltat = DeltaT

	fit_range = 4:(deltat-4);
	Type1_name = ['type1','_deltat_',int2str(deltat),'_sep_',int2str(sep)];
	Type2_name = ['type2','_deltat_',int2str(deltat),'_sep_',int2str(sep)];
	Type3_name = ['type3','_deltat_',int2str(deltat),'_sep_',int2str(sep)];
	Type4_name = ['type4'                                                ];

	%2 I0 kpipi calculation

	%2.1 import data
	FigureVdis = import_corr_bin(place,FigureVdis_name,all_traj,2,bin_size);
	Type1 = tavg_import_corr_bin(place,Type1_name,all_traj,t_trans,2*[1:48]+1,t_size,bin_size);
	Type2 = tavg_import_corr_bin(place,Type2_name,all_traj,t_trans,2*[1:48]+1,t_size,bin_size) * 0.5; % 0.5 because mom permutation
	Type3 = tavg_import_corr_bin(place,Type3_name,all_traj,t_trans,2*[1:80]+1,t_size,bin_size) * 0.5; % 0.5 because mom permutation
	Type4_kqi = tavg_import_corr_bin(place,Type4_name,all_traj,t_trans,2*[1:80]+1,t_size,bin_size); % tavg'ed < Qi K >
	Type4_kqi_full = import_corr_bin(place,Type4_name,all_traj,2*[1:80]+1,bin_size); % not tavg'ed < Qi K >
	Type3S5D = tavg_import_corr_bin(place,Type3_name,all_traj,t_trans,2*(80+[1:2])+1,t_size,bin_size); % tavg'ed < Pi Pi s5d K > < Pi Pi u5s' K >
	Type4S5D_kqi = tavg_import_corr_bin(place,Type4_name,all_traj,t_trans,2*(80+[1:2])+1,t_size,bin_size); % tavg'ed < s5d K > < u5s' K >
	Type4S5D_kqi_full = import_corr_bin(place,Type4_name,all_traj,2*(80+[1:2])+1,bin_size); % not tavg'ed < s5d K > < u5s' K >

	fprintf('Number of FigureVdis file:\t%d\n',size(FigureVdis,2));
	nType1 = size(Type1    ,2)/48; fprintf('Number of type1 file:\t%d\n',nType1);
	nType2 = size(Type2    ,2)/48; fprintf('Number of type2 file:\t%d\n',nType2);
	nType3 = size(Type3    ,2)/80; fprintf('Number of type3 file:\t%d\n',nType3);
	nType4 = size(Type4_kqi,2)/80; fprintf('Number of type4 file:\t%d\n',nType4);

	if isequaln(nType1, nType2, nType3, nType4) == 0
		printf('Number of Type1, Type2, Type3, Type4 files dont match!\n');
	end

	%2.2 combine data
	Qi_type1 = combine_type1(Type1);
	Qi_type2 = combine_type2(Type2);
	Qi_type3 = combine_type3(Type3);
	Qi_type4_kqi = combine_type4(Type4_kqi); % tavg'ed < Qi K >
	Qi_type4_kqi_full = combine_type4(Type4_kqi_full); % not tavg'ed < Qi K >
	Qi_type4 = mult_pipi(Qi_type4_kqi_full, - FigureVdis + mean2(FigureVdis), t_size, t_trans, deltat+sep); % tavg'ed < Pi Pi Qi K > - < Pi Pi > < Qi K >
	type4S5D = mult_pipi_single(Type4S5D_kqi_full(:,2*[1:nType4]-1), - FigureVdis + mean2(FigureVdis), t_size, t_trans, deltat+sep); % tavg'ed < Pi Pi s5d K > - < Pi Pi > < s5d K >
	type3S5D = Type3S5D(:,2*[1:nType4]-1); % tavg'ed < Pi Pi s5d K > 

	%2.3 calculate \alpha = <0| Qi K> / <0| s5d K >
	jackknifed_type4S5D_kqi = jackknife(Type4S5D_kqi(:,2*[1:nType4]-1),fit_range,frozen,correlated); % odd columns of Type4S5D: only consider < s5d K >
	jackknifed_type4_kqi = cell(10,1);
	alpha = cell(10,1);
	for Q = 1:10 
		jackknifed_type4_kqi{Q,1} = jackknife(Qi_type4_kqi{Q,1},fit_range,frozen,correlated);
		alpha{Q,1} = [];
		for conf = 1:nType4
			alpha{Q,1} = [alpha{Q,1}, jackknifed_type4_kqi{Q,1}{conf,1}(:,2) ./ jackknifed_type4S5D_kqi{conf,1}(:,2)];
		end
	end

	alpha_out = [];
	for Q = 1:10
		alpha_out = [alpha_out, mean(alpha{Q,1},2), std(alpha{Q,1},0,2) * sqrt(nType4-1)];
	end
	csvwrite('alpha', [fit_range',alpha_out]);
	system(['./change_csv.sh ','alpha'] ,'-echo');

	%2.4 jackknife data and subtract mixed graph
	jackknifed_type3S5D = jackknife(type3S5D,fit_range,frozen,correlated);
	jackknifed_type4S5D = jackknife(type4S5D,fit_range,frozen,correlated);
	Qi = cell(10,1);
	jackknifed_Qi = cell(10,1);
	for Q = 1:10

		Qi{Q,1} = Qi_type1{Q,1} + Qi_type2{Q,1} + Qi_type3{Q,1} + Qi_type4{Q,1};
		jackknifed_Qi{Q,1} = jackknife(Qi{Q,1},fit_range,frozen,correlated);
		for conf = 1:nType4 
			jackknifed_Qi{Q,1}{conf,1}(:,2) = jackknifed_Qi{Q,1}{conf,1}(:,2) - ...
				alpha{Q,1}(:,conf) .* (jackknifed_type3S5D{conf,1}(:,2) + jackknifed_type4S5D{conf,1}(:,2));

			% jackknifed_Qi{Q,1}{conf,1}(:,2) = jackknifed_Qi{Q,1}{conf,1}(:,2);
		end

		% Qi{Q,1} = Qi_type1{Q,1} + Qi_type2{Q,1} + Qi_type3{Q,1};
		% jackknifed_Qi{Q,1} = jackknife(Qi{Q,1},fit_range,frozen,correlated);
		% for conf = 1:nType4 
		% 	jackknifed_Qi{Q,1}{conf,1}(:,2) = jackknifed_Qi{Q,1}{conf,1}(:,2) - ...
		% 		alpha{Q,1}(:,conf) .* (jackknifed_type3S5D{conf,1}(:,2));
		% end

		% Qi{Q,1} = Qi_type1{Q,1} + Qi_type2{Q,1}; 
		% jackknifed_Qi{Q,1} = jackknife(Qi{Q,1},fit_range,frozen,correlated);

	end

	%2.5 fit Qi
	Qi_result = cell(10,1);
	for Q = 1:10
		for conf = 1:nType4
			deltaE = I0_result(conf,1) - kaon_result(conf,1);
			zPiPi = I0_result(conf,2);
			zK = kaon_result(conf,2);
			zTotal = (zPiPi * zK)^0.5 * exp(- I0_result(conf,1) * deltat) / sqrt(2); % sqrt(2) because only 1/2 of the kaon is contracting
			guess = jackknifed_Qi{Q,1}{conf,1}(1,2);
			Qi_result{Q,1} = [Qi_result{Q,1}; fit(jackknifed_Qi{Q,1}{conf,1},jackknifed_Qi{Q,1}{conf,2},@(x1,x2,x3) exp_with_no_const(x1,x2,x3),guess,deltaE) / zTotal];
		end
	end

	%2.6 display result
	for Q = 1:10
		fprintf('Q%d(I=0):\t%.6e\t std:\t %.6e\n',Q,mean(Qi_result{Q,1}),std(Qi_result{Q,1}) * sqrt(nType4-1));
		Qi_result_glb = [Qi_result_glb; mean(Qi_result{Q,1}), std(Qi_result{Q,1})];
	end

	fn = ['Qi_result/Qi_result_deltat_',int2str(deltat),'_sep_',int2str(sep)];
	save(fn,'Qi_result');
	
end


fprintf('EWA:\n');
% error-weighted avg
err_avg = [];
for Q = 1:10
	fac = 0;
	err_avg_tmp = zeros(1,2);
	for ndelta = 0:(size(DeltaT,2)-1)
		err_avg_tmp = err_avg_tmp + Qi_result_glb(Q+10*ndelta,:) * Qi_result_glb(Q+10*ndelta,2)^(-2);
		fac = fac + Qi_result_glb(Q+10*ndelta,2)^(-2);
	end
	err_avg = [err_avg; err_avg_tmp / fac ];
	fprintf('Q%d(I=0):\t%.6e\t std:\t %.6e\n',Q,err_avg(Q,1),err_avg(Q,2) * sqrt(nType4-1));
end

%% display partial results
%avg_type1 = [];
%avg_type12 = [];
%for i = 1:10
%	avg_type1 = [avg_type1, mean(Qi_type1{i,1},2), std(Qi_type1{i,1},1,2)];
%	avg_type12 = [avg_type12, mean(Qi_type1{i,1}+Qi_type2{i,1},2), std(Qi_type1{i,1}+Qi_type2{i,1},1,2)];
%end
%
%fn_avg_type1  = ['type1_deltat',num2str(deltat)];
%fn_avg_type12 = ['type12_deltat',num2str(deltat)];
%csvwrite(fn_avg_type1 ,[[0:63]',avg_type1]);
%csvwrite(fn_avg_type12,[[0:63]',avg_type12]);
%system(['./change_csv.sh ',fn_avg_type1] ,'-echo');
%system(['./change_csv.sh ',fn_avg_type12] ,'-echo');
