% data analysis for kaon mass

% initialize fitting parameters:
run('conf.m');

kaon_name = 'kaoncorr';
frozen = 0;% 1 if doing frozen jackknife
correlated = 0; % 1 if doing correlated fitting
cal_effMass = 1;
t_trans = 0:(t_size-1);
fit_range = 7:(t_size+2-7); 

% 1 import data
kaon_corr = tavg_import_corr_bin(place,kaon_name,all_traj,t_trans,3,t_size,bin_size);
% 2 jackknife data
jackknifed_kaoncorr = jackknife(kaon_corr,fit_range,frozen,correlated);
jackknifed_kaoncorr_full = jackknife(kaon_corr,1:t_size,frozen,correlated);
% 3 fit
num_kaoncorr = size(jackknifed_kaoncorr,1);
fprintf('Number of kaon_corr file:\t%d\n',num_kaoncorr);
kaon_result = [];
guess = [0.3, 1e+6];
for i = 1:num_kaoncorr
	kaon_result = [kaon_result; abs(fit(jackknifed_kaoncorr{i,1},jackknifed_kaoncorr{i,2},@(x1,x2,x3) cosh_no_const(x1,x2,x3),guess,t_size))];
end
% 4 effective mass
if cal_effMass
	eff_mass(jackknifed_kaoncorr_full,  t_size, 0.,  num_kaoncorr,  'kaon_Meff');
end

% 5 show results
fprintf('kaon mass:\t%.6f\t std:\t %.6f\n',mean(kaon_result(:,1)),std(kaon_result(:,1)) * sqrt(num_kaoncorr-1));
