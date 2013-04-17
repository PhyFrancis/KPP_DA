%data analysis for meson mass

place = '/home/daiqian/A2Adata/L100_11_S1_21_Mu_0.010_Ms_0.099/';
frozen = 0;% 1 if doing frozen jackknife
correlated = 0; % 1 if doing correlated fitting
kaon_name = 'kaoncorr_1.00';
pion_name = 'pioncorr';

% fitting parameters:
t_size = 32;
all_traj = 1000:20:2140;
field = 2; % 2nd colume of corr file is the real part of correlator
fit_range = 11:23;

%1 kaon calculation
%1.1 import data
kaon_corr = import_corr(place,kaon_name,all_traj,field);
%1.2 jackknife data
jackknifed_kaoncorr = jackknife(kaon_corr,fit_range,frozen,correlated);
num_kaoncorr = size(jackknifed_kaoncorr,1);
fprintf('Number of kaon_corr file:\t%d\n',num_kaoncorr);
%1.3 fit
kaon_result = [];
guess = [0.2, 1e+6];
for i = 1:size(jackknifed_kaoncorr,1)
	kaon_result = [kaon_result; abs(fit(jackknifed_kaoncorr{i,1},jackknifed_kaoncorr{i,2},@(x1,x2,x3) cosh_no_const(x1,x2,x3),guess,t_size))];
end

%2 pion calculation
%2.1 import data
pion_corr = import_corr(place,pion_name,all_traj,field);
%2.2 jackknife data
jackknifed_pioncorr = jackknife(pion_corr,fit_range,frozen,correlated);
num_pioncorr = size(jackknifed_pioncorr,1);
fprintf('Number of pion_corr file:\t%d\n',num_pioncorr);
%2.3 fit
pion_result = [];
guess = [0.2, 1e+6];
for i = 1:size(jackknifed_pioncorr,1)
	pion_result = [pion_result; abs(fit(jackknifed_pioncorr{i,1},jackknifed_pioncorr{i,2},@(x1,x2,x3) cosh_no_const(x1,x2,x3),guess,t_size))];
end

if(frozen)
	disp('frozen');
else 
	disp('unfrozen');
end
if(correlated)
	disp('correlated');
else 
	disp('uncorrelated');
end
fprintf('kaon mass:\t%.6f\t std:\t %.6f\n',mean(kaon_result(:,1)),std(kaon_result(:,1)) * sqrt(num_kaoncorr-1));
fprintf('pion mass:\t%.6f\t std:\t %.6f\n',mean(pion_result(:,1)),std(pion_result(:,1)) * sqrt(num_pioncorr-1));
