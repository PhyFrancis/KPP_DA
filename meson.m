%data analysis for meson mass

place = '/home/daiqian/A2Adata/L100_11_S1_21_Mu_0.010_Ms_0.099/';
frozen = 1;% 1 if doing frozen jackknife
correlated = 0; % 1 if doing correlated fitting
kaon_name = 'kaoncorr_1.00';
pion_name = 'pioncorr';

% fitting parameters:
t_size = 32;
all_traj = 1000:20:2140;
field = 2; % 2nd colume of corr file is the real part of correlator
fit_range = 11:22;

%1 kaon calculation
%1.1 import data
kaon_corr = import_corr(place,kaon_name,all_traj,field);
%1.2 jackknife data
jackknifed_kaoncorr = jackknife(kaon_corr,fit_range,frozen,correlated);
%1.3 fit
num_kaoncorr = size(jackknifed_kaoncorr,1);
fprintf('Number of kaon_corr file:\t%d\n',num_kaoncorr);
kaon_result = [];
guess = [0.2, 1e+6];
for i = 1:num_kaoncorr
	kaon_result = [kaon_result; abs(fit(jackknifed_kaoncorr{i,1},jackknifed_kaoncorr{i,2},@(x1,x2,x3) cosh_no_const(x1,x2,x3),guess,t_size))];
end

%2 pion calculation
%2.1 import data
pion_corr = import_corr(place,pion_name,all_traj,field);
%2.2 jackknife data
jackknifed_pioncorr = jackknife(pion_corr,fit_range,frozen,correlated);
%2.3 fit
num_pioncorr = size(jackknifed_pioncorr,1);
fprintf('Number of pion_corr file:\t%d\n',num_pioncorr);
pion_result = [];
guess = [0.2, 1e+6];
for i = 1:num_pioncorr
	pion_result = [pion_result; abs(fit(jackknifed_pioncorr{i,1},jackknifed_pioncorr{i,2},@(x1,x2,x3) cosh_no_const(x1,x2,x3),guess,t_size))];
end

%3 pipi calculation
%3.1 import data
sep = 2;
FigureC_name = ['FigureC_sep',int2str(sep)];
FigureD_name = ['FigureD_sep',int2str(sep)];
FigureR_name = ['FigureR_sep',int2str(sep)];
FigureVdis_name = ['FigureVdis_sep',int2str(sep)];
FigureC = import_corr(place,FigureC_name,all_traj,field);
FigureD = import_corr(place,FigureD_name,all_traj,field);
FigureR = import_corr(place,FigureR_name,all_traj,field);
FigureVdis = import_corr(place,FigureVdis_name,all_traj,field);
I2 = 2 * (FigureD - FigureC);
I0V = 2 * FigureD + FigureC - 6 * FigureR;
%3.2 jackknife data
fit_range = 8:14;
jackknifed_I2 = jackknife(I2,fit_range,frozen,correlated);
fit_range = 7:14;
jackknifed_I0V = jackknife(I0V,fit_range,frozen,correlated);
%3.3 fit
num_I2 = size(jackknifed_I2,1);
num_I0V = size(jackknifed_I0V,1);
fprintf('Number of I2 pipi file:\t%d\n',num_I2);
fprintf('Number of I0V pipi file:\t%d\n',num_I0V);
I2_result = [];
I0V_result = [];
guess = [0.5, I2(1,1), I2(t_size/2,1)];
for i = 1:num_I2
	I2_result = [I2_result; abs(fit(jackknifed_I2{i,1},jackknifed_I2{i,2},@(x1,x2,x3) cosh_with_const(x1,x2,x3),guess,t_size-2*sep))];
end
guess = [0.5, I0V(1,1), I0V(t_size/2,1)];
for i = 1:num_I0V
	I0V_result = [I0V_result; abs(fit(jackknifed_I0V{i,1},jackknifed_I0V{i,2},@(x1,x2,x3) cosh_with_const(x1,x2,x3),guess,t_size-2*sep))];
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
fprintf('pipi_I2 energy:\t%.6f\t std:\t %.6f\n',mean(I2_result(:,1)),std(I2_result(:,1)) * sqrt(num_I2-1));
fprintf('pipi_I0V energy:\t%.6f\t std:\t %.6f\n',mean(I0V_result(:,1)),std(I0V_result(:,1)) * sqrt(num_I0V-1));
