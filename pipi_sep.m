%data analysis for meson spectrum

%1 initialize fitting parameters:
place = '/home/daiqian/BGQ/32nt_kpp/L500_11_S0_11_Mu_0.000_Ms_0.045/';
frozen = 0;% 1 if doing frozen jackknife
correlated = 0; % 1 if doing correlated fitting
cal_effMass = 1;
all_traj = [320:8:336,342:8:400];
t_size = 64;
t_trans = 0:63;
bin_size = 1;
sep = 4;
do_kaon = 1;
pion_name = 'pioncorr';
kaon_name = 'kaoncorr';
FigureC_name = ['FigureC','_sep',int2str(sep)];
FigureD_name = ['FigureD','_sep',int2str(sep)];
FigureR_name = ['FigureR','_sep',int2str(sep)];
FigureVdis_name = ['FigureVdis','_sep',int2str(sep)];
field = 3; % 3rd colume of corr file is the real part of correlator
pion_fit_range = 7:59; 
kaon_fit_range = 7:59; 
I2_fit_range = [8:(32-sep),(34-sep):(58-2*sep)];
I0V_fit_range = [8:(32-sep),(34-sep):(58-2*sep)];
I0_fit_range = [7:(28-sep),(38-sep):(59-2*sep)];

%2 pion calculation
%2.1 import data
pion_corr = tavg_import_corr_bin(place,pion_name,all_traj,t_trans,field,t_size,bin_size);
%2.2 jackknife data
jackknifed_pioncorr = jackknife(pion_corr,pion_fit_range,frozen,correlated);
jackknifed_pioncorr_full = jackknife(pion_corr,1:t_size,frozen,correlated);
%2.3 fit
num_pioncorr = size(jackknifed_pioncorr,1);
fprintf('Number of pion_corr file:\t%d\n',num_pioncorr);
pion_result = [];
guess = [0.3, 1e+6];
for i = 1:num_pioncorr
	pion_result = [pion_result; abs(fit(jackknifed_pioncorr{i,1},jackknifed_pioncorr{i,2},@(x1,x2,x3) cosh_no_const(x1,x2,x3),guess,t_size))];
end
%2.4 effective mass
if cal_effMass
	eff_mass(jackknifed_pioncorr_full, t_size, 0., num_pioncorr, 'pion_Meff');
end

%3 pipi calculation
%3.1 import data
FigureC = tavg_import_corr_bin(place,FigureC_name,all_traj,t_trans,field,t_size,bin_size);
FigureD = tavg_import_corr_bin(place,FigureD_name,all_traj,t_trans,field,t_size,bin_size);
FigureR = tavg_import_corr_bin(place,FigureR_name,all_traj,t_trans,field,t_size,bin_size);
FigureVdis = import_corr_bin(place,FigureVdis_name,all_traj,2,bin_size);
FigureV_A = cal_FigureV_A(FigureVdis,sep);
FigureV_B = cal_FigureV_B(FigureVdis,sep); % This is the term to be subtracted
I2 = 2 * (FigureD - FigureC);
I0V = 2 * FigureD + FigureC - 6 * FigureR;
I0 = 2 * FigureD + FigureC - 6 * FigureR + 3 * FigureV_A;
jackknifed_I2 = jackknife(I2,I2_fit_range,frozen,correlated);
jackknifed_I2_full = jackknife(I2,1:t_size,frozen,correlated);
jackknifed_I0V = jackknife(I0V,I0V_fit_range,frozen,correlated);
jackknifed_I0V_full = jackknife(I0V,1:t_size,frozen,correlated);
jackknifed_I0 = jackknife(I0,I0_fit_range,frozen,correlated);
jackknifed_I0_subtracted = vac_sub(jackknifed_I0,3*FigureV_B);
jackknifed_I0_full = jackknife(I0,1:t_size,frozen,correlated);
jackknifed_I0_full_subtracted = vac_sub(jackknifed_I0_full,3*FigureV_B);
%%3.3.1 fit I2
num_I2 = size(jackknifed_I2,1);
fprintf('Number of I2 pipi file:\t%d\n',num_I2);
I2_result = [];
guess = [0.5, I2(1,1), I2(t_size/2,1)];
for i = 1:num_I2
	I2_result = [I2_result; abs(fit(jackknifed_I2{i,1},jackknifed_I2{i,2},@(x1,x2,x3) cosh_with_const(x1,x2,x3),guess,t_size-2*sep))];
end
%%3.3.2 fit I0V
num_I0V = size(jackknifed_I0V,1);
fprintf('Number of I0V pipi file:\t%d\n',num_I0V);
I0V_result = [];
guess = [0.5, I0V(1,1), I0V(t_size/2,1)];
for i = 1:num_I2
	I0V_result = [I0V_result; abs(fit(jackknifed_I0V{i,1},jackknifed_I0V{i,2},@(x1,x2,x3) cosh_with_const(x1,x2,x3),guess,t_size-2*sep))];
end
%%3.3.3 fit I0
num_I0 = size(jackknifed_I0,1);
fprintf('Number of I0 pipi file:\t%d\n',num_I0);
I0_result = [];
guess = [0.5, I0(1,1), I0(t_size/2,1)];
for i = 1:num_I0
	I0_result = [I0_result; abs(fit(jackknifed_I0_subtracted{i,1},jackknifed_I0_subtracted{i,2},@(x1,x2,x3) cosh_with_const(x1,x2,x3),guess,t_size-2*sep))];
end

if cal_effMass
	eff_mass(jackknifed_I2_full,  t_size - 2 * sep, mean(I2_result(:,3)),  num_I2,  'I2_Meff');
	eff_mass(jackknifed_I0V_full, t_size - 2 * sep, mean(I0V_result(:,3)), num_I0V, 'I0V_Meff');
	eff_mass(jackknifed_I0_full,  t_size - 2 * sep, mean(I0_result(:,3)),  num_I0,  'I0_Meff');
end

if do_kaon
	%4 kaon calculation
	%4.1 import data
	kaon_corr = tavg_import_corr_bin(place,kaon_name,all_traj,t_trans,field,t_size,bin_size);
	%4.2 jackknife data
	jackknifed_kaoncorr = jackknife(kaon_corr,kaon_fit_range,frozen,correlated);
	jackknifed_kaoncorr_full = jackknife(kaon_corr,1:t_size,frozen,correlated);
	%4.3 fit
	num_kaoncorr = size(jackknifed_kaoncorr,1);
	fprintf('Number of kaon_corr file:\t%d\n',num_kaoncorr);
	kaon_result = [];
	guess = [0.3, 1e+6];
	for i = 1:num_kaoncorr
		kaon_result = [kaon_result; abs(fit(jackknifed_kaoncorr{i,1},jackknifed_kaoncorr{i,2},@(x1,x2,x3) cosh_no_const(x1,x2,x3),guess,t_size))];
	end
	%4.4 effective mass
	if cal_effMass
		eff_mass(jackknifed_kaoncorr_full,  t_size, 0.,  num_kaoncorr,  'kaon_Meff');
	end
end

%5 display result
fileID = fopen('results','a');
fprintf(fileID,'%s\n',datestr(clock));
fprintf(fileID,'Data file:\t%s\n',place);

if(frozen)
	fprintf('frozen ');
	fprintf(fileID,'frozen ');
else 
	fprintf('unfrozen ');
	fprintf(fileID, 'unfrozen ');
end
if(correlated)
	fprintf('correlated\n');
	fprintf(fileID, 'correlated\n');
else 
	fprintf('uncorrelated\n');
	fprintf(fileID, 'uncorrelated\n');
end

fprintf(fileID,'Number of pion_corr file:\t%d\n',num_pioncorr);
fprintf(fileID,'Number of I2 pipi file:\t%d\n',num_I2);
fprintf(fileID,'Number of I0V pipi file:\t%d\n',num_I0V);
fprintf(fileID,'Number of I0 pipi file:\t%d\n',num_I0);
fprintf(fileID,'Pi Pi separation is:\t%d\n',sep);

fprintf('pion mass:\t%.6f\t std:\t %.6f\n',mean(pion_result(:,1)),std(pion_result(:,1)) * sqrt(num_pioncorr-1));
if do_kaon 
	fprintf('kaon mass:\t%.6f\t std:\t %.6f\n',mean(kaon_result(:,1)),std(kaon_result(:,1)) * sqrt(num_kaoncorr-1));
end
fprintf('pipi_I2 energy:\t%.6f\t std:\t %.6f\n',mean(I2_result(:,1)),std(I2_result(:,1)) * sqrt(num_I2-1));
fprintf('pipi_I0V energy:\t%.6f\t std:\t %.6f\n',mean(I0V_result(:,1)),std(I0V_result(:,1)) * sqrt(num_I0V-1));
fprintf('pipi_I0 energy:\t%.6f\t std:\t %.6f\n',mean(I0_result(:,1)),std(I0_result(:,1)) * sqrt(num_I0-1));

%fprintf(fileID,'pion mass:\t%.6f\t std:\t %.6f\t fit range:\t[%d:%d]\n',mean(pion_result(:,1)),std(pion_result(:,1)) * sqrt(num_pioncorr-1), pion_fit_range(1),pion_fit_range(size(pion_fit_range,2)));
%fprintf(fileID,'pipi_I2 energy:\t%.6f\t std:\t %.6f\t fit range:\t[%d:%d]\n',mean(I2_result(:,1)),std(I2_result(:,1)) * sqrt(num_I2-1), I2_fit_range(1),I2_fit_range(size(I2_fit_range,2)));
%fprintf(fileID,'pipi_I0V energy:\t%.6f\t std:\t %.6f\t fit range:\t[%d:%d]\n',mean(I0V_result(:,1)),std(I0V_result(:,1)) * sqrt(num_I0V-1), I0V_fit_range(1),I0V_fit_range(size(I0V_fit_range,2)));
%fprintf(fileID,'pipi_I0 energy:\t%.6f\t std:\t %.6f\t fit range:\t[%d:%d]\n',mean(I0_result(:,1)),std(I0_result(:,1)) * sqrt(num_I0-1), I0_fit_range(1),I0_fit_range(size(I0_fit_range,2)));
%fprintf(fileID,'\n\n\n');
