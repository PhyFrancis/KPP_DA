%data analysis for meson spectrum

%1 initialize fitting parameters:
run('conf.m');

nGBC = GparityX + GparityY + GparityZ;
frozen = 0;% 1 if doing frozen jackknife
correlated = 0; % 1 if doing correlated fitting
cal_effMass = 1;
t_trans = 0:(t_size-1);
sep = 4;
do_kaon = 1;
pion_name = 'pioncorr';
pion_name_2 = 'pioncorr_offDiag';
kaon_name = 'kaoncorr';
FigureC_name = ['FigureC','_sep',int2str(sep)];
FigureD_name = ['FigureD','_sep',int2str(sep)];
field = 3; % 3rd colume of corr file is the real part of correlator
pion_fit_range = 7:(t_size+2-7); 
kaon_fit_range = 7:(t_size+2-7); 
I2_fit_range = [8:(t_size/2-sep),(t_size/2+2-sep):(t_size+2-8-2*sep)];

%2 pion calculation
%2.1 import data
pion_corr   = tavg_import_corr_bin(place,pion_name  ,all_traj,t_trans,field,t_size,bin_size);
pion_corr_2 = tavg_import_corr_bin(place,pion_name_2,all_traj,t_trans,field,t_size,bin_size);
%2.2 jackknife data
jackknifed_pioncorr = jackknife(pion_corr,pion_fit_range,frozen,correlated);
jackknifed_pioncorr_full = jackknife(pion_corr,1:t_size,frozen,correlated);
jackknifed_pioncorr_2 = jackknife(pion_corr_2,pion_fit_range,frozen,correlated);
jackknifed_pioncorr_2_full = jackknife(pion_corr_2,1:t_size,frozen,correlated);
%2.3 fit
num_pioncorr = size(jackknifed_pioncorr,1);
if (nGBC > 1) && (isequaln(num_pioncorr,size(jackknifed_pioncorr_2,1)) == 0)
	printf('Number of pioncorr, pioncorr_offDiag files dont match!\n');
end
fprintf('Number of pion_corr file:\t%d\n',num_pioncorr);
pion_result = [];
pion_result_2 = [];
guess = [0.3, 1e+6];
for i = 1:num_pioncorr
	pion_result   = [pion_result  ; abs(fit(jackknifed_pioncorr{i,1}  ,jackknifed_pioncorr{i,2}  ,@(x1,x2,x3) cosh_no_const(x1,x2,x3),guess,t_size))];
	if nGBC > 1
		pion_result_2 = [pion_result_2; abs(fit(jackknifed_pioncorr_2{i,1},jackknifed_pioncorr_2{i,2},@(x1,x2,x3) cosh_no_const(x1,x2,x3),guess,t_size))];
	end
end
%2.4 effective mass
if cal_effMass
	eff_mass(jackknifed_pioncorr_full, t_size, 0., num_pioncorr, 'pion_Meff');
	if nGBC > 1
		eff_mass(jackknifed_pioncorr_2_full, t_size, 0., num_pioncorr, 'pion_offDiag_Meff');
	end
end

if nGBC > 1
	ZpiRatio = mean(pion_result(:,2)) / mean(pion_result_2(:,2)); % The norm difference due to different p directions
end

%3 pipi calculation
%3.1 import data
switch nGBC 
	case 1
		FigureC = 0.5 * (tavg_import_corr_bin(place,[FigureC_name,'_A'],all_traj,t_trans,field,t_size,bin_size) + ...
			tavg_import_corr_bin(place,[FigureC_name,'_B'],all_traj,t_trans,field,t_size,bin_size));
		FigureD = 0.5 * (tavg_import_corr_bin(place,[FigureD_name,'_A'],all_traj,t_trans,field,t_size,bin_size) + ...
			tavg_import_corr_bin(place,[FigureD_name,'_B'],all_traj,t_trans,field,t_size,bin_size));
	case 2
		FigureC = 0.5 * (tavg_import_corr_bin(place,[FigureC_name,'_A'],all_traj,t_trans,field,t_size,bin_size) + ...
			ZpiRatio * tavg_import_corr_bin(place,[FigureC_name,'_B'],all_traj,t_trans,field,t_size,bin_size));
		FigureD = 0.5 * (tavg_import_corr_bin(place,[FigureD_name,'_A'],all_traj,t_trans,field,t_size,bin_size) + ...
			ZpiRatio * tavg_import_corr_bin(place,[FigureD_name,'_B'],all_traj,t_trans,field,t_size,bin_size));
	case 3
		FigureC = 0.25 * tavg_import_corr_bin(place,[FigureC_name,'_A'],all_traj,t_trans,field,t_size,bin_size) + ...
			0.75 * ZpiRatio * tavg_import_corr_bin(place,[FigureC_name,'_B'],all_traj,t_trans,field,t_size,bin_size);
		FigureD = 0.25 * tavg_import_corr_bin(place,[FigureD_name,'_A'],all_traj,t_trans,field,t_size,bin_size) + ...
			0.75 * ZpiRatio * tavg_import_corr_bin(place,[FigureD_name,'_B'],all_traj,t_trans,field,t_size,bin_size);
	otherwise 
		fprintf('Number of G-parity twist: %d, not supported', nGBC);
end
I2 = 2 * (FigureD - FigureC);
jackknifed_I2 = jackknife(I2,I2_fit_range,frozen,correlated);
jackknifed_I2_full = jackknife(I2,1:t_size,frozen,correlated);
%%3.3.1 fit I2
num_I2 = size(jackknifed_I2,1);
fprintf('Number of I2 pipi file:\t%d\n',num_I2);
I2_result = [];
guess = [0.5, I2(1,1), I2(t_size/2,1)];
for i = 1:num_I2
	I2_result = [I2_result; abs(fit(jackknifed_I2{i,1},jackknifed_I2{i,2},@(x1,x2,x3) cosh_with_const(x1,x2,x3),guess,t_size-2*sep))];
end
if cal_effMass
	eff_mass(jackknifed_I2_full,  t_size - 2 * sep, mean(I2_result(:,3)),  num_I2,  'I2_Meff');
end
%%3.3.4 phase shift
q_I2 = (I2_result(:,1).^2 / 4 - mpi_rest^2).^0.5 * L_size / (2 * 3.1415926);
phase_shift_I2 = [];
for i = 1:num_I2
	phase_shift_I2 = [phase_shift_I2; -phi_q_twist(q_I2(i),GparityX,GparityY,GparityZ)];
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

%6 LL factor
FF_I2 = [];
for i = 1:num_I2
	FF_I2 = [FF_I2; 
	(phase_shift_I2(i) / q_I2(i) + d_phi_d_q_twist(q_I2(i),GparityX,GparityY,GparityZ))^0.5 * kaon_result(i,1)^0.5 * I2_result(i,1) / 2^0.5 / 3.1415926 / q_I2(i) * L_size^1.5];
end
fn = ['Qi_result/FF_I2'];
save(fn,'FF_I2');

%6 display result
fprintf('pion mass:\t%.6f\t std:\t %.6f\n',mean(pion_result(:,1)),std(pion_result(:,1)) * sqrt(num_pioncorr-1));
if do_kaon 
	fprintf('kaon mass:\t%.6f\t std:\t %.6f\n',mean(kaon_result(:,1)),std(kaon_result(:,1)) * sqrt(num_kaoncorr-1));
end
fprintf('pipi_I2 energy:\t%.6f\t std:\t %.6f\n',mean(I2_result(:,1)),std(I2_result(:,1)) * sqrt(num_I2-1));
fprintf('pipi_I2 phase shift:\t%.6f\t std:\t %.6f\n',mean(phase_shift_I2(:,1)),std(phase_shift_I2(:,1)) * sqrt(num_I2-1));
