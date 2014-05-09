%data analysis for meson mass
%
%1 initialize fitting parameters:
place = '/home/daiqian/BGQ/32nt_kpp/L500_11_S0_11_Mu_0.000_Ms_0.045/';
frozen = 0;% 1 if doing frozen jackknife
correlated = 0; % 1 if doing correlated fitting
cal_effMass = 0;
all_traj = [328:8:336,342:8:390];
%all_traj = [320:8:328];
t_size = 64;
t_trans = 0:63;
deltat = 18;
sep = 4;
Type1_name = ['type1','_deltat_',int2str(deltat),'_sep_',int2str(sep)];
Type2_name = ['type2','_deltat_',int2str(deltat),'_sep_',int2str(sep)];
Type3_name = ['type3','_deltat_',int2str(deltat),'_sep_',int2str(sep)];
FigureVdis_name = ['FigureVdis_sep',int2str(sep)];
bin_size = 1;
field = 3; % 3rd colume of corr file is the real part of correlator

%2 I0 kpipi calculation

%2.1 import data
FigureVdis = import_corr_bin(place,FigureVdis_name,all_traj,2,bin_size);
Type1_data = tavg_import_corr_bin(place,Type1_name,all_traj,t_trans,2*[1:48]+1,t_size,bin_size);
Type2_data = tavg_import_corr_bin(place,Type2_name,all_traj,t_trans,2*[1:48]+1,t_size,bin_size) * 0.5; % 0.5 because mom permutation
Type3_data = tavg_import_corr_bin(place,Type3_name,all_traj,t_trans,2*[1:80]+1,t_size,bin_size) * 0.5; % 0.5 because mom permutation

fprintf('Number of type1 file:\t%d\n',size(Type1_data,2)/48);
fprintf('Number of type2 file:\t%d\n',size(Type2_data,2)/48);
fprintf('Number of type3 file:\t%d\n',size(Type3_data,2)/80);

%2.2 combine data
Qi_type1 = combine_type1(Type1_data);
Qi_type2 = combine_type2(Type2_data);
Qi_type3 = combine_type3(Type3_data);

avg_type1 = [];
avg_type12 = [];
for i = 1:10
	avg_type1 = [avg_type1, mean(Qi_type1{i,1},2), std(Qi_type1{i,1},1,2)];
	avg_type12 = [avg_type12, mean(Qi_type1{i,1}+Qi_type2{i,1},2), std(Qi_type1{i,1}+Qi_type2{i,1},1,2)];
end

fn_avg_type1  = ['type1_deltat',num2str(deltat)];
fn_avg_type12 = ['type12_deltat',num2str(deltat)];
csvwrite(fn_avg_type1 ,[[0:63]',avg_type1]);
csvwrite(fn_avg_type12,[[0:63]',avg_type12]);
system(['./change_csv.sh ',fn_avg_type1] ,'-echo');
system(['./change_csv.sh ',fn_avg_type12] ,'-echo');


%Q1 = - (1/3)^0.5*(-Type1{1,1}+Type1{2,1}+Type1{5,1}+Type1{6,1}+Type1{7,1}-Type1{8,1});
%Q7 = - (0.5*3^0.5)*(-Type1{1,1}-Type1{2,1}+Type1{5,1}-Type1{6,1}+Type1{7,1}+Type1{8,1});
%Q8 = - (0.5*3^0.5)*(-Type1{3,1}+Type1{4,1}-Type1{9,1}+Type1{10,1}+Type1{11,1}+Type1{12,1});
%zk = (kaon_result(:,2)/32/2).^0.5;
%zpipi = I2_result(:,2).^0.5;
%%Q1_avg = mean(Q1,2);
%%Q7_avg = mean(Q7,2);
%%Q8_avg = mean(Q8,2);
%%2.2 jackknife data
%jackknifed_Q1 = jackknife(Q1,I2_kpp_fit_range,frozen,correlated);
%jackknifed_Q7 = jackknife(Q7,I2_kpp_fit_range,frozen,correlated);
%jackknifed_Q8 = jackknife(Q8,I2_kpp_fit_range,frozen,correlated);
%for i = 1:num_I2_kpp
%	jackknifed_Q1{i,1}(:,2) = jackknifed_Q1{i,1}(:,2)/(zk(i)*zpipi(i))*2*2^0.5;
%	jackknifed_Q7{i,1}(:,2) = jackknifed_Q7{i,1}(:,2)/(zk(i)*zpipi(i))*2*2^0.5;
%	jackknifed_Q8{i,1}(:,2) = jackknifed_Q8{i,1}(:,2)/(zk(i)*zpipi(i))*2*2^0.5;
%end
%%2.3 fit
%Q1_result = [];
%Q7_result = [];
%Q8_result = [];
%guess = [0.3];
%for i = 1:num_I2_kpp
%	Q1_result = [Q1_result; fit(jackknifed_Q1{i,1},jackknifed_Q1{i,2},@(x1,x2,x3) exp_with_no_const(x1,x2,x3),guess,I2_result(i,1)-kaon_result(i,1))];
%	Q7_result = [Q7_result; fit(jackknifed_Q7{i,1},jackknifed_Q7{i,2},@(x1,x2,x3) exp_with_no_const(x1,x2,x3),guess,I2_result(i,1)-kaon_result(i,1))];
%	Q8_result = [Q8_result; fit(jackknifed_Q8{i,1},jackknifed_Q8{i,2},@(x1,x2,x3) exp_with_no_const(x1,x2,x3),guess,I2_result(i,1)-kaon_result(i,1))];
%end
%
%%3 display result
%fileID = fopen('result','a');
%fprintf(fileID,'%s\n',datestr(clock));
%fprintf(fileID,'Data file:\t%s\n',place);
%
%if(frozen)
%	fprintf('frozen ');
%	fprintf(fileID,'frozen ');
%else 
%	fprintf('unfrozen ');
%	fprintf(fileID, 'unfrozen ');
%end
%if(correlated)
%	fprintf('correlated\n');
%	fprintf(fileID, 'correlated\n');
%else 
%	fprintf('uncorrelated\n');
%	fprintf(fileID, 'uncorrelated\n');
%end
%
%fprintf(fileID,'Number of pion_corr file:\t%d\n',num_pioncorr);
%fprintf(fileID,'Number of I2 pipi file:\t%d\n',num_I2);
%fprintf(fileID,'Pi Pi separation is:\t%d\n',sep);
%
%fprintf('pion mass:\t%.6e\t std:\t %.6e\n',mean(pion_result(:,1)),std(pion_result(:,1)) * sqrt(num_pioncorr-1));
%fprintf('kaon mass:\t%.6e\t std:\t %.6e\n',mean(kaon_result(:,1)),std(kaon_result(:,1)) * sqrt(num_kaoncorr-1));
%fprintf('pipi_I2 energy:\t%.6e\t std:\t %.6e\n',mean(I2_result(:,1)),std(I2_result(:,1)) * sqrt(num_I2-1));
%fprintf('Q1(I=2):\t%.6e\t std:\t %.6e\n',mean(Q1_result(:,1)),std(Q1_result(:,1)) * sqrt(num_I2_kpp-1));
%fprintf('Q7(I=2):\t%.6e\t std:\t %.6e\n',mean(Q7_result(:,1)),std(Q7_result(:,1)) * sqrt(num_I2_kpp-1));
%fprintf('Q8(I=2):\t%.6e\t std:\t %.6e\n',mean(Q8_result(:,1)),std(Q8_result(:,1)) * sqrt(num_I2_kpp-1));
