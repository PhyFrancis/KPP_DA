% this file specifies the basic parameters of fitting

% place = '/home/daiqian/BGQ/32nt_kpp/L900_11_S0_11_Mu_0.000_Ms_0.045/';
% % all_traj = [286:4:600];
% all_traj = [286:4:338,346:8:510,518:4:566];
% L_size = 32;
% t_size = 64;
% GparityX = 1; % 1 if using Gparity
% GparityY = 1;
% GparityZ = 1;
% bin_size = 2;
% % place = '/home/daiqian/BGQ/32nt_kpp/L500_11_S0_11_Mu_0.000_Ms_0.045/';
% % all_traj = [328,336,342:8:398];
% % place = '/home/daiqian/BGQ/32nt_kpp/L500_L900_avg/';
% % all_traj = [342:8:398];
% mpi_rest = 0.10421; % from 32nt64 DSDR paper, for phase shift calculation

place = '/home/daiqian/BGQ/32nt_kpp/L900_11_S0_11_Mu_0.000_Ms_0.045/';
all_traj = [286:4:500];
L_size = 32;
t_size = 64;
GparityX = 1; % 1 if using Gparity
GparityY = 1;
GparityZ = 1;
mpi_rest = 0.10421; % from 32nt64 DSDR paper, for phase shift calculation
bin_size = 2;

% place = '/home/daiqian/BGQ/16nt_kpp/L100_11_S0_11_Mu_0.010_Ms_0.099_2twist/';
% all_traj = [1000:10:1560];
% L_size = 16;
% t_size = 32;
% GparityX = 1; % 1 if using Gparity
% GparityY = 1;
% GparityZ = 0;
% mpi_rest = 0.2437; 
% bin_size = 1;
