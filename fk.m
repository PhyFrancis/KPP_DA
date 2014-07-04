% this file analysis fk : <K|\bar{s}\gamma_5 d|0> 
% TODO: consider sqrt(2) due to G-parity

% initialize parameters
run('conf.m');

% run kaon analysis
run('kaon.m');

Type4_name = 'type4';
frozen = 0;
correlated = 0;
t_trans = 0:(t_size-1);
fit_range = 1:t_size;

% get <K|A|0>
KA0 = tavg_import_corr_bin(place,Type4_name,all_traj,t_trans,2*(80+[1:2])+1,t_size,bin_size); % tavg'ed < s5d K > < u5s' K >
nKA0 = size(KA0,2)/2; fprintf('Number of KA0 file:\t%d\n',nKA0);

jackknifed_KA0 = jackknife(KA0(:,2*[1:nKA0]-1),fit_range,frozen,correlated); % odd columns of KA0: only consider < s5d K >

% normalize
for i = 1:nKA0
	coshMk = zeros(t_size,1);
	for t = 0:(t_size-1)
		coshMk(t+1) = exp(-t*kaon_result(i,1)) + exp((t-t_size)*kaon_result(i,1));
	end
	jackknifed_KA0{i,1}(:,2) = (jackknifed_KA0{i,1}(:,2) / (kaon_result(i,2))^0.5)./coshMk;
end

% 'fit' fk
fk_result = zeros(nKA0,1);
for i = 1:nKA0
	fk_result(i) = mean(jackknifed_KA0{i,1}(5:(t_size+2-5),2));
end

% show result
fprintf('fk: %.6f  std: %.6f\n',mean(fk_result),std(fk_result) * (nKA0-1)^0.5);
