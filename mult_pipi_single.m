function KQiPiPi = mult_pipi_single(KQi, PiPi, t_size, t_trans, deltat)
	% almost same as mult_pipi except KQi contains 1 operator instead of 10
	% KQi is < K(t) Q_i(t') >
	% PiPi is < PiPi(t) >
	% KQiPiPi is < K(t) Q_i(t') PiPi(t+deltat) >
	
	conf = size(PiPi, 2);
	nTslice = 0;
	KQiPiPi = zeros(t_size,conf);
	for tk = 0:t_size-1
		if ~any(tk == t_trans)
			continue;
		end
		pipi_tmp = rotate(PiPi, rem(tk+deltat,t_size));
		kqi_tmp = KQi(tk*t_size+1:(tk+1)*t_size,:);
		KQiPiPi = KQiPiPi + (kqi_tmp .* pipi_tmp);
		nTslice = nTslice + 1;
	end
	KQiPiPi = KQiPiPi / nTslice;
end
