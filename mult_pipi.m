function KQiPiPi = mult_pipi(KQi, PiPi, t_size, t_trans, deltat)
	% KQi is < K(t) Q_i(t') >
	% PiPi is < PiPi(t) >
	% KQiPiPi is < K(t) Q_i(t') PiPi(t+deltat) >
	
	conf = size(PiPi, 2);

	KQiPiPi = cell(10,1);
	for Q = 1:10
		KQiPiPi{Q,1} = zeros(t_size,conf);
		nTslice = 0;
		for tk = 0:t_size-1
			if ~any(tk == t_trans)
				continue;
			end
			pipi_tmp = rotate(PiPi, rem(tk+deltat,t_size));
			kqi_tmp = KQi{Q,1}(tk*t_size+1:(tk+1)*t_size,:);
			KQiPiPi{Q,1} = KQiPiPi{Q,1} + (kqi_tmp .* pipi_tmp);
			nTslice = nTslice + 1;
		end
		KQiPiPi{Q,1} = KQiPiPi{Q,1} / nTslice;
	end
end

