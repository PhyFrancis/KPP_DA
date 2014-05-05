function FigureV_B = cal_FigureV_B(FigureVdis, sep)

	t_size = size(FigureVdis,1);
	count = size(FigureVdis,2);
	FigureV_B = [];

	for i = 1:count

		tmp = [];
		for j = 1:count
			if j==i
				continue; % omitting i's piece
			end
			tmp = [tmp, FigureVdis(:,j)];
		end
		tmp_mean = mean(tmp,2);

		tmp_col = [];
		for dis = 0:t_size-1
			tmp = 0;
			for src = 1:t_size
				snk = src + dis + sep;
				while snk > t_size
					snk = snk - t_size;
				end
				tmp = tmp + tmp_mean(src) * tmp_mean(snk);
			end
			tmp_col = [tmp_col; tmp / t_size];
		end
		FigureV_B = [FigureV_B,tmp_col];

	end

	clear t_size
	clear count;
	clear tmp;
	clear tmp_mean;
	clear tmp_col;
end
