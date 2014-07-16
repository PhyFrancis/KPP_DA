function FigureV_A = cal_FigureV_A(FigureVdis1,FigureVdis2,sep)
	if isequaln(size(FigureVdis1),size(FigureVdis2)) == 0
		fprintf('cal_FigureV_A::size of two disconnected piece dont agree.\n');
	end
	t_size = size(FigureVdis1,1);
	count = size(FigureVdis1,2);
	FigureV_A = [];
	for i = 1:count
		tmp_col = [];
		for dis = 0:t_size-1
			tmp = 0;
			for src = 1:t_size
				snk = src + dis + sep;
			  while snk > t_size
					snk = snk - t_size;
				end
				tmp = tmp + FigureVdis1(src,i) * FigureVdis2(snk,i);
			end
			tmp_col = [tmp_col; tmp / t_size];
		end
		FigureV_A = [FigureV_A, tmp_col];
	end

	clear t_size;
	clear count;
	clear tmp_col;
	clear tmp;
end
