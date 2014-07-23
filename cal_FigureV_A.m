function FigureV_A = cal_FigureV_A(FigureVdis1,FigureVdis2,sep)
	if isequaln(size(FigureVdis1),size(FigureVdis2)) == 0
		fprintf('cal_FigureV_A::size of two disconnected piece dont agree.\n');
	end

	t_size = size(FigureVdis1,1);
	FigureV_A = zeros(size(FigureVdis1));

	for t = 1:t_size
		offset = rem(t-1+sep,t_size);
		FigureV_A(t,:) = mean(FigureVdis1 .* rotate(FigureVdis2, offset));
	end

	clearvars -except FigureV_A;
end
