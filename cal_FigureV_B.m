function FigureV_B = cal_FigureV_B(FigureVdis1, FigureVdis2, sep)
	if isequaln(size(FigureVdis1),size(FigureVdis2)) == 0
		fprintf('cal_FigureV_B::size of two disconnected piece dont agree.\n');
	end

	FigureV_B = zeros(size(FigureVdis1));
	t_size = size(FigureVdis1,1);
	count = size(FigureVdis1,2);

	sum1 = sum(FigureVdis1,2);
	sum2 = sum(FigureVdis2,2);
	jackVdis1 = (repmat(sum1,1,count) - FigureVdis1) / (count-1);
	jackVdis2 = (repmat(sum2,1,count) - FigureVdis2) / (count-1);
	

	for t = 1:t_size
		offset = rem(t-1+sep,t_size);
		FigureV_B(t,:) = mean(jackVdis1 .* rotate(jackVdis2, offset));
	end

	clearvars -except FigureV_B;
end
