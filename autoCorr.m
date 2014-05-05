function y = autoCorr(data) 
	y = [];
	lenth = size(data,1);
	avg = mean(data);
	var = std(data,1)^2;
	data = data - avg;
	for sep = 0:1:lenth-1
		y = [y; (data(1:lenth-sep)' * data(sep+1:lenth)) / (lenth-sep) / var];
	end
end
