function data_bin = bin(data,bin_size)
	% the 1:bin_size columes in data are bin-ed into one colume,
	% and the same for bin_size+1:bin_size*2 columes, and so on
	
	count_tot = size(data,2);
	count = count_tot / bin_size;
	data_bin = [];

	for i = 1:count
		data_bin = [data_bin, mean(data(:,(1+(i-1)*bin_size):(i*bin_size)),2)];
	end

	clear count_tot;
	clear count;

end
