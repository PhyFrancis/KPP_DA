function data_bin = bin(data,bin_size,block_size)
	% the 1:bin_size blocks in data are bin-ed into one block,
	% and the same for bin_size+1:bin_size*2 blocks, and so on
	
	count_tot = size(data,2);
	count = count_tot / (block_size * bin_size); % number of blocks finally
	data_bin = [];

	for i = 1:count
		block_tmp = zeros(size(data,1),block_size);
		for b = 1:bin_size
			block_tmp = block_tmp + data(:,((i-1)*bin_size+b-1)*block_size+[1:block_size]);
		end
		block_tmp = block_tmp / bin_size;
		data_bin = [data_bin, block_tmp];
	end

	clearvars -except data_bin;
end
