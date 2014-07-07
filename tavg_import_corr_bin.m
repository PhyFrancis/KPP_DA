function corr = tavg_import_corr_bin(place,name,conf,t_range,field,T,bin_size)

	% import as usual
	% corr = tavg_import_corr_2(place,name,conf,t_range,field,T); % sum
	corr = tavg_import_corr(place,name,conf,t_range,field,T); % avg

	% now bin the data
	corr = bin(corr,bin_size,size(field,2));

end
