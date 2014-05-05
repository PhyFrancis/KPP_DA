function corr = import_corr_bin(place,name,conf,field,bin_size)

	% import as usual
  corr = import_corr(place,name,conf,field);

	% bin the data
	corr = bin(corr,bin_size);

end
