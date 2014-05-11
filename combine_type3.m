function y = combine_type3(type3_data) 
	% y is a cell(10,1)
	% type3 consistants of C13 to C22
	
	y = cell(10,1);

	nConf = size(type3_data,2) / 80;

	y{1,1} = 1.0 / (6^0.5) * (...
		(+3) * type3_data(:,map_C(13,0,1,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(13,0,1,'AV')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(16,1,0,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(16,1,0,'AV')+[0:nConf-1]*80) );

	y{2,1} = 1.0 / (6^0.5) * (...
		(+3) * type3_data(:,map_C(14,0,1,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(14,0,1,'AV')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(17,1,0,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(17,1,0,'AV')+[0:nConf-1]*80) );

	y{3,1} = 1.0 / (6^0.5) * (...
		(+3) * type3_data(:,map_C(13,0,1,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(13,0,1,'AV')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(13,0,0,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(13,0,0,'AV')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(16,1,0,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(16,1,0,'AV')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(16,0,0,'VA')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(16,0,0,'AV')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(19,0,0,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(19,0,0,'AV')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(21,0,0,'VA')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(21,0,0,'AV')+[0:nConf-1]*80) );

	y{4,1} = 1.0 / (6^0.5) * (...
		(+3) * type3_data(:,map_C(14,0,1,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(14,0,1,'AV')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(15,0,0,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(15,0,0,'AV')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(17,1,0,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(17,1,0,'AV')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(18,0,0,'VA')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(18,0,0,'AV')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(20,0,0,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(20,0,0,'AV')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(22,0,0,'VA')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(22,0,0,'AV')+[0:nConf-1]*80) );

	y{5,1} = 1.0 / (6^0.5) * (...
		(-3) * type3_data(:,map_C(13,0,1,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(13,0,1,'AV')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(13,0,0,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(13,0,0,'AV')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(16,1,0,'VA')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(16,1,0,'AV')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(16,0,0,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(16,0,0,'AV')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(19,0,0,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(19,0,0,'AV')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(21,0,0,'VA')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(21,0,0,'AV')+[0:nConf-1]*80) );

	y{6,1} = 1.0 / (6^0.5) * (...
		(-3) * type3_data(:,map_C(14,0,1,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(14,0,1,'AV')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(15,0,0,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(15,0,0,'AV')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(17,1,0,'VA')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(17,1,0,'AV')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(18,0,0,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(18,0,0,'AV')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(20,0,0,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(20,0,0,'AV')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(22,0,0,'VA')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(22,0,0,'AV')+[0:nConf-1]*80) );

	y{7,1} = 1.0 / (6^0.5) * (...
		(-3) * type3_data(:,map_C(13,0,1,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(13,0,1,'AV')+[0:nConf-1]*80) +...
		(-3/2) * type3_data(:,map_C(13,0,0,'VA')+[0:nConf-1]*80) +...
		(+3/2) * type3_data(:,map_C(13,0,0,'AV')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(16,1,0,'VA')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(16,1,0,'AV')+[0:nConf-1]*80) +...
		(-3/2) * type3_data(:,map_C(16,0,0,'VA')+[0:nConf-1]*80) +...
		(+3/2) * type3_data(:,map_C(16,0,0,'AV')+[0:nConf-1]*80) +...
		(-3/2) * type3_data(:,map_C(19,0,0,'VA')+[0:nConf-1]*80) +...
		(+3/2) * type3_data(:,map_C(19,0,0,'AV')+[0:nConf-1]*80) +...
		(+3/2) * type3_data(:,map_C(21,0,0,'VA')+[0:nConf-1]*80) +...
		(-3/2) * type3_data(:,map_C(21,0,0,'AV')+[0:nConf-1]*80) );

	y{8,1} = 1.0 / (6^0.5) * (...
		(-3) * type3_data(:,map_C(14,0,1,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(14,0,1,'AV')+[0:nConf-1]*80) +...
		(-3/2) * type3_data(:,map_C(15,0,0,'VA')+[0:nConf-1]*80) +...
		(+3/2) * type3_data(:,map_C(15,0,0,'AV')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(17,1,0,'VA')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(17,1,0,'AV')+[0:nConf-1]*80) +...
		(-3/2) * type3_data(:,map_C(18,0,0,'VA')+[0:nConf-1]*80) +...
		(+3/2) * type3_data(:,map_C(18,0,0,'AV')+[0:nConf-1]*80) +...
		(-3/2) * type3_data(:,map_C(20,0,0,'VA')+[0:nConf-1]*80) +...
		(+3/2) * type3_data(:,map_C(20,0,0,'AV')+[0:nConf-1]*80) +...
		(+3/2) * type3_data(:,map_C(22,0,0,'VA')+[0:nConf-1]*80) +...
		(-3/2) * type3_data(:,map_C(22,0,0,'AV')+[0:nConf-1]*80) );

	y{9,1} = 1.0 / (6^0.5) * (...
		(+3) * type3_data(:,map_C(13,0,1,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(13,0,1,'AV')+[0:nConf-1]*80) +...
		(+3/2) * type3_data(:,map_C(13,0,0,'VA')+[0:nConf-1]*80) +...
		(+3/2) * type3_data(:,map_C(13,0,0,'AV')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(16,1,0,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(16,1,0,'AV')+[0:nConf-1]*80) +...
		(-3/2) * type3_data(:,map_C(16,0,0,'VA')+[0:nConf-1]*80) +...
		(-3/2) * type3_data(:,map_C(16,0,0,'AV')+[0:nConf-1]*80) +...
		(+3/2) * type3_data(:,map_C(19,0,0,'VA')+[0:nConf-1]*80) +...
		(+3/2) * type3_data(:,map_C(19,0,0,'AV')+[0:nConf-1]*80) +...
		(-3/2) * type3_data(:,map_C(21,0,0,'VA')+[0:nConf-1]*80) +...
		(-3/2) * type3_data(:,map_C(21,0,0,'AV')+[0:nConf-1]*80) );

	y{10,1} = 1.0 / (6^0.5) * (...
		(+3) * type3_data(:,map_C(14,0,1,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(14,0,1,'AV')+[0:nConf-1]*80) +...
		(+3/2) * type3_data(:,map_C(15,0,0,'VA')+[0:nConf-1]*80) +...
		(+3/2) * type3_data(:,map_C(15,0,0,'AV')+[0:nConf-1]*80) +...
		(+3) * type3_data(:,map_C(17,1,0,'VA')+[0:nConf-1]*80) +...
		(-3) * type3_data(:,map_C(17,1,0,'AV')+[0:nConf-1]*80) +...
		(-3/2) * type3_data(:,map_C(18,0,0,'VA')+[0:nConf-1]*80) +...
		(-3/2) * type3_data(:,map_C(18,0,0,'AV')+[0:nConf-1]*80) +...
		(+3/2) * type3_data(:,map_C(20,0,0,'VA')+[0:nConf-1]*80) +...
		(+3/2) * type3_data(:,map_C(20,0,0,'AV')+[0:nConf-1]*80) +...
		(-3/2) * type3_data(:,map_C(22,0,0,'VA')+[0:nConf-1]*80) +...
		(-3/2) * type3_data(:,map_C(22,0,0,'AV')+[0:nConf-1]*80) );
end
