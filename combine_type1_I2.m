function y = combine_type1(type1_data) 
	% y is a cell(10,1)
	% type1 consistants of C1 to C6
	
	y = cell(10,1);

	nConf = size(type1_data,2) / 48;

	y{1,1} = 1.0 / (3^0.5) * (...
		(+1) * type1_data(:,map_C(1,1,0,'VA')+[0:nConf-1]*48) +...
		(-1) * type1_data(:,map_C(1,1,0,'AV')+[0:nConf-1]*48) +...
		(+1) * type1_data(:,map_C(4,0,1,'VA')+[0:nConf-1]*48) +...
		(-1) * type1_data(:,map_C(4,0,1,'AV')+[0:nConf-1]*48) +...
		(+1) * type1_data(:,map_C(4,0,0,'VA')+[0:nConf-1]*48) +...
		(+1) * type1_data(:,map_C(4,0,0,'AV')+[0:nConf-1]*48) ); 

	y{7,1} = (3^0.5) / 2.0 * (...
		(+1) * type1_data(:,map_C(1,1,0,'VA')+[0:nConf-1]*48) +...
		(+1) * type1_data(:,map_C(1,1,0,'AV')+[0:nConf-1]*48) +...
		(-1) * type1_data(:,map_C(4,0,1,'VA')+[0:nConf-1]*48) +...
		(-1) * type1_data(:,map_C(4,0,1,'AV')+[0:nConf-1]*48) +...
		(-1) * type1_data(:,map_C(4,0,0,'VA')+[0:nConf-1]*48) +...
		(+1) * type1_data(:,map_C(4,0,0,'AV')+[0:nConf-1]*48) );

	y{8,1} = (3^0.5) / 2.0 * (...
		(-1) * type1_data(:,map_C(3,0,0,'VA')+[0:nConf-1]*48) +...
		(+1) * type1_data(:,map_C(3,0,0,'AV')+[0:nConf-1]*48) +...
		(-1) * type1_data(:,map_C(6,0,0,'VA')+[0:nConf-1]*48) +...
		(+1) * type1_data(:,map_C(6,0,0,'AV')+[0:nConf-1]*48) +...
		(-1) * type1_data(:,map_C(5,0,1,'VA')+[0:nConf-1]*48) +...
		(-1) * type1_data(:,map_C(5,0,1,'AV')+[0:nConf-1]*48) );

end
