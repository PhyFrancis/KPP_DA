function y = map_C(C_id, M_left_id, M_right_id, VA) 
	% map the index to type[1-4] offset
	
	VA_id = 0;
	if VA == 'AV'
		VA_id = 1;
	elseif VA == 'VA'
		VA_id = 0;
	else
		fprintf('Warning: referring neither VA nor AV contraction!\n');
	end

	if C_id > 22 
		C_id = C_id - 22;
	elseif C_id > 12
		C_id = C_id - 12;
	elseif C_id > 6
		C_id = C_id - 6;
	end

  y = (C_id - 1) * 8 + (M_left_id * 2 + M_right_id) * 2 + VA_id + 1;
end
