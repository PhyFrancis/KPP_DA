function subtracted = vac_sub(jackknifed_data, sub_term)
  % This is to subtract the vacuum
	count = size(jackknifed_data,1);
	if size(sub_term,2) ~= count
		fprintf('The size of vacuum subtraction doesnt match!\n');
	end

	subtracted = cell(count,2);
	for i = 1:count
		subtracted{i,1} = jackknifed_data{i,1};
		subtracted{i,2} = jackknifed_data{i,2};
		t = jackknifed_data{i,1}(:,1) + 1;
		t_off = 1:size(jackknifed_data{i,1},1);
		subtracted{i,1}(t_off,2) = subtracted{i,1}(t_off,2) - sub_term(t,i);
	end

	clear count;
	clear t;
	clear t_off;
	clear i;
end
