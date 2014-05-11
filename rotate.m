function rData = rotate(data, row_offset) 
	% 'data' is original matrix
	% 'rData' is rotated one
	
	row = size(data,1);
	rData = data((row_offset+1):row,:);
	rData = [rData ; data(1:row_offset,:)];
end
