function eff = eff_mass(jackknifed_data, t_size, c_avg, conf, output_name) 
	fp = fopen('meff_tmp','w');
	for i = 1:conf
		fprintf(fp,'%d %.6e\n',jackknifed_data{i,1}(1:t_size,:)');
	end
	
	% use Meff.x to output into file
	comm = ['./Meff.x ',num2str(t_size),' ',num2str(c_avg),' meff_tmp ',num2str(conf),' > ', 'Meff/',output_name];
	system(comm,'-echo');
	
	delete('meff_tmp');
end
