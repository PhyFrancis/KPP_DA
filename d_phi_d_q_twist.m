function dphidq = d_phi_d_q_twist(q, dx, dy, dz)
	dq = 0.002;
	phi_p = phi_q_twist(q+dq,dx,dy,dz);
	phi_m = phi_q_twist(q-dq,dx,dy,dz);
	dphidq = (phi_p-phi_m) / 2.0 / dq;
end
