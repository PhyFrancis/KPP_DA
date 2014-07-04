%This calculaton method is based on Physical review D 70 074513

function phi = phi_q_twist(q, dx, dy, dz) %answe is expressed in degree
%parameters needed {q, dx, dy, dz}
d=[dx, dy, dz];
%here d stand for whether there is a twist or not in the x,y,z direction
%If there is a twist, use 1 or else use 0
%relative momentum in unite of 2*pi/L 

N=5; %number of summation and intrgations

result=0.0;
q2=q.^2;
if norm(d)==0
	dnorm=[0,0,0];
else
	dnorm=d/norm(d);
end

for nx=-N:N
	Ny=floor(sqrt(N*N-nx*nx));
for ny=-Ny:Ny
	Nz=round(sqrt(N*N-nx*nx-ny*ny));
for nz=-Nz:Nz
n=[nx,ny,nz];

%summation part (1)
r2 = (nx+dx/2.0)^2+(ny+dy/2.0)^2+(nz+dz/2.0)^2;
result = result + exp(q2-r2)./(r2-q2);
if nx==0 && ny==0 &&nz==0
	continue
end
% integration part (2)
np=n-(n*dnorm')*d; % n pependicular
gn2=(n*dnorm')^2+np*np';
int_n = quad(@(t)zeta_func(t, q, gn2), 1e-6, 1);
result = result + int_n*(-1)^(n*d');
end
end
end

%constant part (do not depend on n) (3)
const = 0.0;
no_warn = false;
for l=0:9
	c_aux = q2.^l/factorial(l)/(l-0.5);
	const = const + c_aux;
end
for l=10:100
	c_aux = q2.^l/factorial(l)/(l-0.5);
	if c_aux < 1e-8 * abs(const)
		no_warn = true;
		break;
	end
	const = const + c_aux;
end

if ~no_warn
	warning(' reaches the maximum loop number when doing the summation');
end

result = result + pi^1.5*const;

zeta_00 = result/sqrt(4*pi);
pi^1.5/zeta_00;
phi = atan(-q*pi^1.5/zeta_00);
