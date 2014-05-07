#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>

using namespace std;

float M_eff(int, int, float);

int main(int argc, char *argv[])
{
	int T = atoi(argv[1]);
	double c_ave = atof(argv[2]);
	ifstream p(argv[3]);
	int conf = atoi(argv[4]);

	//ofstream p2("now/Meff.txt");
	int skip_int;
	double skip_double;
	double corr[T];
	double mass[T][conf];
	for(int j=0;j<conf;j++)
	{
		for(int i=0;i<T;i++)
			p>>skip_int>>corr[i];
		for(int i=0;i<T-1;i++)
			mass[i][j] = M_eff(T, i, (corr[i+1]-c_ave)/(corr[i]-c_ave));
	}
	p.close();
	//p2.close();

	double ave;
  double err;
	for(int i=0;i<T-1;i++)
	{
		ave=0;
		err=0;
		for(int j=0;j<conf;j++)
		{
			ave+=mass[i][j];
			err+=mass[i][j]*mass[i][j];
		}
		ave/=conf;
		err=sqrt((err/conf-ave*ave))*sqrt(conf);
		cout<<i<<"\t"<<ave<<"\t"<<err<<endl;
	}
}

float M_eff(int Tot_T, int T_sites, float goal)
{
	float up=2.0, down=0;
	while(up-down>0.00001)
	{
		if((cosh((up+down)/2*(Tot_T/2-T_sites-1))/cosh((up+down)/2*(Tot_T/2-T_sites))-goal)*(cosh((down)/2*(Tot_T/2-T_sites-1))/cosh((down)/2*(Tot_T/2-T_sites))-goal)<0)
			up=(up+down)/2;
		else down=(up+down)/2;
	}
	return up;
}
