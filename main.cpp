<<<<<<< Updated upstream
//Authors : Hanlin Lin(linhl24@mails.tsinghua.edu.cn) , Qian Huang(qian.huang@mathematik.uni-stuttgart.de; hqqh91@qq.com)
//the function of the code is to perform numerical simulation of particle clusters settling in the rare flow on Mars. The output of the code is a txt file, containing the position of the particles, the slip velocity, and the flow field velocity.
=======
/*
 *  main.cpp
 *  
 *  
 *  All rights are retained by the authors and Tsinghua University and University of Stuttgart.
 *  Please contact linhl24@mails.tsinghua.edu.cn or qian.huang@mathematik.uni-stuttgart.de for licensing inquiries.
 *  
 *  Authors: Hanlin Lin and Qian Huang
 *  Contact: linhl24@mails.tsinghua.edu.cn; qian.huang@mathematik.uni-stuttgart.de; hqqh91@qq.com
 */
>>>>>>> Stashed changes

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <ctime>
#include "function_DustFall.h"

using namespace std; 

const double PI = acos(-1.0); 
<<<<<<< Updated upstream
const int N = 512; 
double TIME_TOTAL = 300;
const double rp = 0.0003202; 
const double R0 = 0.01008;
=======
const int N = 256; 
double TIME_TOTAL = 600;
const double rp = 0.0006404; 
const double R0 = 0.008;
>>>>>>> Stashed changes
const double Re = 5.4961e-5; 
int Nq = 10000;
double q = 1.6e-19;
double m = 1.40816e-12;
double g = 3.72;
double ke = 0;
double E = 0;




double* x, * y, * z;			
double* vpx, * vpy, * vpz;		
double* usx, * usy, * usz;		
double* ufx, * ufy, * ufz;		
double* smx, * smy, * smz;		
double* Re_now;

double* dismin;					
int* mindex;					

double x_mean, y_mean, z_mean; 
double vx_mean;

double* create_vec_d(int length)
{
	double* m;

	m = (double*)malloc((unsigned)length * sizeof(double));
	if (!m) cout << "Allocation failure in create_vec_d!" << endl;

	return m;
}

double rand_real()
{
	double rd;
	rd = (double)rand() / (double)RAND_MAX;

	return rd;
}

double func_mars(double r, double kk) 
{
	double a;

	
	a = kk / r / r;
	return a;
}

double func_r(double r) 
{
	double b;
	double rnd = r / rp;
	 

	if(1 <= rnd && rnd < 1.01)
	{
		 b = 0 + (0.0086 - 0)/(1.01 - 1) * (rnd - 1);
	}
	else if(1.01 <= rnd && rnd < 1.1)
	{
		b = 0.0086 + (0.0973 - 0.0086)/(1.1 - 1.01) * (rnd - 1.01);
	}
	else if(1.1 <= rnd && rnd < 1.2)
	{
		b = 0.0973 + (0.1914 - 0.0973)/(1.2 - 1.1) * (rnd - 1.1);
	}
	else if(1.2 <= rnd && rnd < 1.3)
	{
		b = 0.1914 + (0.2734 - 0.1914)/(1.3 - 1.2) * (rnd - 1.2);
	}
	else if(1.3 <= rnd && rnd < 1.4)
	{
		b = 0.2734 + (0.3435 - 0.2734)/(1.4 - 1.3) * (rnd - 1.3);
	}
	else if(1.4 <= rnd && rnd < 1.5)
	{
		b = 0.3435 + (0.4034 - 0.3435)/(1.5 - 1.4) * (rnd - 1.4);
	}
	else if(1.5 <= rnd && rnd < 1.7)
	{
		b = 0.4034 + (0.4985 - 0.4034)/(1.7 - 1.5) * (rnd - 1.5);
	}
	else if(1.7 <= rnd && rnd < 2)
	{
		b = 0.4985 + (0.5988 - 0.4985)/(2 - 1.7) * (rnd - 1.7);
	}
	else if(2 <= rnd && rnd < 2.5)
	{
		b = 0.5988 + (0.7021 - 0.5988)/(2.5 - 2) * (rnd - 2);
	}
	else if(2.5 <= rnd && rnd < 3)
	{
		b = 0.7021 + (0.7641 - 0.7021)/(3 - 2.5) * (rnd - 2.5);
	}
	else if(3 <= rnd && rnd < 4)
	{
		b = 0.7641 + (0.8339 - 0.7641)/(4 - 3) * (rnd - 3);
	}
	else if(4 <= rnd && rnd < 5)
	{
		b = 0.8339 + (0.8717 - 0.8339)/(5 - 4) * (rnd - 4);
	}
	else if(5 <= rnd && rnd < 7)
	{
		b = 0.8717 + (0.9115 - 0.8717)/(7 - 5) * (rnd - 5);
	}
	else if(7 <= rnd && rnd < 10)
	{
		b = 0.9115 + (0.9393 - 0.9115)/(10 - 7) * (rnd - 7);
	}
    else if(10 <= rnd && rnd < 15)
	{
		b = 0.9393 + (0.9600 - 0.9393)/(15 - 10) * (rnd - 10);
	}
    else if(15 <= rnd && rnd < 20)
	{
		b = 0.9600 + (0.9700 - 0.9600)/(20 - 15) * (rnd - 15);
	}
    else if(20 <= rnd && rnd < 25)
	{
		b = 0.9700 + (0.9760 - 0.9700)/(25 - 20) * (rnd - 20);
	}
	else 
	{
		b = 1 - 1/(1.5443 * rnd); 
	}

	return b;
}

double func_theta(double r) 
{
	double c;
	double rnd = r / rp;
	
	if(1 <= rnd && rnd < 1.01)
	{
		 c = 0.3805 + (0.4658 - 0.3805)/(1.01 - 1) * rnd;
	}
	else if(1.01 <= rnd && rnd < 1.1)
	{
		c = 0.4658 + (0.6425 - 0.4658)/(1.1 - 1.01) * (rnd - 1.01);
	}
	else if(1.1 <= rnd && rnd < 1.2)
	{
		c = 0.6425 + (0.7215 - 0.6425)/(1.2 - 1.1) * (rnd - 1.1);
	}
	else if(1.2 <= rnd && rnd < 1.3)
	{
		c = 0.7215 + (0.7668 - 0.7215)/(1.3 - 1.2) * (rnd - 1.2);
	}
	else if(1.3 <= rnd && rnd < 1.4)
	{
		c = 0.7668 + (0.7965 - 0.7668)/(1.4 - 1.3) * (rnd - 1.3);
	}
	else if(1.4 <= rnd && rnd < 1.5)
	{
		c = 0.7965 + (0.8177 - 0.7965)/(1.5 - 1.4) * (rnd - 1.4);
	}
	else if(1.5 <= rnd && rnd < 1.7)
	{
		c = 0.8177 + (0.8460 - 0.8177)/(1.7 - 1.5) * (rnd - 1.5);
	}
	else if(1.7 <= rnd && rnd < 2)
	{
		c = 0.8460 + (0.8717 - 0.8460)/(2 - 1.7) * (rnd - 1.7);
	}
	else if(2 <= rnd && rnd < 2.5)
	{
		c = 0.8717 + (0.8967 - 0.8717)/(2.5 - 2) * (rnd - 2);
	}
	else if(2.5 <= rnd && rnd < 3)
	{
		c = 0.8967 + (0.9124 - 0.8967)/(3 - 2.5) * (rnd - 2.5);
	}
	else if(3 <= rnd && rnd < 4)
	{
		c = 0.9124 + (0.9321 - 0.9124)/(4 - 3) * (rnd - 3);
	}
	else if(4 <= rnd && rnd < 5)
	{
		c = 0.9321 + (0.9444 - 0.9321)/(5 - 4) * (rnd - 4);
	}
	else if(5 <= rnd && rnd < 7)
	{
		c = 0.9444 + (0.9592 - 0.9444)/(7 - 5) * (rnd - 5);
	}
	else if(7 <= rnd && rnd < 10)
	{
		c = 0.9592 + (0.9708 - 0.9592)/(10 - 7) * (rnd - 7);
	}
	else if(10 <= rnd && rnd < 15)
	{
		c = 0.9708 + (0.9802 - 0.9708)/(15 - 10) * (rnd - 10);
	}
    else if(15 <= rnd && rnd < 20)
	{
		c = 0.9802 + (0.9850 - 0.9802)/(20 - 15) * (rnd - 15);
	}
    else if(20 <= rnd && rnd < 25)
	{
		c = 0.9850 + (0.9879 - 0.9850)/(25 - 20) * (rnd - 20);
	}
	else 
	{
		c = 1 - 1/(2.8775 * rnd); 
	}

	return c;
}

void initialize(int seed)
{
	double r, ang1, ang2;

	srand((unsigned)time(NULL)); 
	x = create_vec_d(N); 
	y = create_vec_d(N);
	z = create_vec_d(N);

	
	for (int i = 0; i < N; i++)
	{
		r = pow(rand_real(), 1.0 / 3.0);
		ang1 = rand_real() * 2.0 * PI;
		ang2 = acos(rand_real() * 2.0 - 1.0);

		x[i] = r * cos(ang1) * sin(ang2);
		y[i] = r * sin(ang1) * sin(ang2);
		z[i] = r * cos(ang2);
	}

	

	dismin = create_vec_d(N); 
	mindex = (int*)create_vec_d(N); 

	vpx = create_vec_d(N);
	vpy = create_vec_d(N);
	vpz = create_vec_d(N);

	usx = create_vec_d(N);
	usy = create_vec_d(N);
	usz = create_vec_d(N);
	Re_now = create_vec_d(N);

	smx = create_vec_d(N);
	smy = create_vec_d(N);
	smz = create_vec_d(N);

	ufx = create_vec_d(N);
	ufy = create_vec_d(N);
	ufz = create_vec_d(N);

	for (int i = 0; i < N; i++) {
		ufx[i] = 0.0;
		ufy[i] = 0.0;
		ufz[i] = 0.0;

		vpx[i] = 0.0;
		vpy[i] = 0.0;
		vpz[i] = 0.0;

		smx[i] = 0.0;
		smy[i] = 0.0;
		smz[i] = 0.0;


	}

	return;
}

void check_overlap()
{
	int i, j;
	double dx, dy, dz;
	double r_ij;
	double dist_mv;
	bool bl_touch;

	while (true)
	{
		bl_touch = false;

		for (i = 0; i < N - 1; i++)
		{
			for (j = i + 1; j < N; j++)
			{
				dx = x[i] - x[j];
				dy = y[i] - y[j];
				dz = z[i] - z[j];
				r_ij = sqrt(dx * dx + dy * dy + dz * dz);
				if (r_ij <= 2.0 * rp)
				{
					bl_touch = true;
					
					dist_mv = 3.0 * rp / r_ij - 1;
					x[i] += 2.0 * dx * dist_mv;
					y[i] += 2.0 * dy * dist_mv;
					z[i] += 2.0 * dz * dist_mv;
				}
			}
		}
		if (!bl_touch) break;
	}

	x_mean = 0.0;
	y_mean = 0.0;
	z_mean = 0.0;
	for (i = 0; i < N; i++)
	{
		x_mean += x[i];
		y_mean += y[i];
		z_mean += z[i];
	}
	x_mean /= double(N);
	y_mean /= double(N);
	z_mean /= double(N);

	return;
}


void cal_mars()
{
	int i, j;
	double dx, dy, dz;
	double a; 
	double r_ij;

	for (i = 0; i < N; i++)
	{
		
		smx[i] = 0.0;
		smy[i] = 0.0;
		smz[i] = 0.0;
	}

	
	

	
	for (i = 0; i < N - 1; i++)
	{
		for(j = i + 1; j < N ; j++)
		{
			
				dx = x[i] - x[j];
				dy = y[i] - y[j];
				dz = z[i] - z[j];
				r_ij = sqrt(dx * dx + dy * dy + dz * dz);
				
				dx = (x[i] - x[j]) / r_ij;
				dy = (y[i] - y[j]) / r_ij;
				dz = (z[i] - z[j]) / r_ij;

				a = func_mars(r_ij , ke); 
				smx[i] = smx[i] + a * dx;
				smy[i] = smy[i] + a * dy;
				smz[i] = smz[i] + a * dz;

                smx[j] = smx[j] - a * dx;
				smy[j] = smy[j] - a * dy;
				smz[j] = smz[j] - a * dz;


			
		}
	}

	return;
}


void cal_mars_Oseen()
{
	int i, j;
	double dx, dy, dz;
	double a,b,c;
	double r_ij;
	double usx_i, usy_i, usz_i;
	double usx_j, usy_j, usz_j;
	double us_i, us_j;
	double Res_i, Res_j;
	double rnd, rnd2, rnd3;
	double cosj, sinj, efj, cosi, sini, efi;
	double ur_i, uth_i, ur_j, uth_j;
	double thx1, thy1, thz1, thx2, thy2, thz2, thdrt; 

	for (i = 0; i < N; i++)
	{
		

		ufx[i] = 0.0;
		ufy[i] = 0.0;
		ufz[i] = 0.0;
		smx[i] = 0.0;
		smy[i] = 0.0;
		smz[i] = 0.0;
	}

	for (i = 0; i < N - 1; i++)
	{
		for (j = i + 1; j < N; j++)
		{
			dx = x[i] - x[j];
			dy = y[i] - y[j];
			dz = z[i] - z[j];
			r_ij = sqrt(dx * dx + dy * dy + dz * dz);

			

			
			usx_j = usx[j];
			usy_j = usy[j];
			usz_j = usz[j];
			us_j = sqrt(usx_j * usx_j + usy_j * usy_j + usz_j * usz_j);
			Res_j = Re_now[j];

			cosj = (dx * usx_j + dy * usy_j + dz * usz_j) / (us_j * r_ij);
			sinj = sqrt(abs(1.0 - cosj * cosj));
			
			b = func_r(r_ij);
			ur_i =  cosj * (1 - b);

			c = func_theta(r_ij);
			uth_i = sinj * (c - 1);

			
			thx1 = usy_j * dz - usz_j * dy;
			thy1 = usz_j * dx - usx_j * dz;
			thz1 = usx_j * dy - usy_j * dx;

			thx2 = thy1 * dz - thz1 * dy;
			thy2 = thz1 * dx - thx1 * dz;
			thz2 = thx1 * dy - thy1 * dx;
			thdrt = sqrt(thx2 * thx2 + thy2 * thy2 + thz2 * thz2);

			
			ufx[i] += us_j * (ur_i * dx / r_ij + uth_i * thx2 / thdrt);
			ufy[i] += us_j * (ur_i * dy / r_ij + uth_i * thy2 / thdrt);
			ufz[i] += us_j * (ur_i * dz / r_ij + uth_i * thz2 / thdrt);

			
			usx_i = usx[i];
			usy_i = usy[i];
			usz_i = usz[i];
			us_i = sqrt(usx_i * usx_i + usy_i * usy_i + usz_i * usz_i);
			Res_i = Re_now[i];

			cosi = -(dx * usx_i + dy * usy_i + dz * usz_i) / (us_i * r_ij); 
			sini = sqrt(abs(1.0 - cosi * cosi));
			

			
			b = func_r(r_ij);
			ur_j =  cosi * (1 - b);

			c = func_theta(r_ij);
			uth_j = sini * (c - 1);

			
			thx1 = usy_i * (-dz) - usz_i * (-dy);
			thy1 = usz_i * (-dx) - usx_i * (-dz);
			thz1 = usx_i * (-dy) - usy_i * (-dx);

			thx2 = thy1 * (-dz) - thz1 * (-dy);
			thy2 = thz1 * (-dx) - thx1 * (-dz);
			thz2 = thx1 * (-dy) - thy1 * (-dx);
			thdrt = sqrt(thx2 * thx2 + thy2 * thy2 + thz2 * thz2);

		
			ufx[j] += us_i * (ur_j * (-dx) / r_ij + uth_j * thx2 / thdrt);
			ufy[j] += us_i * (ur_j * (-dy) / r_ij + uth_j * thy2 / thdrt);
			ufz[j] += us_i * (ur_j * (-dz) / r_ij + uth_j * thz2 / thdrt);
		}
	}


	for (i = 0; i < N - 1; i++)
	{
		for(j = i + 1;j < N; j++)
		{
				dx = x[i] - x[j];
				dy = y[i] - y[j];
				dz = z[i] - z[j];
				r_ij = sqrt(dx * dx + dy * dy + dz * dz);
				
				dx = (x[i] - x[j]) / r_ij;
				dy = (y[i] - y[j]) / r_ij;
				dz = (z[i] - z[j]) / r_ij;

				a = func_mars(r_ij, ke); 

				smx[i] = smx[i] + a * dx;
				smy[i] = smy[i] + a * dy;
				smz[i] = smz[i] + a * dz;
                smx[j] = smx[j] - a * dx;
				smy[j] = smy[j] - a * dy;
				smz[j] = smz[j] - a * dz;
			}
		}
	
	return;
}

void cal_slip_velocity()
{
	double us;

	for (int i = 0; i < N; i++)
	{
		usx[i] = smx[i]; 
		usy[i] = smy[i];
		usz[i] = smz[i] + 1;
		us = sqrt(usx[i] * usx[i] + usy[i] * usy[i] + usz[i] * usz[i]);

		Re_now[i] = us * Re;
	}
	return;
}

double update_velocity() 
{
	double vp_temp;
	double vpmax = 0.0;
	double vtot_x = 0.0;

	for (int i = 0; i < N; i++)
	{
		vpx[i] = ufx[i] + usx[i];
		vpy[i] = ufy[i] + usy[i];
		vpz[i] = ufz[i] + usz[i];

		vtot_x += vpx[i];
	}
	vx_mean = vtot_x / (double)N;
	for (int i = 0; i < N; i++)
	{
		vp_temp = sqrt((vpx[i] - vx_mean) * (vpx[i] - vx_mean) + vpy[i] * vpy[i] + vpz[i] * vpz[i]);
		vpmax = vp_temp > vpmax ? vp_temp : vpmax;
	}
	return vpmax;
}

void update_position(double dt, double* vpx_old, double* vpy_old, double* vpz_old)
{
	for (int i = 0; i < N; i++)
	{
		x[i] += (1.5 * vpx[i] - 0.5 * vpx_old[i]) * dt;
		y[i] += (1.5 * vpy[i] - 0.5 * vpy_old[i]) * dt;
		z[i] += (1.5 * vpz[i] - 0.5 * vpz_old[i]) * dt;

		
		vpx_old[i] = vpx[i];
		vpy_old[i] = vpy[i];
		vpz_old[i] = vpz[i];
	}
	check_overlap();
};

int main()
{
	double seed;
	double vpmax;					
	double t_total = 0.0, dt = 0.0;		

	double tout = 0.0, dtout = 0.2; 
	double dfact = 0.1;				

	
	seed = 0.5976;

	ke= Nq * Nq / R0 / R0 * 8.988e9 * q * q/ (m * g + E* Nq * q);
	
	double* vpx_old = create_vec_d(N); 
	double* vpy_old = create_vec_d(N);
	double* vpz_old = create_vec_d(N);

	
	char filename[150];
	char str_N[5], str_seed[12],  str_Re[10], str_TIME_TOTAL[8],str_Nq[9],str_R0[6];
	sprintf(str_seed, "%.4f", seed);
	sprintf(str_N, "%d", N);
	sprintf(str_Re, "%.6f", Re);
	sprintf(str_Nq, "%d", Nq);
	sprintf(str_R0, "%.3f", R0);
	
	
	
	sprintf(str_TIME_TOTAL, "%.1f", TIME_TOTAL);
	

	sprintf(filename, "rare_25_N%sTIME%sseed%sRe%sNq%sR0%s.txt", str_N,str_TIME_TOTAL, str_seed, str_Re,str_Nq,str_R0);
	
	ofstream cfile(filename);
	cfile.precision(8);
	

	
	initialize(seed);
	check_overlap();
	cal_mars(); 
	cal_slip_velocity();

	
	time_t now = time(0);				
	char* time_string = ctime(&now);	
	cout << "Start calculation: " << time_string << endl;  
	
	for (int k = 1; t_total < TIME_TOTAL; t_total += dt)
	{
		cal_mars_Oseen(); 
		vpmax = update_velocity();
		cal_slip_velocity(); 

		if (k == 1)
		{
			k = 0;
			for (int i = 0; i < N; i++)
			{
				vpx_old[i] = vpx[i];
				vpy_old[i] = vpy[i];
				vpz_old[i] = vpz[i];
			}
		}

		
		if (t_total >= tout)
		{
			tout += dtout;

			cout << "Time = " << t_total << "; dt = " << dt << endl;
			
			for (int i = 0; i < N; i++)
			{
				cfile << x[i] - x_mean << "\t" << y[i] - y_mean << "\t" << z[i] - z_mean << "\t"
					<< ufx[i] << "\t" << ufy[i] << "\t" << ufz[i] << "\t"  << usx[i] << "\t" << usy[i] << "\t" << usz[i] << "\t"  << endl;
				
			}
		}

		
		dt = dfact * rp / vpmax;
		dt = dt < 0.1 ? dt : 0.1;

		update_position(dt, vpx_old, vpy_old, vpz_old); 

	}

	cfile.close();

	now = time(0);
	time_string = ctime(&now);
	cout << "Calculation completed: " << time_string << endl;

	return 0;
}
