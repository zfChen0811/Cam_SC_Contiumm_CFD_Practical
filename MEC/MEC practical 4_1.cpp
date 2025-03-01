#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>
#include <fstream>

double Minbee(double r);
std::vector<double> f(std::vector<double>con, double gama, double p_inf);
std::vector<double> con2pri(std::vector<double>con, double gama, double p_inf);
std::vector<double> pri2con(std::vector<double>pri, double gama, double p_inf);
std::vector<double> Fri_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama, double p_inf);
std::vector<double> RI_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama, double p_inf);
std::vector<double> FORCE_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama, double p_inf);
double calecon(std::vector<double> con);
double calepri(std::vector<double> pri, double gama, double p_inf);
double calcs(std::vector<double> pri, double gama, double p_inf);

double Minbee(double r){
	double result;
	if(r <= 0){
		result = 0;
	}else if(r <= 1){
		result = r;
	}else{
		if(1 < (2/(1 + r))){
			result = 1;
		}else{
			result = 2/(1 + r);
		}
	}
	return result;
}

std::vector<double> f(std::vector<double>con, double gama, double p_inf){
	// calculate f
	// based on conservative value roh rohv E
	double rho = con[0];
	double vx = con[1] / con[0];
	double E = con[2];
	std::vector<double> pri = con2pri(con, gama, p_inf);
	double p = pri[2];
	std::vector<double> f;
	f.resize(3);
	f[0] = rho * vx;
	f[1] = rho * vx * vx + p;
	f[2] = (E + p) * vx;
	return f;
}

std::vector<double> con2pri(std::vector<double>con, double gama, double p_inf){
	// convert conservative (rho rhov E) to primitive (rho v p)
	std::vector<double> pri;
	pri.resize(3);
	// rho v p
	double e = calecon(con);
	pri[0] = con[0];
	pri[1] = con[1] / con[0];
	pri[2] = (gama - 1) * con[0] * e - gama * p_inf;
	return pri;
}

std::vector<double> pri2con(std::vector<double>pri, double gama, double p_inf){
	// convert primitive (rho v p) to conservative  (rho rhov E)
	std::vector<double> con;
	con.resize(3);
	double e = calepri(pri, gama, p_inf);
	// rho rhov E
	con[0] = pri[0];
	con[1] = pri[0] * pri[1];
	con[2] = pri[0] * e + 0.5 * pri[0] * (pri[1] * pri[1]);
	return con;
}

std::vector<double> Fri_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama, double p_inf){
	// calculate the flux for lax-friedrich
	std::vector<double> Fri;
	Fri.resize(3);
	std::vector<double> f0 = f(con0, gama, p_inf);
	std::vector<double> f1 = f(con1, gama, p_inf);
	for(int i = 0; i < 3; i++){
		Fri[i] = 0.5 * (dx/dt) * (con0[i] - con1[i]) + 0.5 * (f1[i] + f0[i]);
	}
	return Fri;
}

std::vector<double> RI_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama, double p_inf){
	std::vector<double> RI;
	std::vector<double> u;
	u.resize(3);
	std::vector<double> f0 = f(con0, gama, p_inf);
	std::vector<double> f1 = f(con1, gama, p_inf);
	for(int i = 0; i < 3; i++){
		u[i] = 0.5 * (con0[i] + con1[i]) - 0.5 * (dt/dx) * (f1[i] - f0[i]);
	}
	RI = f(u, gama, p_inf);
	return RI;
}

std::vector<double> FORCE_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama, double p_inf){
	std::vector<double> flux;
	flux.resize(3);
	std::vector<double> RI = RI_flux(con0, con1, dt, dx, gama, p_inf);
	std::vector<double> Fri = Fri_flux(con0, con1, dt, dx, gama, p_inf);
	for(int i = 0; i < 3 ; i++){
		flux[i] = 0.5 * (Fri[i] + RI[i]);
	}
	return flux;
}

double calecon(std::vector<double> con){
	double rho = con[0];
	double vx = con[1] / con[0];
	double E = con[2];
	double e = (E - 0.5 * rho * (vx * vx))/rho;
	return e;
}

double calepri(std::vector<double> pri, double gama, double p_inf){
	double e = (pri[2] + gama * p_inf)/((gama - 1) * pri[0]);
	return e;
}

double calcs(std::vector<double> pri, double gama, double p_inf){
	double cs = sqrt(gama *  (pri[2] + p_inf) / pri[0]);
	return cs;
}

using namespace std;
int main(int argc, char *argv[]) {
	double x0 = 0;
	double x1 = 1;
	double gama = 7.15;
	double p_inf = 3 * pow(10,8);
	double c = 0.8;
	int nCells = 100;
	double nxPoints = nCells + 1;
	double time = 1e-4;
	double t = 0;
	double dx = (x1 - x0)/(nxPoints-1);
	
	std::vector<std::vector<double>> u_pri;
	u_pri.resize(nCells + 4,vector<double>(3));

	std::vector<std::vector<double>> u_con;
	u_con.resize(nCells + 4,vector<double>(3));
	
	std::vector<double> wx;
	wx.resize(3,0);
	
	std::vector<std::vector<double>> deltax0;
	deltax0.resize(nCells + 4,vector<double>(3));
	
	std::vector<std::vector<double>> deltax1;
	deltax1.resize(nCells + 4,vector<double>(3));
	
	std::vector<std::vector<double>> deltax;
	deltax.resize(nCells + 4, vector<double>(3));
	
	for(int i = 2; i < nCells + 2; i++){
		double x = x0 + (i - 1.5) * dx;
		if(x <= 0.5){
			u_pri[i][0] = 1000;
			u_pri[i][1] = -100;
			u_pri[i][2] = 2 * 101325;
		}else if(x > 0.5){
			u_pri[i][0] = 1000;
			u_pri[i][1] = 100;
			u_pri[i][2] = 2 * 101325;
		}
	}
	
//	for(int i = 2; i < nCells + 2; i++){
//		double x = x0 + (i - 1.5) * dx;
//		if(x <= 0.5){
//			u_pri[i][0] = 1;
//			u_pri[i][1] = 0;
//			u_pri[i][2] = 1;
//		}else if(x > 0.5){
//			u_pri[i][0] = 0.125;
//			u_pri[i][1] = 0;
//			u_pri[i][2] = 0.1;
//		}
//	}

	do{
		for(int k = 0; k < 3; k++){
			u_pri[1][k] = u_pri[2][k];
			u_pri[0][k] = u_pri[2][k];
			u_pri[nCells + 2][k] = u_pri[nCells + 1][k];
			u_pri[nCells + 3][k] = u_pri[nCells + 1][k];
		}
		
		for(int i = 0; i < nCells + 4; i++){
			u_con[i] = pri2con(u_pri[i], gama, p_inf);
		}
		
		double a_max = -1e8;
		
		for(int i = 0; i < nCells + 4; i++){
			double cs = calcs(u_pri[i], gama, p_inf);
			if(a_max < fabs(u_pri[i][1]) + cs){
				a_max = fabs(u_pri[i][1]) + cs;
			}
		}
		
		double dt = c * dx / a_max;
		t = t + dt;
		
		for(int i = 1; i < nCells + 3; i++){
			for(int k = 0; k < 3; k++){
				deltax0[i][k] = u_con[i][k] - u_con[i - 1][k];
				deltax1[i][k] = u_con[i + 1][k] - u_con[i][k];
			}
		}
		
		for(int i = 1; i < nCells + 3; i++){
			for(int k = 0; k < 3; k++){
				deltax[i][k] = 0.5 * (1 + wx[k]) * deltax0[i][k] + 0.5 * (1 - wx[k]) * deltax1[i][k];
			}
		}

		std::vector<std::vector<double>> xL;
		xL.resize(nCells + 4,vector<double>(3));
		
		std::vector<std::vector<double>> xR;
		xR.resize(nCells + 4,vector<double>(3));
			
		for(int i = 1; i < nCells + 3; i++){
			std::vector<double> rx;
			rx.resize(3);
			for(int k = 0; k < 3; k++){
				rx[k] = (u_con[i][k] - u_con[i - 1][k])/(u_con[i + 1][k] - u_con[i][k]);
				if(u_con[i + 1][k] - u_con[i][k] == 0){
					rx[k] = 0;
				}
				xL[i][k] = u_con[i][k] - 0.5 * Minbee(rx[k]) * deltax[i][k];
				xR[i][k] = u_con[i][k] + 0.5 * Minbee(rx[k]) * deltax[i][k];
			}
		}
		
		std::vector<std::vector<double>> u_bar;
		u_bar.resize(nCells + 4,vector<double>(3));
		
		std::vector<std::vector<double>> halfxL;
		halfxL.resize(nCells + 4,vector<double>(3));
		
		std::vector<std::vector<double>> halfxR;
		halfxR.resize(nCells + 4,vector<double>(3));
		
		for(int i = 1; i < nCells + 3; i++){
			std::vector<double> fxL = f(xL[i], gama, p_inf);
			std::vector<double> fxR = f(xR[i], gama, p_inf);
				
			for(int k = 0; k < 3; k++){
				halfxL[i][k] = xL[i][k] - 0.5 * (dt/dx) * (fxR[k] - fxL[k]);
				halfxR[i][k] = xR[i][k] - 0.5 * (dt/dx) * (fxR[k] - fxL[k]);
			}
		}
		
		std::vector<std::vector<double>> flux_x;
		flux_x.resize(nCells + 3,vector<double>(3));
		
		for(int i = 1; i < nCells + 2; i++){
			flux_x[i] = FORCE_flux(halfxR[i], halfxL[i + 1], dt, dx, gama, p_inf);
		}
		
		for(int i = 2; i < nCells + 2; i++){
			for(int k = 0; k < 3; k++){
				u_bar[i][k] = u_con[i][k] - (dt/dx) * (flux_x[i][k] - flux_x[i - 1][k]);
			}
		}
		
		u_con = u_bar;
		
		for(int i = 2; i < nCells + 2; i++){
			u_pri[i] = con2pri(u_con[i], gama, p_inf);
		}
		std::cout << "stop" << std::endl;
	}while(t < time);
	
	ofstream outFile("/Users/chenzefeng/Desktop/velocity.dat");
	for(int i = 2; i < nCells + 2; i++){
		outFile << x0 + dx * (i - 0.5 - 2) << " " << u_pri[i][0] << endl;
	}
	outFile.close();
}