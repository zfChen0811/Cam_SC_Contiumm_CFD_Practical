#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>
#include <fstream>

std::vector<double> fr(std::vector<double> con, double gama);
double Minbee(double r);
std::vector<double> f(std::vector<double>con, double gama);
std::vector<double> g(std::vector<double>con, double gama);
std::vector<double> con2pri(std::vector<double>con, double gama);
std::vector<double> pri2con(std::vector<double>pri, double gama);
std::vector<double> Fri_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama);
std::vector<double> RI_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama);
std::vector<double> FORCE_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama);
double calecon(std::vector<double> con);
double calepri(std::vector<double> pri, double gama);
double calcs(std::vector<double> pri, double gama);


std::vector<double> fr(std::vector<double> con, double gama){
	std::vector<double> pri = con2pri(con, gama);
	std::vector<double> value;
	value.resize(3);
	value[0] = con[1];
	value[1] = con[1] * pri[1];
	value[2] = (con[2] + pri[2]) * pri[1];
	return value;
}

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

std::vector<double> f(std::vector<double>con, double gama){
	// calculate f
	// based on conservative value roh rohv E
	double rho = con[0];
	double v = con[1] / con[0];
	double E = con[2];
	std::vector<double> pri = con2pri(con, gama);
	double p = pri[2];
	std::vector<double> f;
	f.resize(3);
	f[0] = rho * v;
	f[1] = rho * v * v + p;
	f[2] = (E + p) * v;
	return f;
}

std::vector<double> con2pri(std::vector<double>con, double gama){
	// convert conservative (rho rhov E) to primitive (rho v p)
	std::vector<double> pri;
	pri.resize(3);
	// rho v p
	double e = calecon(con);
	pri[0] = con[0];
	pri[1] = con[1] / con[0];
	pri[2] = (gama - 1) * con[0] * e;
	return pri;
}

std::vector<double> pri2con(std::vector<double>pri, double gama){
	// convert primitive (rho v p) to conservative  (rho rhov E)
	std::vector<double> con;
	con.resize(3);
	double e = calepri(pri, gama);
	// rho rhov E
	con[0] = pri[0];
	con[1] = pri[0] * pri[1];
	con[2] = pri[0] * e + 0.5 * pri[0] * pri[1] * pri[1];
	return con;
}

std::vector<double> Fri_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama){
	// calculate the flux for lax-friedrich
	std::vector<double> Fri;
	Fri.resize(3);
	std::vector<double> f0 = f(con0, gama);
	std::vector<double> f1 = f(con1, gama);
	for(int i = 0; i < 3; i++){
		Fri[i] = 0.5 * (dx/dt) * (con0[i] - con1[i]) + 0.5 * (f1[i] + f0[i]);
	}
	
	return Fri;
}

std::vector<double> RI_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama){
	std::vector<double> RI;
	std::vector<double> u;
	u.resize(3);
	std::vector<double> f0 = f(con0, gama);
	std::vector<double> f1 = f(con1, gama);
	for(int i = 0; i < 3; i++){
		u[i] = 0.5 * (con0[i] + con1[i]) - 0.5 * (dt/dx) * (f1[i] - f0[i]);
	}
	RI = f(u, gama);
	return RI;
}

std::vector<double> FORCE_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama){
	std::vector<double> flux;
	flux.resize(3);
	std::vector<double> RI = RI_flux(con0, con1, dt, dx, gama);
	std::vector<double> Fri = Fri_flux(con0, con1, dt, dx, gama);
	for(int i = 0; i < 3 ; i++){
		flux[i] = 0.5 * (Fri[i] + RI[i]);
	}
	return flux;
}

double calecon(std::vector<double> con){
	double rho = con[0];
	double v = con[1] / con[0];
	double E = con[2];
	double e = (E - 0.5 * rho * v * v)/rho;
	return e;
}

double calepri(std::vector<double> pri, double gama){
	double e = pri[2]/((gama - 1) * pri[0]);
	return e;
}

double calcs(std::vector<double> pri, double gama){
	double cs = sqrt(gama *  pri[2]/ pri[0]);
	return cs;
}

using namespace std;
int main(int argc, char *argv[]) {
	double x0 = 0;
	double x1 = 1;
	double gama = 1.4;
	double c = 0.9;
	double nCells = 100;
	double nPoints = nCells + 1;
	double time = 0.25;
	double t = 0;
	double D = 3;
	double nxPoints = nCells + 1;
	double alpha = D - 1;
	double dx = (x1 - x0)/(nxPoints-1);
	
	std::vector<std::vector<double>> u_pri;
	u_pri.resize(nCells + 3, vector<double>(3));

	std::vector<std::vector<double>> u_con;
	u_con.resize(nCells + 3, vector<double>(3));
	
	std::vector<double> wx;
	wx.resize(3,0);
	
	
	std::vector<std::vector<double>> deltax0;
	deltax0.resize(nCells + 3, vector<double>(3));
	
	std::vector<std::vector<double>> deltax1;
	deltax1.resize(nCells + 3, vector<double>(3));
	
	std::vector<std::vector<double>> deltax;
	deltax.resize(nCells + 3, vector<double>(3));
		
	for(int i = 2; i < nCells + 2; i++){
		double x = x0 + (i - 0.5 - 1) * dx;
		if(x < 0.4){
			u_pri[i][0] = 1;
			u_pri[i][1] = 0;
			u_pri[i][2] = 1;
		}else if(x >= 0.4){
			u_pri[i][0] = 0.125;
			u_pri[i][1] = 2;
			u_pri[i][2] = 0.1;
		}
	}

	do{
		u_pri[0][0] = u_pri[3][0];
		u_pri[1][0] = u_pri[2][0];
		u_pri[0][1] = -u_pri[3][1];
		u_pri[1][1] = -u_pri[2][1];
		u_pri[0][2] = u_pri[3][2];
		u_pri[1][2] = u_pri[2][2];
		
		u_pri[nCells + 2] = u_pri[nCells + 1];

		double a_max;
		for(int i = 0; i < nCells + 3; i++){
			if(i == 0){
				a_max = fabs(u_pri[i][1]) + calcs(u_pri[i], gama);
			}else if(a_max < fabs(u_pri[i][1]) + calcs(u_pri[i], gama)){
				a_max = fabs(u_pri[i][1]) + calcs(u_pri[i], gama);
			}
		}
		
		double dt = c * dx / a_max;
		t = t + dt;
		
		for(int i = 0; i < nCells + 3; i++){
			u_con[i] = pri2con(u_pri[i], gama);
		}

		
		for(int i = 2; i < nCells + 2; i++){
			for(int j = 0; j < 3; j++){
				deltax0[i][j] = u_con[i][j] - u_con[i - 1][j];
				deltax1[i][j] = u_con[i + 1][j] - u_con[i][j];
			}
		}

		for(int i = 2; i < nCells + 2; i++){
			for(int j = 0; j < 3; j++){
				deltax[i][j] = 0.5 * (1 + wx[j]) * deltax0[i][j] + 0.5 * (1 - wx[j]) * deltax1[i][j];
			}
		}

		deltax[0][0] = deltax[3][0];
		deltax[1][0] = deltax[2][0];
		deltax[0][1] = -deltax[3][1];
		deltax[1][1] = -deltax[2][1];
		deltax[0][2] = deltax[3][2];
		deltax[1][2] = deltax[2][2];
		
		deltax[nCells + 2] = deltax[nCells + 1];
		
		std::vector<std::vector<double>> xL;
		xL.resize(nCells + 3,vector<double>(3));
		
		std::vector<std::vector<double>> xR;
		xR.resize(nCells + 3,vector<double>(3));
			
		for(int i = 2; i < nCells + 2; i++){
			std::vector<double> rx;
			rx.resize(3);
			for(int j = 0; j < 3; j++){
				rx[j] = (u_con[i][j] - u_con[i - 1][j])/(u_con[i + 1][j] - u_con[i][j]);
				if(u_con[i + 1][j] - u_con[i][j] == 0){
					rx[j] = 0;
				}
				xL[i][j] = u_con[i][j] - 0.5 * Minbee(rx[j]) * deltax[i][j];
				xR[i][j] = u_con[i][j] + 0.5 * Minbee(rx[j]) * deltax[i][j];
			}
		}
		
		xL[0][0] = xL[3][0];
		xL[1][0] = xL[2][0];
		xL[0][1] = -xL[3][1];
		xL[1][1] = -xL[2][1];
		xL[0][2] = xL[3][2];
		xL[1][2] = xL[2][2];
		
		xL[nCells + 2] = xL[nCells + 1];
		
		xR[0][0] = xR[3][0];
		xR[1][0] = xR[2][0];
		xR[0][1] = -xR[3][1];
		xR[1][1] = -xR[2][1];
		xR[0][2] = xR[3][2];
		xR[1][2] = xR[2][2];
		
		xR[nCells + 2] = xR[nCells + 1];
		
		
		std::vector<std::vector<double>> u_con_bar;
		u_con_bar.resize(nCells + 3,vector<double>(3));
		
		std::vector<std::vector<double>> halfxL;
		halfxL.resize(nCells + 3,vector<double>(3));
		
		std::vector<std::vector<double>> halfxR;
		halfxR.resize(nCells + 3,vector<double>(3));
		
		for(int i = 2; i < nCells + 2; i++){
			std::vector<double> fxL = f(xL[i], gama);
			std::vector<double> fxR = f(xR[i], gama);
			
			for(int j = 0; j < 3; j++){
				halfxL[i][j] = xL[i][j] - 0.5 * (dt/dx) * (fxR[j] - fxL[j]);
				halfxR[i][j] = xR[i][j] - 0.5 * (dt/dx) * (fxR[j] - fxL[j]);
			}
		}
		
		halfxL[0][0] = halfxL[3][0];
		halfxL[1][0] = halfxL[2][0];
		halfxL[0][1] = -halfxL[3][1];
		halfxL[1][1] = -halfxL[2][1];
		halfxL[0][2] = halfxL[3][2];
		halfxL[1][2] = halfxL[2][2];
		
		halfxL[nCells + 2] = halfxL[nCells + 1];
		
		halfxR[0][0] = halfxR[3][0];
		halfxR[1][0] = halfxR[2][0];
		halfxR[0][1] = -halfxR[3][1];
		halfxR[1][1] = -halfxR[2][1];
		halfxR[0][2] = halfxR[3][2];
		halfxR[1][2] = halfxR[2][2];
		
		halfxR[nCells + 2] = halfxR[nCells + 1];
		
		std::vector<std::vector<double>> flux_x;
		flux_x.resize(nCells + 3,vector<double>(3));
		
		for(int i = 0; i < nCells + 2; i++){
			flux_x[i] = FORCE_flux(halfxR[i], halfxL[i + 1], dt, dx, gama);
		}
		
		for(int i = 2; i < nCells + 2; i++){
			for(int j = 0; j < 3; j++){
				u_con_bar[i][j] = u_con[i][j] - (dt/dx) * (flux_x[i][j] - flux_x[i - 1][j]);
			}
		}
		
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		std::vector<std::vector<double>> K1;
		std::vector<std::vector<double>> K2;
		std::vector<std::vector<double>> K3;
		std::vector<std::vector<double>> K4;
		K1.resize(nCells + 3, std::vector<double> (3)); 
		K2.resize(nCells + 3, std::vector<double> (3)); 
		K3.resize(nCells + 3, std::vector<double> (3)); 
		K4.resize(nCells + 3, std::vector<double> (3)); 
		
		std::vector<std::vector<double>> u_con_plus;
		u_con_plus.resize(nCells + 3, std::vector<double> (3));
		
		for(int i = 2; i < nCells + 2; i++){
			double x = x0 + dx * (i - 0.5 - 1);
			std::vector<double> u;
			std::vector<double> s;
			s.resize(3);
			
			u = u_con_bar[i];
			
			std::vector<double> fr_value = fr(u, gama);
			
			for(int j = 0; j < 3; j++){
				s[j] = (-alpha/x) * fr_value[j];
				K1[i][j] = dt * s[j];
			}
		}
		
		for(int i = 2; i < nCells + 2; i++){
			double x = x0 + dx * (i - 0.5 - 1);
			std::vector<double> u;
			u.resize(3);
			std::vector<double> s;
			s.resize(3);
			
			for(int j = 0; j < 3; j++){
				u[j] = u_con_bar[i][j] + 0.5 * K1[i][j];
			}
			
			std::vector<double> fr_value = fr(u, gama);
			
			for(int j = 0; j < 3; j++){
				s[j] = (-alpha/x) * fr_value[j];
				K2[i][j] = dt * s[j];
			}
		}
		
		for(int i = 2; i < nCells + 2; i++){
			double x = x0 + dx * (i - 0.5 - 1);
			std::vector<double> u;
			u.resize(3);
			std::vector<double> s;
			s.resize(3);
			
			for(int j = 0; j < 3; j++){
				u[j] = u_con_bar[i][j] + 0.5 * K2[i][j];
			}
			
			std::vector<double> fr_value = fr(u, gama);
			
			for(int j = 0; j < 3; j++){
				s[j] = (-alpha/x) * fr_value[j];
				K3[i][j] = dt * s[j];
			}
		}
		
		for(int i = 2; i < nCells + 2; i++){
			double x = x0 + dx * (i - 0.5 - 1);
			std::vector<double> u;
			u.resize(3);
			std::vector<double> s;
			s.resize(3);
			
			for(int j = 0; j < 3; j++){
				u[j] = u_con_bar[i][j] + K3[i][j];
			}
			
			std::vector<double> fr_value = fr(u, gama);
			
			for(int j = 0; j < 3; j++){
				s[j] = (-alpha/x) * fr_value[j];
				K4[i][j] = dt * s[j];
			}
		}
		
		for(int i = 2; i < nCells + 2; i++){
			for(int j = 0; j < 3; j++){
				u_con_plus[i][j] = u_con_bar[i][j] + (1.0/6.0) * (K1[i][j] + 2 * K2[i][j] + 2 * K3[i][j] + K4[i][j]);
			}
		}
		
		
		u_con = u_con_plus;
		
		for(int i = 2; i < nCells + 2; i++){
			u_pri[i] = con2pri(u_con[i], gama);
		}
		
	}while(t < time);
	
	ofstream outFile("/Users/chenzefeng/Desktop/rho.dat");
	for(int i = 2; i < nCells + 2; i++){
		outFile << x0 + dx * (i - 0.5 - 1) << " " << " " << u_pri[i][0] << endl;
	}
	outFile.close();
}