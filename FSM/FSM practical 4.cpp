#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>
#include <fstream>

double Minbee(double r);
std::vector<double> f(std::vector<double>con, double gama);
std::vector<double> con2pri(std::vector<double>con, double gama);
std::vector<double> pri2con(std::vector<double>pri, double gama);
std::vector<double> Fri_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama);
std::vector<double> RI_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama);
std::vector<double> FORCE_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama);
double calecon(std::vector<double> con);
double calepri(std::vector<double> pri, double gama);
double calcs(std::vector<double> pri, double gama);

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
	double vx = con[1] / con[0];
	double E = con[2];
	std::vector<double> pri = con2pri(con, gama);
	double p = pri[2];
	std::vector<double> f;
	f.resize(3);
	f[0] = rho * vx;
	f[1] = rho * vx * vx + p;
	f[2] = (E + p) * vx;
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
	pri[2] = (gama - 1) * pri[0] * e;
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
	con[2] = pri[0] * e + 0.5 * pri[0] * (pri[1] * pri[1]);
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
	double vx = con[1] / con[0];
	double E = con[2];
	double e = (E - 0.5 * rho * (vx * vx))/rho;
	return e;
}

double calepri(std::vector<double> pri, double gama){
	double e = pri[2]/((gama - 1) * pri[0]);
	return e;
}

double calcs(std::vector<double> pri, double gama){
	double cs = sqrt(gama *  pri[2] / pri[0]);
	return cs;
}

using namespace std;
int main(int argc, char *argv[]) {
	double x0 = 0;
	double x1 = 1;
	double gama1 = 1.4;
	double gama2 = 1.67;
	double c = 1;
	int nCells = 300;
	double nxPoints = nCells + 1;
	double time = 0.3;
	double t = 0;
	double dx = (x1 - x0)/(nxPoints-2);
	
	std::vector<double> phi;
	phi.resize(nCells + 6,0);
	
	std::vector<std::vector<double>> u_pri;
	u_pri.resize(nCells + 6,vector<double>(3));
	
	std::vector<std::vector<double>> u_con;
	u_con.resize(nCells + 6,vector<double>(3));
	
	std::vector<std::vector<double>> u1_pri;
	u1_pri.resize(nCells + 6,vector<double>(3));

	std::vector<std::vector<double>> u1_con;
	u1_con.resize(nCells + 6,vector<double>(3));
	
	std::vector<std::vector<double>> u2_pri;
	u2_pri.resize(nCells + 6,vector<double>(3));
	
	std::vector<std::vector<double>> u2_con;
	u2_con.resize(nCells + 6,vector<double>(3));
	
	std::vector<double> w1lx;
	w1lx.resize(3,0);
	
	std::vector<std::vector<double>> delta1lx0;
	delta1lx0.resize(nCells + 6,vector<double>(3));
	
	std::vector<std::vector<double>> delta1lx1;
	delta1lx1.resize(nCells + 6,vector<double>(3));
	
	std::vector<std::vector<double>> delta1lx;
	delta1lx.resize(nCells + 6, vector<double>(3));
	
	std::vector<double> w1rx;
	w1rx.resize(3,0);
	
	std::vector<std::vector<double>> delta1rx0;
	delta1rx0.resize(nCells + 6,vector<double>(3));
	
	std::vector<std::vector<double>> delta1rx1;
	delta1rx1.resize(nCells + 6,vector<double>(3));
	
	std::vector<std::vector<double>> delta1rx;
	delta1rx.resize(nCells + 6, vector<double>(3));
	
	std::vector<double> w2x;
	w2x.resize(3,0);
	
	std::vector<std::vector<double>> delta2x0;
	delta2x0.resize(nCells + 6,vector<double>(3));
	
	std::vector<std::vector<double>> delta2x1;
	delta2x1.resize(nCells + 6,vector<double>(3));
	
	std::vector<std::vector<double>> delta2x;
	delta2x.resize(nCells + 6, vector<double>(3));
	
	double interface_left_bar, interface_right_bar;
	
	
	for(int i = 3; i < nCells + 3; i++){
		double x = x0 + (i - 2.5) * dx;
		if(x <= 0.25){
			u1_pri[i][0] = 1.3765;
			u1_pri[i][1] = 0.3948;
			u1_pri[i][2] = 1.57;
		}
		if(x >= 0.25){
			u1_pri[i][0] = 1;
			u1_pri[i][1] = 0;
			u1_pri[i][2] = 1;
		}
		if(x > 0.4 && x < 0.6){
			u2_pri[i][0] = 0.138;
			u2_pri[i][1] = 0;
			u2_pri[i][2] = 1;
		}
		phi[i] = 0.1 - fabs(x - 0.5);
	}
	
	std::vector<double> phi_bar;
	phi_bar.resize(nCells + 6, 0);
	
	std::vector<double> phi_bar_1;
	phi_bar_1.resize(nCells + 6, 0);
	
	
	std::vector<double> phi_bar_2;
	phi_bar_2.resize(nCells + 6, 0);
	
	std::vector<double> phi_bar_3;
	phi_bar_3.resize(nCells + 6, 0);
	
	do{
		double interface_left;
		double interface_right;
		
		for(int i = 3; i < nCells + 2; i++){
			if(phi[i] * phi[i + 1] <= 0){
				interface_left = i;
				break;
			}
		}
		
		for(int i = interface_left + 1; i < nCells + 2; i++){
			if(phi[i] * phi[i + 1] <= 0){
				interface_right = i;
				break;
			}
		}
		
		phi[0] = phi[3];
		phi[1] = phi[3];
		phi[2] = phi[3];
		phi[nCells + 3] = phi[nCells + 2];
		phi[nCells + 4] = phi[nCells + 2];
		phi[nCells + 5] = phi[nCells + 2];
		
		// interface_index + 1
		for(int i = interface_left + 1; i < interface_left + 4; i++){
			u1_pri[i][1] = u2_pri[i][1];
			u1_pri[i][2] = u2_pri[i][2];
			u1_pri[i][0] = u1_pri[interface_left][0] * pow(u1_pri[i][2]/u1_pri[interface_left][2],1.0/gama1);
		}
		
		for(int i = interface_left - 2; i < interface_left + 1; i++){
			u2_pri[i][1] = u1_pri[i][1];
			u2_pri[i][2] = u1_pri[i][2];
			u2_pri[i][0] = u2_pri[interface_left + 1][0] * pow(u2_pri[i][2]/u2_pri[interface_left + 1][2],1.0/gama2);
		}
		
		for(int i = interface_right + 1; i < interface_right + 4; i++){
			u2_pri[i][1] = u1_pri[i][1];
			u2_pri[i][2] = u1_pri[i][2];
			u2_pri[i][0] = u2_pri[interface_right][0] * pow(u2_pri[i][2]/u2_pri[interface_right][2],1.0/gama2);
		}
		
		for(int i = interface_right - 2; i < interface_right + 1; i++){
			u1_pri[i][1] = u2_pri[i][1];
			u1_pri[i][2] = u2_pri[i][2];
			u1_pri[i][0] = u1_pri[interface_right + 1][0] * pow(u1_pri[i][2]/u1_pri[interface_right + 1][2],1.0/gama1);
		}
		
		for(int k = 0; k < 3; k++){
			u1_pri[2][k] = u1_pri[3][k];
			u1_pri[1][k] = u1_pri[4][k];
			u1_pri[0][k] = u1_pri[5][k];
			
			u1_pri[nCells + 3][k] = u1_pri[nCells + 2][k];
			u1_pri[nCells + 4][k] = u1_pri[nCells + 1][k];
			u1_pri[nCells + 5][k] = u1_pri[nCells][k];
		}
		
		for(int i = 0; i < interface_left + 4; i++){
			u1_con[i] = pri2con(u1_pri[i], gama1);
		}
		
		for(int i = interface_left - 2; i < interface_right + 4; i++){
			u2_con[i] = pri2con(u2_pri[i], gama2);
		}
		
		for(int i = interface_right - 2; i < nCells + 6; i++){
			u1_con[i] = pri2con(u1_pri[i], gama1);
		}
		
		double a_max = -1e8;
		
		for(int i = 0; i < interface_left + 4; i++){
			double cs1 = calcs(u1_pri[i], gama1);
			if(a_max < fabs(u1_pri[i][1]) + cs1){
				a_max = fabs(u1_pri[i][1]) + cs1;
			}
		}
		
		for(int i = interface_left - 2; i < interface_right + 4; i++){
			double cs2 = calcs(u2_pri[i], gama2);
			if(a_max < fabs(u2_pri[i][1]) + cs2){
				a_max = fabs(u2_pri[i][1]) + cs2;
			}
		}
		
		for(int i = interface_right - 2; i < nCells + 6; i++){
			double cs1 = calcs(u1_pri[i], gama2);
			if(a_max < fabs(u1_pri[i][1]) + cs1){
				a_max = fabs(u1_pri[i][1]) + cs1;
			}
		}
		
		double dt = c * dx / a_max;
		t = t + dt;

		std::vector<double> phi_bar;
		phi_bar.resize(nCells + 6, 0);

		for(int i = 3; i < nCells + 3; i++){
			if(phi[i] < 0){
				if(u1_con[i][1] > 0){
					phi_bar[i] = phi[i] - (dt / dx) * u1_pri[i][1] * (phi[i] - phi[i - 1]);
				}else{
					phi_bar[i] = phi[i] - (dt / dx) * u1_pri[i][1] * (phi[i + 1] - phi[i]);
				}
			}else{
				if(u2_con[i][1] > 0){
					phi_bar[i] = phi[i] - (dt / dx) * u2_pri[i][1] * (phi[i] - phi[i - 1]);
				}else{
					phi_bar[i] = phi[i] - (dt / dx) * u2_pri[i][1] * (phi[i + 1] - phi[i]);
				}
			}
		}
		
		phi_bar[0] = phi_bar[3];
		phi_bar[1] = phi_bar[3];
		phi_bar[2] = phi_bar[3];
		phi_bar[nCells + 3] = phi_bar[nCells + 2];
		phi_bar[nCells + 4] = phi_bar[nCells + 2];
		phi_bar[nCells + 5] = phi_bar[nCells + 2];
		
		for(int i = 1; i < interface_left + 3; i++){
			for(int k = 0; k < 3; k++){
				delta1lx0[i][k] = u1_con[i][k] - u1_con[i - 1][k];
				delta1lx1[i][k] = u1_con[i + 1][k] - u1_con[i][k];
			}
		}
		
		for(int i = interface_left - 1; i < interface_right + 3; i++){
			for(int k = 0; k < 3; k++){
				delta2x0[i][k] = u2_con[i][k] - u2_con[i - 1][k];
				delta2x1[i][k] = u2_con[i + 1][k] - u2_con[i][k];
			}
		}
		
		for(int i = interface_right - 1; i < nCells + 5; i++){
			for(int k = 0; k < 3; k++){
				delta1rx0[i][k] = u1_con[i][k] - u1_con[i - 1][k];
				delta1rx1[i][k] = u1_con[i + 1][k] - u1_con[i][k];
			}
		}
		
		for(int i = 1; i < interface_left + 3; i++){
			for(int k = 0; k < 3; k++){
				delta1lx[i][k] = 0.5 * (1 + w1lx[k]) * delta1lx0[i][k] + 0.5 * (1 - w1lx[k]) * delta1lx1[i][k];
			}
		}
		
		for(int i = interface_left - 1; i < interface_right + 3; i++){
			for(int k = 0; k < 3; k++){
				delta2x[i][k] = 0.5 * (1 + w2x[k]) * delta2x0[i][k] + 0.5 * (1 - w2x[k]) * delta2x1[i][k];
			}
		}
		
		for(int i = interface_right - 1; i < nCells + 5; i++){
			for(int k = 0; k < 3; k++){
				delta1rx[i][k] = 0.5 * (1 + w1rx[k]) * delta1rx0[i][k] + 0.5 * (1 - w1rx[k]) * delta1rx1[i][k];
			}
		}
		
		std::vector<std::vector<double>> x1lL;
		x1lL.resize(nCells + 6,vector<double>(3));
		
		std::vector<std::vector<double>> x1lR;
		x1lR.resize(nCells + 6,vector<double>(3));
		
		std::vector<std::vector<double>> x2L;
		x2L.resize(nCells + 6,vector<double>(3));
		
		std::vector<std::vector<double>> x2R;
		x2R.resize(nCells + 6,vector<double>(3));
		
		std::vector<std::vector<double>> x1rL;
		x1rL.resize(nCells + 6,vector<double>(3));
		
		std::vector<std::vector<double>> x1rR;
		x1rR.resize(nCells + 6,vector<double>(3));
		
		for(int i = 1; i < interface_left + 3; i++){
			std::vector<double> r1x;
			r1x.resize(3);
			for(int k = 0; k < 3; k++){
				r1x[k] = (u1_con[i][k] - u1_con[i - 1][k])/(u1_con[i + 1][k] - u1_con[i][k]);
				if(u1_con[i + 1][k] - u1_con[i][k] == 0){
					r1x[k] = 0;
				}
				x1lL[i][k] = u1_con[i][k] - 0.5 * Minbee(r1x[k]) * delta1lx[i][k];
				x1lR[i][k] = u1_con[i][k] + 0.5 * Minbee(r1x[k]) * delta1lx[i][k];
			}
		}
		
		for(int i = interface_left - 1; i < interface_right + 3; i++){
			std::vector<double> r2x;
			r2x.resize(3);
			for(int k = 0; k < 3; k++){
				r2x[k] = (u2_con[i][k] - u2_con[i - 1][k])/(u2_con[i + 1][k] - u2_con[i][k]);
				if(u2_con[i + 1][k] - u2_con[i][k] == 0){
					r2x[k] = 0;
				}
				x2L[i][k] = u2_con[i][k] - 0.5 * Minbee(r2x[k]) * delta2x[i][k];
				x2R[i][k] = u2_con[i][k] + 0.5 * Minbee(r2x[k]) * delta2x[i][k];
			}
		}
		
		for(int i = interface_right - 1; i < nCells + 5; i++){
			std::vector<double> r1x;
			r1x.resize(3);
			for(int k = 0; k < 3; k++){
				r1x[k] = (u1_con[i][k] - u1_con[i - 1][k])/(u1_con[i + 1][k] - u1_con[i][k]);
				if(u1_con[i + 1][k] - u1_con[i][k] == 0){
					r1x[k] = 0;
				}
				x1rL[i][k] = u1_con[i][k] - 0.5 * Minbee(r1x[k]) * delta1rx[i][k];
				x1rR[i][k] = u1_con[i][k] + 0.5 * Minbee(r1x[k]) * delta1rx[i][k];
			}
		}
		
		
		std::vector<std::vector<double>> u_bar;
		u_bar.resize(nCells + 6,vector<double>(3));
		
		std::vector<std::vector<double>> u1_bar;
		u1_bar.resize(nCells + 6,vector<double>(3));
		
		std::vector<std::vector<double>> half1lxL;
		half1lxL.resize(nCells + 6,vector<double>(3));
		
		std::vector<std::vector<double>> half1lxR;
		half1lxR.resize(nCells + 6,vector<double>(3));
		
		std::vector<std::vector<double>> half1rxL;
		half1rxL.resize(nCells + 6,vector<double>(3));
		
		std::vector<std::vector<double>> half1rxR;
		half1rxR.resize(nCells + 6,vector<double>(3));
		
		std::vector<std::vector<double>> u2_bar;
		u2_bar.resize(nCells + 6,vector<double>(3));
		
		std::vector<std::vector<double>> half2xL;
		half2xL.resize(nCells + 6,vector<double>(3));
		
		std::vector<std::vector<double>> half2xR;
		half2xR.resize(nCells + 6,vector<double>(3));
		
		for(int i = 1; i < interface_left + 3; i++){
			std::vector<double> f1xL = f(x1lL[i], gama1);
			std::vector<double> f1xR = f(x1lR[i], gama1);
			for(int k = 0; k < 3; k++){
				half1lxL[i][k] = x1lL[i][k] - 0.5 * (dt/dx) * (f1xR[k] - f1xL[k]);
				half1lxR[i][k] = x1lR[i][k] - 0.5 * (dt/dx) * (f1xR[k] - f1xL[k]);
			}
		}
		
		for(int i = interface_left - 1; i < interface_right + 3; i++){
			std::vector<double> f2xL = f(x2L[i], gama2);
			std::vector<double> f2xR = f(x2R[i], gama2);
			for(int k = 0; k < 3; k++){
				half2xL[i][k] = x2L[i][k] - 0.5 * (dt/dx) * (f2xR[k] - f2xL[k]);
				half2xR[i][k] = x2R[i][k] - 0.5 * (dt/dx) * (f2xR[k] - f2xL[k]);
			}
		}
		
		for(int i = interface_right - 1; i < nCells + 5; i++){
			std::vector<double> f1xL = f(x1rL[i], gama1);
			std::vector<double> f1xR = f(x1rR[i], gama1);
			for(int k = 0; k < 3; k++){
				half1rxL[i][k] = x1rL[i][k] - 0.5 * (dt/dx) * (f1xR[k] - f1xL[k]);
				half1rxR[i][k] = x1rR[i][k] - 0.5 * (dt/dx) * (f1xR[k] - f1xL[k]);
			}
		}
		
		std::vector<std::vector<double>> flux1l_x;
		flux1l_x.resize(nCells + 5,vector<double>(3));
		
		std::vector<std::vector<double>> flux2_x;
		flux2_x.resize(nCells + 5,vector<double>(3));
		
		std::vector<std::vector<double>> flux1r_x;
		flux1r_x.resize(nCells + 5,vector<double>(3));
		
		for(int i = 1; i < interface_left + 3; i++){
			flux1l_x[i] = FORCE_flux(half1lxR[i], half1lxL[i + 1], dt, dx, gama1);
		}
		
		for(int i = interface_left - 1; i < interface_right + 3; i++){
			flux2_x[i] = FORCE_flux(half2xR[i], half2xL[i + 1], dt, dx, gama2);
		}

		for(int i = interface_right - 1; i < nCells + 5; i++){
			flux1r_x[i] = FORCE_flux(half1rxR[i], half1rxL[i + 1], dt, dx, gama1);
		}
		
		for(int i = 2; i < interface_left + 2; i++){
			for(int k = 0; k < 3; k++){
				u1_bar[i][k] = u1_con[i][k] - (dt/dx) * (flux1l_x[i][k] - flux1l_x[i - 1][k]);
			}
		}
		
		for(int i = interface_left; i < interface_right + 2; i++){
			for(int k = 0; k < 3; k++){
				u2_bar[i][k] = u2_con[i][k] - (dt/dx) * (flux2_x[i][k] - flux2_x[i - 1][k]);
			}
		}
		
		for(int i = interface_right; i < nCells + 4; i++){
			for(int k = 0; k < 3; k++){
				u1_bar[i][k] = u1_con[i][k] - (dt/dx) * (flux1r_x[i][k] - flux1r_x[i - 1][k]);
			}
		}
		
		
		u1_con = u1_bar;
		u2_con = u2_bar;
		
		
		for(int i = 3; i < nCells + 2; i++){
			if(phi_bar[i] * phi_bar[i + 1] <= 0){
				interface_left_bar = i;
				break;
			}
		}
		
		for(int i = interface_left_bar + 1; i < nCells + 2; i++){
			if(phi_bar[i] * phi_bar[i + 1] <= 0){
				interface_right_bar = i;
				break;
			}
		}
		
		phi_bar_1 = phi_bar;
		phi_bar_2 = phi_bar;
		
		for(int j = 3; j < interface_left_bar + 1; j++){
			phi_bar_1[j] = (x0 + dx * (j - 3)) - (x0 + dx * (interface_left_bar - 3)) + phi_bar[interface_left_bar];
			phi_bar[j] = phi_bar_1[j];
		}
		for(int j = interface_left_bar + 1; j < interface_right_bar + 1; j++){
			phi_bar_1[j] = (x0 + dx * (j - 3)) - (x0 + dx * (interface_left_bar + 1 - 3)) + phi_bar[interface_left_bar + 1];
			phi_bar_2[j] = -((x0 + dx * (j - 3)) - (x0 + dx * (interface_right_bar - 3))) + phi_bar[interface_right_bar];
			phi_bar[j] = fmin(phi_bar_1[j], phi_bar_2[j]);
		}
		
		for(int j = interface_right_bar + 1; j < nCells + 3; j++){
			phi_bar_2[j] = -((x0 + dx * (j - 3)) - (x0 + dx * (interface_right_bar + 1 - 3))) + phi_bar[interface_right_bar + 1];
			phi_bar[j] = phi_bar_2[j];
		}

		phi_bar[0] = phi_bar[3];
		phi_bar[1] = phi_bar[3];
		phi_bar[2] = phi_bar[3];
		phi_bar[nCells + 3] = phi_bar[nCells + 2];
		phi_bar[nCells + 4] = phi_bar[nCells + 2];
		phi_bar[nCells + 5] = phi_bar[nCells + 2];
		
		phi = phi_bar;
		
		for(int i = 3; i < interface_left_bar + 1; i++){
			u1_pri[i] = con2pri(u1_con[i], gama1);
		}
		
		for(int i = interface_left + 1; i < interface_right_bar + 1; i++){
			u2_pri[i] = con2pri(u2_con[i], gama2);
		}
		
		for(int i = interface_right_bar + 1; i < nCells + 3; i++){
			u1_pri[i] = con2pri(u1_con[i], gama1);
		}

		for(int i = 3;  i < nCells + 3; i++){
			if(phi[i] < 0){
				u_pri[i] = u1_pri[i];
			}else{
				u_pri[i] = u2_pri[i];
			}
		}
		
	}while(t < time);
	
	
	ofstream outFile("/Users/chenzefeng/Desktop/velocity.dat");
	for(int i = 3; i < nCells + 3; i++){
		outFile << x0 + dx * (i - 3) << " " << u_pri[i][0] << endl;
	}
	outFile.close();
	
}