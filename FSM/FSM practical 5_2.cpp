#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>
#include <fstream>

std::vector<double> f(std::vector<double>con, double gama, double p_inf);
std::vector<double> con2pri(std::vector<double>con, double gama, double p_inf);
std::vector<double> pri2con(std::vector<double>pri, double gama, double p_inf);
double calecon(std::vector<double> con);
double calepri(std::vector<double> pri, double gama, double p_inf);
double calcs(std::vector<double> pri, double gama, double p_inf);

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

std::vector<double> pri2con(std::vector<double>pri,double gama, double p_inf){
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
	double cs = sqrt(gama *  (pri[2] + p_inf)/ pri[0]);
	return cs;
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


using namespace std;
int main(int argc, char *argv[]) {
	double x0 = 0;
	double x1 = 1;
	double gama1 = 4.4;
	double gama2 = 1.4;
	double p_inf1 = 6 * pow(10, 8);
	double p_inf2 = 0;
	double c = 0.8;
	int nCells = 200;
	double nxPoints = nCells + 1;
	double time = 237.44 * pow(10, -6);
	double t = 0;
	double dx = (x1 - x0)/(nxPoints-1);
	double p_old, p_new;
	double eplison = 10e-8;
	
	std::vector<double> phi;
	phi.resize(nCells + 6,0);
	
	std::vector<std::vector<double>> u_pri;
	u_pri.resize(nCells + 6,vector<double>(3));
	
	std::vector<std::vector<double>> u1_pri;
	u1_pri.resize(nCells + 6,vector<double>(3));
	
	std::vector<std::vector<double>> u1_con;
	u1_con.resize(nCells + 6,vector<double>(3));
	
	std::vector<std::vector<double>> flux1;
	flux1.resize(nCells + 5,vector<double>(3));
	
	std::vector<std::vector<double>> u2_pri;
	u2_pri.resize(nCells + 6,vector<double>(3));
	
	std::vector<std::vector<double>> u2_con;
	u2_con.resize(nCells + 6,vector<double>(3));
	
	std::vector<std::vector<double>> flux2;
	flux2.resize(nCells + 5,vector<double>(3));
	
	std::vector<double> w1x;
	w1x.resize(3,0);
	
	std::vector<std::vector<double>> delta1x0;
	delta1x0.resize(nCells + 6,vector<double>(3));
	
	std::vector<std::vector<double>> delta1x1;
	delta1x1.resize(nCells + 6,vector<double>(3));
	
	std::vector<std::vector<double>> delta1x;
	delta1x.resize(nCells + 6, vector<double>(3));
	
	std::vector<double> w2x;
	w2x.resize(3,0);
	
	std::vector<std::vector<double>> delta2x0;
	delta2x0.resize(nCells + 6,vector<double>(3));
	
	std::vector<std::vector<double>> delta2x1;
	delta2x1.resize(nCells + 6,vector<double>(3));
	
	std::vector<std::vector<double>> delta2x;
	delta2x.resize(nCells + 6, vector<double>(3));
		
	for(int i = 3; i < nCells + 3; i++){
		double x = x0 + (i - 2.5) * dx;
		if(x <= 0.7){
			u1_pri[i][0] = 1000;
			u1_pri[i][1] = 0;
			u1_pri[i][2] = pow(10,9);
		}
		if(x >= 0.7){
			u2_pri[i][0] = 50;
			u2_pri[i][1] = 0;
			u2_pri[i][2] = pow(10,5);
		}
		phi[i] = x - 0.7;
	}

	
	do{
		double interface_index;
		
		for(int i = 3; i < nCells + 2; i++){
			if(phi[i] * phi[i + 1] <= 0){
				interface_index = i;
				break;
			}
		}
		
		for(int i = 3; i < nCells + 2; i++){
			u1_con[i] = pri2con(u1_pri[i], gama1, p_inf1);
			u2_con[i] = pri2con(u2_pri[i], gama2, p_inf2);
		}
		
		std::vector<double> uL = u1_con[interface_index];
		std::vector<double> uR = u2_con[interface_index + 1];
		std::vector<double> priL = u1_pri[interface_index];
		std::vector<double> priR = u2_pri[interface_index + 1];
		std::vector<double> fL = f(uL, gama1, p_inf1);
		std::vector<double> fR = f(uR, gama2, p_inf2);
		double pL = u1_pri[interface_index][2];
		double pR = u2_pri[interface_index + 1][2];
		double rhoL = uL[0];
		double rhoR = uR[0];
		double vL = priL[1];
		double vR = priR[1];
		double EL = uL[2];
		double ER = uR[2];
		double csL = calcs(priL, gama1, p_inf1);
		double csR = calcs(priR, gama2, p_inf2);
		double S_plus = fmax(fabs(vL) + csL, fabs(vR) + csR);
		double SL = -S_plus;
		double SR = S_plus;
		std::vector<double> uL_HLLC;
		uL_HLLC.resize(3,0);
		std::vector<double> uR_HLLC;
		uR_HLLC.resize(3,0);
		
		double S_star = (pR - pL + rhoL * vL * (SL - vL) - rhoR * vR * (SR - vR))/(rhoL * (SL - vL) - rhoR * (SR - vR));
		double v_star = S_star;

		uL_HLLC[0] = rhoL * ((SL - vL)/(SL - S_star)) * 1;
		uL_HLLC[1] = rhoL * ((SL - vL)/(SL - S_star)) * S_star;
		uL_HLLC[2] = rhoL * ((SL - vL)/(SL - S_star)) * ((EL / rhoL) + (S_star - vL) * (S_star + pL / (rhoL * (SL - vL))));
		
		uR_HLLC[0] = rhoR * ((SR - vR)/(SR - S_star)) * 1;
		uR_HLLC[1] = rhoR * ((SR - vR)/(SR - S_star)) * S_star;
		uR_HLLC[2] = rhoR * ((SR - vR)/(SR - S_star)) * ((ER / rhoR) + (S_star - vR) * (S_star + pR / (rhoR * (SR - vR))));
		
		for(int i = interface_index + 1; i < interface_index + 4; i++){
			u1_con[i] = uL_HLLC;
			u1_pri[i] = con2pri(u1_con[i], gama1, p_inf1);
		}
		
		for(int i = interface_index - 2; i < interface_index + 1; i++){
			u2_con[i] = uR_HLLC;
			u2_pri[i] = con2pri(u2_con[i], gama2, p_inf2);
		}
		
		for(int k = 0; k < 3; k++){
			u1_pri[2][k] = u1_pri[3][k];
			u1_pri[1][k] = u1_pri[4][k];
			u1_pri[0][k] = u1_pri[5][k];
			
			u2_pri[nCells + 5][k] = u2_pri[nCells][k];
			u2_pri[nCells + 4][k] = u2_pri[nCells + 1][k];
			u2_pri[nCells + 3][k] = u2_pri[nCells + 2][k];
		}
		
		phi[0] = phi[3];
		phi[1] = phi[3];
		phi[2] = phi[3];
		phi[nCells + 3] = phi[nCells + 2];
		phi[nCells + 4] = phi[nCells + 2];
		phi[nCells + 5] = phi[nCells + 2];
		
		for(int i = 0; i < interface_index + 4; i++){
			u1_con[i] = pri2con(u1_pri[i], gama1, p_inf1);
		}
		
		for(int i = interface_index - 2; i < nCells + 6; i++){
			u2_con[i] = pri2con(u2_pri[i], gama2, p_inf2);
		}
		
		double a_max = -1e8;
		
		for(int i = 0; i < interface_index + 4; i++){
			double cs1 = calcs(u1_pri[i], gama1, p_inf1);
			if(a_max < fabs(u1_pri[i][1]) + cs1){
				a_max = fabs(u1_pri[i][1]) + cs1;
			}
		}
		
		for(int i = interface_index - 2; i < nCells + 6; i++){
			double cs2 = calcs(u2_pri[i], gama2, p_inf2);
			if(a_max < fabs(u2_pri[i][1]) + cs2){
				a_max = fabs(u2_pri[i][1]) + cs2;
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
		
		for(int i = 3; i < nCells + 3; i++){
			if(phi_bar[i] * phi_bar[i + 1] <= 0){
				for(int j = 3; j < i + 1; j++){
					phi_bar[j] = (x0 + dx * (j - 2)) - (x0 + dx * (i - 2)) + phi_bar[i];
				}
				for(int j = i + 1; j < nCells + 3; j++){
					phi_bar[j] = (x0 + dx * (j - 2)) - (x0 + dx * (i + 1 - 2)) + phi_bar[i + 1];
				}
				break;
			}
		}
		
		phi_bar[0] = phi_bar[3];
		phi_bar[1] = phi_bar[3];
		phi_bar[2] = phi_bar[3];
		phi_bar[nCells + 3] = phi_bar[nCells + 2];
		phi_bar[nCells + 4] = phi_bar[nCells + 2];
		phi_bar[nCells + 5] = phi_bar[nCells + 2];
		
		phi = phi_bar;
		
		for(int i = 1; i < interface_index + 3; i++){
			for(int k = 0; k < 3; k++){
				delta1x0[i][k] = u1_con[i][k] - u1_con[i - 1][k];
				delta1x1[i][k] = u1_con[i + 1][k] - u1_con[i][k];
			}
		}
		
		for(int i = interface_index - 1; i < nCells + 5; i++){
			for(int k = 0; k < 3; k++){
				delta2x0[i][k] = u2_con[i][k] - u2_con[i - 1][k];
				delta2x1[i][k] = u2_con[i + 1][k] - u2_con[i][k];
			}
		}
		
		for(int i = 1; i < interface_index + 3; i++){
			for(int k = 0; k < 3; k++){
				delta1x[i][k] = 0.5 * (1 + w1x[k]) * delta1x0[i][k] + 0.5 * (1 - w1x[k]) * delta1x1[i][k];
			}
		}
		
		for(int i = interface_index - 1; i < nCells + 5; i++){
			for(int k = 0; k < 3; k++){
				delta2x[i][k] = 0.5 * (1 + w2x[k]) * delta2x0[i][k] + 0.5 * (1 - w2x[k]) * delta2x1[i][k];
			}
		}
		
		std::vector<std::vector<double>> x1L;
		x1L.resize(nCells + 6,vector<double>(3));
		
		std::vector<std::vector<double>> x1R;
		x1R.resize(nCells + 6,vector<double>(3));
		
		std::vector<std::vector<double>> x2L;
		x2L.resize(nCells + 6,vector<double>(3));
		
		std::vector<std::vector<double>> x2R;
		x2R.resize(nCells + 6,vector<double>(3));
		
		for(int i = 1; i < interface_index + 3; i++){
			std::vector<double> r1x;
			r1x.resize(3);
			for(int k = 0; k < 3; k++){
				r1x[k] = (u1_con[i][k] - u1_con[i - 1][k])/(u1_con[i + 1][k] - u1_con[i][k]);
				if(u1_con[i + 1][k] - u1_con[i][k] == 0){
					r1x[k] = 0;
				}
				x1L[i][k] = u1_con[i][k] - 0.5 * Minbee(r1x[k]) * delta1x[i][k];
				x1R[i][k] = u1_con[i][k] + 0.5 * Minbee(r1x[k]) * delta1x[i][k];
			}
		}
		
		for(int i = interface_index - 1; i < nCells + 5; i++){
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
		std::vector<std::vector<double>> u_bar;
		u_bar.resize(nCells + 6,vector<double>(3));
		
		std::vector<std::vector<double>> u1_bar;
		u1_bar.resize(nCells + 6,vector<double>(3));
		
		std::vector<std::vector<double>> half1xL;
		half1xL.resize(nCells + 6,vector<double>(3));
		
		std::vector<std::vector<double>> half1xR;
		half1xR.resize(nCells + 6,vector<double>(3));
		
		std::vector<std::vector<double>> u2_bar;
		u2_bar.resize(nCells + 6,vector<double>(3));
		
		std::vector<std::vector<double>> half2xL;
		half2xL.resize(nCells + 6,vector<double>(3));
		
		std::vector<std::vector<double>> half2xR;
		half2xR.resize(nCells + 6,vector<double>(3));
		
		for(int i = 1; i < interface_index + 3; i++){
			std::vector<double> f1xL = f(x1L[i], gama1, p_inf1);
			std::vector<double> f1xR = f(x1R[i], gama1, p_inf1);
			for(int k = 0; k < 3; k++){
				half1xL[i][k] = x1L[i][k] - 0.5 * (dt/dx) * (f1xR[k] - f1xL[k]);
				half1xR[i][k] = x1R[i][k] - 0.5 * (dt/dx) * (f1xR[k] - f1xL[k]);
			}
		}
		
		for(int i = interface_index - 1; i < nCells + 5; i++){
			std::vector<double> f2xL = f(x2L[i], gama2, p_inf2);
			std::vector<double> f2xR = f(x2R[i], gama2, p_inf2);
			for(int k = 0; k < 3; k++){
				half2xL[i][k] = x2L[i][k] - 0.5 * (dt/dx) * (f2xR[k] - f2xL[k]);
				half2xR[i][k] = x2R[i][k] - 0.5 * (dt/dx) * (f2xR[k] - f2xL[k]);
			}
		}
		
		std::vector<std::vector<double>> flux1_x;
		flux1_x.resize(nCells + 5,vector<double>(3));
		
		std::vector<std::vector<double>> flux2_x;
		flux2_x.resize(nCells + 5,vector<double>(3));
		
		for(int i = 1; i < interface_index + 3; i++){
			flux1_x[i] = FORCE_flux(half1xR[i], half1xL[i + 1], dt, dx, gama1, p_inf1);
		}
		
		for(int i = interface_index - 1; i < nCells + 5; i++){
			flux2_x[i] = FORCE_flux(half2xR[i], half2xL[i + 1], dt, dx, gama2, p_inf2);
		}
		
		
		for(int i = 2; i < interface_index + 3; i++){
			for(int k = 0; k < 3; k++){
				u1_bar[i][k] = u1_con[i][k] - (dt/dx) * (flux1_x[i][k] - flux1_x[i - 1][k]);
			}
		}
		
		for(int i = interface_index; i < nCells + 5; i++){
			for(int k = 0; k < 3; k++){
				u2_bar[i][k] = u2_con[i][k] - (dt/dx) * (flux2_x[i][k] - flux2_x[i - 1][k]);
			}
		}
		
		u1_con = u1_bar;
		u2_con = u2_bar;
		
		for(int i = 2; i < interface_index + 3; i++){
			u1_pri[i] = con2pri(u1_con[i], gama1, p_inf1);
		}
		
		for(int i = interface_index; i < nCells + 5; i++){
			u2_pri[i] = con2pri(u2_con[i], gama2, p_inf2);
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
		outFile << x0 + dx * (i - 2) << " " << u_pri[i][0] << endl;
	}
	outFile.close();
}