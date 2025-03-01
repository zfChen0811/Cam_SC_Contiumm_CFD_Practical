#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>
#include <fstream>
using namespace std;

double fK(double p_star, std::vector<double>uK, double gama, double p_inf);
double fK_d(double p_star, std::vector<double>uK, double gama, double p_inf);
double f_star(double p_star, std::vector<double>uL, std::vector<double>uR, double gama1, double gama2, double p_inf1, double p_inf2);
double f_star_d(double p_star, std::vector<double>uL, std::vector<double>uR, double gama1, double gama2, double p_inf1, double p_inf2);
std::vector<double> f(std::vector<double>con, double gama, double p_inf);
double Minbee(double r);
std::vector<double> con2pri(std::vector<double>con, double gama, double p_inf);
std::vector<double> pri2con(std::vector<double>pri,double gama, double p_inf);
double calecon(std::vector<double> con);
double calepri(std::vector<double> pri, double gama, double p_inf);
double calcs(std::vector<double> pri, double gama, double p_inf);
double caldt(std::vector<std::vector<double>> u_pri, double gama, double p_inf, double c, double dx, double start_index, double end_index);
std::vector<double> Fri_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama, double p_inf);
std::vector<double> RI_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama, double p_inf);
std::vector<double> FORCE_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama, double p_inf);
std::vector<std::vector<double>> exact_solver(std::vector<double> priL, std::vector<double> priR, double gamaL, double gamaR, double p_infL, double p_infR);
std::vector<double> phi_update(std::vector<double> phi, std::vector<std::vector<double>> u1_pri, std::vector<std::vector<double>> u2_pri, std::vector<std::vector<double>> u1_con, std::vector<std::vector<double>> u2_con, int nCells, int ghost_cells, double dt, double dx);
std::vector<double> phi_reinit(std::vector<double> phi, int nCells, int ghost_cells, double x0, double dx);
std::vector<std::vector<double>> trans_bound(std::vector<std::vector<double>> u, int nCells, int ghost_cells);
double find_interface(std::vector<double>phi, int nCells, int ghost_cells);
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> halftime_update(std::vector<std::vector<double>>u_con, double start_index, double end_index, double gama, double p_inf, int nCells, int ghost_cells, double dt, double dx);
std::vector<std::vector<double>> SLIC(std::vector<std::vector<double>>u_con, double start_index, double end_index, double gama, double p_inf, int nCells, int ghost_cells, double dt, double dx);

double fK(double p_star, std::vector<double>uK, double gama, double p_inf){
	std::vector<double> pri = con2pri(uK, gama, p_inf);
	double rhoK = pri[0];
	double pK = pri[2];
	double f_value;
	double AK = 2 / ((gama + 1) * rhoK);
	double BK = ((gama - 1) / (gama + 1)) * pK + (2 * gama * p_inf) / (gama + 1);
	double csK = calcs(pri, gama, p_inf);
	if(p_star > pK){
		f_value = (p_star - pK) * sqrt(AK/(p_star + BK));
	}else{
		f_value = ((2 * csK) / (gama - 1)) * (pow(((p_star + p_inf)/ (pK + p_inf)), (gama - 1)/(2 * gama)) - 1);
	}
	return f_value;
}

double fK_d(double p_star, std::vector<double>uK, double gama, double p_inf){
	std::vector<double> pri = con2pri(uK, gama, p_inf);
	double rhoK = pri[0];
	double pK = pri[2];
	double f_d;
	double AK = 2 / ((gama + 1) * rhoK);
	double BK = ((gama - 1) / (gama + 1)) * pK + (2 * gama * p_inf) / (gama + 1);
	double csK = calcs(pri, gama, p_inf);
	if(p_star > pK){
		f_d = sqrt(AK / (BK + p_star)) * (1 - ((p_star - pK) / (2 * (BK + p_star))));
	}else{
		f_d = (1 / (rhoK * csK)) * pow((p_star + p_inf) / (pK + p_inf), -(gama + 1) / (2 * gama)) ;
	}
	return f_d;
}

double f_star(double p_star, std::vector<double>uL, std::vector<double>uR, double gama1, double gama2, double p_inf1, double p_inf2){
	double vL = uL[1] / uL[0];
	double vR = uR[1] / uR[0];
	double dv = vR - vL;
	double f_value;
	f_value = fK(p_star, uR, gama2, p_inf2) + fK(p_star, uL, gama1, p_inf1) + dv;
	return f_value;
}

double f_star_d(double p_star, std::vector<double>uL, std::vector<double>uR, double gama1, double gama2, double p_inf1, double p_inf2){
	double f_d_value;
	f_d_value = fK_d(p_star, uR, gama2, p_inf2) + fK_d(p_star, uL, gama1, p_inf1);
	return f_d_value;
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

double caldt(std::vector<std::vector<double>> u_pri, double gama, double p_inf, double c, double dx, double start_index, double end_index){
	double a_max = -1e20;
	for(int i = start_index; i < end_index + 1; i++){
		double cs = calcs(u_pri[i], gama, p_inf);
		if(a_max < fabs(u_pri[i][1]) + cs){
			a_max = fabs(u_pri[i][1]) + cs;
		}
	}
	double dt = c * dx / a_max;
	return dt;
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

std::vector<std::vector<double>> exact_solver(std::vector<double> priL, std::vector<double> priR, double gamaL, double gamaR, double p_infL, double p_infR){
	double rhoL = priL[0];
	double rhoR = priR[0];
	double vL = priL[1];
	double vR = priR[1];
	double pL = priL[2];
	double pR = priR[2];
	double csL = calcs(priL, gamaL, p_infL);
	double csR = calcs(priR, gamaR, p_infR);
	double p_new, p_old;
	std::vector<double> u_pri_half;
	u_pri_half.resize(3,0);
	std::vector<double> uL = pri2con(priL, gamaL, p_infL);
	std::vector<double> uR = pri2con(priR, gamaR, p_infR);
	double eplison = 10e-16;
	p_new = pow((csL + csR - 0.5 * (0.5 * (gamaL + gamaR) - 1) * (vR - vL)) / (csL/pow(pL,(gamaL - 1) / (2 * gamaL)) + (csR/pow(pR, (gamaR - 1) / (2 * gamaR)))), (2 * 0.5 * (gamaL + gamaR)) / (0.5 * (gamaL + gamaR) - 1));
	
	do{
		p_old = p_new;
		p_new = p_old - f_star(p_old, uL, uR, gamaL, gamaR, p_infL, p_infR) / f_star_d(p_old, uL, uR, gamaL, gamaR, p_infL, p_infR);
	}while(fabs(p_new - p_old) / fabs(p_old) > eplison);
	
	double p_star = p_new;
	double v_star = 0.5 * (vL + vR) + 0.5 * (fK(p_star, uR, gamaR, p_infR) - fK(p_star, uL, gamaL, p_infL));
	double rhoL_star_shock = rhoL * ((((p_star + p_infL) / (pL+ p_infL)) +  ((gamaL - 1) / (gamaL + 1))) / (((gamaL - 1) / (gamaL + 1)) * ((p_star + p_infL)/ (pL + p_infL)) + 1));
	double rhoL_star_rare = rhoL * pow(((p_star + p_infL)/(pL+ p_infL)), (1/gamaL));
	double rhoR_star_shock = rhoR * ((((p_star + p_infR)/ (pR + p_infR)) +  ((gamaR - 1) / (gamaR + 1))) / (((gamaR - 1) / (gamaR + 1)) * ((p_star + p_infR)/ (pL + p_infR)) + 1));
	double rhoR_star_rare = rhoR * pow(((p_star + p_infR)/(pR + p_infR)), (1/gamaR));
	std::vector<double> priL_star, priR_star;
	priL_star.resize(3,0);
	priR_star.resize(3,0);
	
	if(p_star > pL){
		// left shock
		priL_star[0] = rhoL_star_shock;
		priL_star[1] = v_star;
		priL_star[2] = p_star;
	}else{
		// left rarefaction
		priL_star[0] = rhoL_star_rare;
		priL_star[1] = v_star;
		priL_star[2] = p_star;
	}
	
	if(p_star > pR){
		// right shock
		priR_star[0] = rhoR_star_shock;
		priR_star[1] = v_star;
		priR_star[2] = p_star;
	}else{
		// right rarefaction
		priR_star[0] = rhoR_star_rare;
		priR_star[1] = v_star;
		priR_star[2] = p_star;
	}
	std::vector<std::vector<double>> result;
	result.resize(2,vector<double>(3));
	result[0] = priL_star;
	result[1] = priR_star;
	return result;
}

std::vector<double> phi_update(std::vector<double> phi, std::vector<std::vector<double>> u1_pri, std::vector<std::vector<double>> u2_pri, std::vector<std::vector<double>> u1_con, std::vector<std::vector<double>> u2_con, int nCells, int ghost_cells, double dt, double dx){
	
	std::vector<double> phi_bar;
	phi_bar.resize(nCells + ghost_cells * 2, 0);

	for(int i = ghost_cells; i < nCells + ghost_cells; i++){
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
	return phi_bar;
}
std::vector<double> phi_reinit(std::vector<double> phi, int nCells, int ghost_cells, double x0, double dx){
	for(int i = ghost_cells; i < nCells + ghost_cells - 1; i++){
		if(phi[i] * phi[i + 1] <= 0){
			for(int j = ghost_cells; j < i + 1; j++){
				phi[j] = (x0 + dx * (j - 2)) - (x0 + dx * (i - 2)) + phi[i];
			}
			for(int j = i + 1; j < nCells + ghost_cells; j++){
				phi[j] = (x0 + dx * (j - 2)) - (x0 + dx * (i + 1 - 2)) + phi[i + 1];
			}
			break;
		}
	}
	return phi;
}

std::vector<std::vector<double>> trans_bound(std::vector<std::vector<double>> u, int nCells, int ghost_cells){
	for(int k = 0; k < 3; k++){
		for(int i = 0; i < ghost_cells; i++){
			u[i][k] = u[ghost_cells][k];
		}
		
		for(int i = nCells + 2 * ghost_cells - 1; i > nCells + ghost_cells - 1; i--){
			u[i][k] = u[nCells + ghost_cells - 1][k];
		}
	}
	return u;
}

double find_interface(std::vector<double>phi, int nCells, int ghost_cells){
	double interface_index;
	for(int i = ghost_cells; i < nCells + 2 * ghost_cells - 1; i++){
		if(phi[i] * phi[i + 1] <= 0){
			interface_index = i;
			break;
		}
	}
	return interface_index;
}

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> halftime_update(std::vector<std::vector<double>>u_con, double start_index, double end_index, double gama, double p_inf, int nCells, int ghost_cells, double dt, double dx){
	
	std::vector<double> wx;
	wx.resize(3,0);
	std::vector<std::vector<double>> deltax0, deltax1, deltax, halfxL, halfxR, xL, xR;
	deltax0.resize(nCells + 2 * ghost_cells,vector<double>(3));
	deltax1.resize(nCells + 2 * ghost_cells,vector<double>(3));
	deltax.resize(nCells + 2 * ghost_cells, vector<double>(3));
	halfxL.resize(nCells + 2 * ghost_cells,vector<double>(3));
	halfxR.resize(nCells + 2 * ghost_cells,vector<double>(3));
	xL.resize(nCells + 2 * ghost_cells,vector<double>(3));
	xR.resize(nCells + 2 * ghost_cells,vector<double>(3));
	
	for(int i = start_index + 1; i < (end_index - 1) + 1; i++){
		for(int k = 0; k < 3; k++){
			deltax0[i][k] = u_con[i][k] - u_con[i - 1][k];
			deltax1[i][k] = u_con[i + 1][k] - u_con[i][k];
			deltax[i][k] = 0.5 * (1 + wx[k]) * deltax0[i][k] + 0.5 * (1 - wx[k]) * deltax1[i][k];
		}
	}
	
	for(int i = start_index + 1; i < (end_index - 1) + 1; i++){
		std::vector<double> r1x;
		r1x.resize(3);
		for(int k = 0; k < 3; k++){
			r1x[k] = (u_con[i][k] - u_con[i - 1][k])/(u_con[i + 1][k] - u_con[i][k]);
			if(u_con[i + 1][k] - u_con[i][k] == 0){
				r1x[k] = 0;
			}
			xL[i][k] = u_con[i][k] - 0.5 * Minbee(r1x[k]) * deltax[i][k];
			xR[i][k] = u_con[i][k] + 0.5 * Minbee(r1x[k]) * deltax[i][k];
		}
	}
	
	for(int i = start_index + 1; i < (end_index - 1) + 1; i++){
		std::vector<double> f1xL = f(xL[i], gama, p_inf);
		std::vector<double> f1xR = f(xR[i], gama, p_inf);
		for(int k = 0; k < 3; k++){
			halfxL[i][k] = xL[i][k] - 0.5 * (dt/dx) * (f1xR[k] - f1xL[k]);
			halfxR[i][k] = xR[i][k] - 0.5 * (dt/dx) * (f1xR[k] - f1xL[k]);
		}
	}
	return std::make_pair(halfxL, halfxR);
}

std::vector<std::vector<double>> SLIC(std::vector<std::vector<double>>u_con, double start_index, double end_index, double gama, double p_inf, int nCells, int ghost_cells, double dt, double dx){
	std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> halfx = halftime_update(u_con, start_index,  end_index, gama, p_inf, nCells, ghost_cells, dt, dx);
	std::vector<std::vector<double>> halfxL = halfx.first;
	std::vector<std::vector<double>> halfxR = halfx.second;
	std::vector<std::vector<double>> flux;
	flux.resize(nCells + 2 * ghost_cells - 1,vector<double>(3));
	
	for(int i = start_index + 1; i < (end_index - 1) + 1; i++){
		flux[i] = FORCE_flux(halfxR[i], halfxL[i + 1], dt, dx, gama, p_inf);
	}
	return flux;
}

int main(int argc, char *argv[]) {
	double x0 = 0;
	double x1 = 1;
	double gama1 = 4.4;
	double gama2 = 1.4;
	double p_inf1 = 6 * pow(10, 8);
	double p_inf2 = 0;
	double c = 0.8;
	int nCells = 1000;
	int ghost_cells = 3;
	double nxPoints = nCells + 1;
	double time = 237.44 * pow(10, -6);
	double t = 0;
	double dx = (x1 - x0)/(nxPoints-1);
	std::vector<double> phi, phi_bar;
	phi.resize(nCells + 6,0);
	std::vector<std::vector<double>> u_pri, u1_pri, u1_con, u2_pri, u2_con, flux1, flux2, star_state;
	u_pri.resize(nCells + 6,vector<double>(3));
	u1_pri.resize(nCells + 6,vector<double>(3));
	u1_con.resize(nCells + 6,vector<double>(3));
	u2_pri.resize(nCells + 6,vector<double>(3));
	u2_con.resize(nCells + 6,vector<double>(3));
		
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
		double interface_index = find_interface(phi, nCells, ghost_cells);
			
		for(int i = 3; i < nCells + 3; i++){
			u1_con[i] = pri2con(u1_pri[i], gama1, p_inf1);
			u2_con[i] = pri2con(u2_pri[i], gama2, p_inf2);
		}
		
		u1_pri = trans_bound(u1_pri, nCells, ghost_cells);
		u2_pri = trans_bound(u2_pri, nCells, ghost_cells);
		
		star_state = exact_solver(u1_pri[interface_index], u2_pri[interface_index + 1], gama1, gama2, p_inf1, p_inf2);
		
		for(int i = interface_index + 1; i < interface_index + 4; i++){
			u1_pri[i] = star_state[0];
		}
		
		for(int i = interface_index - 2; i < interface_index + 1; i++){
			u2_pri[i] = star_state[1];
		}
		
		for(int i = 0; i < interface_index + 4; i++){
			u1_con[i] = pri2con(u1_pri[i], gama1, p_inf1);
		}
		
		for(int i = interface_index - 2; i < nCells + 6; i++){
			u2_con[i] = pri2con(u2_pri[i], gama2, p_inf2);
		}
		
		double dt1 = caldt(u1_pri, gama1, p_inf1, c, dx, 0, interface_index + 3);
		double dt2 = caldt(u2_pri, gama2, p_inf2, c, dx, interface_index - 2, nCells + 5);
		double dt = fmin(dt1, dt2);
		t = t + dt;
		
		std::vector<double> phi_bar;
		phi_bar =  phi_update(phi, u1_pri, u2_pri, u1_con, u2_con, nCells, ghost_cells, dt, dx);
		phi = phi_reinit(phi_bar, nCells, ghost_cells, x0, dx);
		
		flux1 = SLIC(u1_con, 0, interface_index + 3, gama1, p_inf1, nCells, ghost_cells, dt, dx);
		flux2 = SLIC(u2_con, interface_index - 2, nCells + 5, gama2, p_inf2, nCells, ghost_cells, dt, dx);
		
		for(int i = 2; i < interface_index + 3; i++){
			for(int k = 0; k < 3; k++){
				u1_con[i][k] = u1_con[i][k] - (dt/dx) * (flux1[i][k] - flux1[i - 1][k]);
			}
		}
		
		for(int i = interface_index; i < nCells + 5; i++){
			for(int k = 0; k < 3; k++){
				u2_con[i][k] = u2_con[i][k] - (dt/dx) * (flux2[i][k] - flux2[i - 1][k]);
			}
		}
		
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