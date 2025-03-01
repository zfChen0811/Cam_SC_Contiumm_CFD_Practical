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
std::vector<double> Fri_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama, char dir);
std::vector<double> RI_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama, char dir);
std::vector<double> FORCE_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama, char dir);
double calecon(std::vector<double> con);
double calepri(std::vector<double> pri, double gama);
double calcs(std::vector<double> pri, double gama);

std::vector<double> fr(std::vector<double> con, double gama){
	std::vector<double> pri = con2pri(con, gama);
	std::vector<double> value;
	value.resize(4);
	value[0] = con[1];
	value[1] = con[1] * pri[1];
	value[2] = con[1] * pri[2];
	value[3] = (con[3] + pri[3]) * pri[1];
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
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double E = con[3];
	std::vector<double> pri = con2pri(con, gama);
	double p = pri[3];
	std::vector<double> f;
	f.resize(4);
	f[0] = rho * vx;
	f[1] = rho * vx * vx + p;
	f[2] = rho * vx * vy;
	f[3] = (E + p) * vx;
	return f;
}

std::vector<double> g(std::vector<double>con, double gama){
	// calculate f
	// based on conservative value roh rohv E
	std::vector<double> g;
	g.resize(4);
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double E = con[3];
	std::vector<double> pri = con2pri(con, gama);
	double p = pri[3];
	g[0] = rho * vy;
	g[1] = rho * vx * vy;
	g[2] = rho * vy * vy + p;
	g[3] = (E + p) * vy;
	return g;
}

std::vector<double> con2pri(std::vector<double>con, double gama){
	// convert conservative (rho rhov E) to primitive (rho v p)
	std::vector<double> pri;
	pri.resize(4);
	// rho v p
	double e = calecon(con);
	pri[0] = con[0];
	pri[1] = con[1] / con[0];
	pri[2] = con[2] / con[0];
	pri[3] = (gama - 1) * con[0] * e;
	return pri;
}

std::vector<double> pri2con(std::vector<double>pri, double gama){
	// convert primitive (rho v p) to conservative  (rho rhov E)
	std::vector<double> con;
	con.resize(4);
	double e = calepri(pri, gama);
	// rho rhov E
	con[0] = pri[0];
	con[1] = pri[0] * pri[1];
	con[2] = pri[0] * pri[2];
	con[3] = pri[0] * e + 0.5 * pri[0] * (pri[1] * pri[1] + pri[2] * pri[2]);
	return con;
}

std::vector<double> Fri_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama, char dir){
	// calculate the flux for lax-friedrich
	std::vector<double> Fri;
	Fri.resize(4);
	
	if(dir == 'x'){
		std::vector<double> f0 = f(con0, gama);
		std::vector<double> f1 = f(con1, gama);
		for(int i = 0; i < 4; i++){
			Fri[i] = 0.5 * (dx/dt) * (con0[i] - con1[i]) + 0.5 * (f1[i] + f0[i]);
		}
	}else if(dir == 'y'){
		std::vector<double> g0 = g(con0, gama);
		std::vector<double> g1 = g(con1, gama);
		for(int i = 0; i < 4; i++){
			Fri[i] = 0.5 * (dx/dt) * (con0[i] - con1[i]) + 0.5 * (g1[i] + g0[i]);
		}
	}else{
		std::cout << "Error in Lax-Friendrich Flux" << std::endl;
	}
	return Fri;
}

std::vector<double> RI_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama, char dir){
	std::vector<double> RI;
	std::vector<double> u;
	u.resize(4);
	
	if(dir == 'x'){
		std::vector<double> f0 = f(con0, gama);
		std::vector<double> f1 = f(con1, gama);
		for(int i = 0; i < 4; i++){
			u[i] = 0.5 * (con0[i] + con1[i]) - 0.5 * (dt/dx) * (f1[i] - f0[i]);
		}
		RI = f(u, gama);
	}else if(dir == 'y'){
		std::vector<double> g0 = g(con0, gama);
		std::vector<double> g1 = g(con1, gama);
		for(int i = 0; i < 4; i++){
			u[i] = 0.5 * (con0[i] + con1[i]) - 0.5 * (dt/dx) * (g1[i] - g0[i]);
		}
		RI = g(u, gama);
	}else{
		std::cout << "Error in RI Flux" << std::endl;
	}
	return RI;
}

std::vector<double> FORCE_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama, char dir){
	std::vector<double> flux;
	flux.resize(4);
	std::vector<double> RI = RI_flux(con0, con1, dt, dx, gama, dir);
	std::vector<double> Fri = Fri_flux(con0, con1, dt, dx, gama, dir);
	for(int i = 0; i < 4 ; i++){
		flux[i] = 0.5 * (Fri[i] + RI[i]);
	}
	return flux;
}

double calecon(std::vector<double> con){
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double E = con[3];
	double e = (E - 0.5 * rho * (vx * vx + vy * vy))/rho;
	return e;
}

double calepri(std::vector<double> pri, double gama){
	double e = pri[3]/((gama - 1) * pri[0]);
	return e;
}

double calcs(std::vector<double> pri, double gama){
	double cs = sqrt(gama *  pri[3]/ pri[0]);
	return cs;
}

using namespace std;
int main(int argc, char *argv[]) {
	double x0 = 0;
	double x1 = 1;
	double y0 = -1;
	double y1 = 1;
	double gama = 1.4;
	double c = 0.9;
	int nCells = 200;
	double nxPoints = nCells + 1;
	double nyPoints = nCells + 1;
	double time = 0.25;
	double t = 0;
	double dx = (x1 - x0)/(nxPoints-1);
	double dy = (y1 - y0)/(nyPoints-1);
	
	std::vector<std::vector< std::vector<double>>> u_pri;
	u_pri.resize(nCells + 3,vector<vector<double>>(nCells + 3,vector<double>(4)));

	std::vector<std::vector< std::vector<double>>> u_con;
	u_con.resize(nCells + 3,vector<vector<double>>(nCells + 3,vector<double>(4)));
	
	std::vector<double> wx;
	wx.resize(4,0);
	
	std::vector<double> wy;
	wy.resize(4,0);
	
	std::vector<std::vector< std::vector<double>>> deltax0;
	deltax0.resize(nCells + 3,vector<vector<double>>(nCells + 3, vector<double>(4)));
	
	std::vector<std::vector< std::vector<double>>> deltax1;
	deltax1.resize(nCells + 3,vector<vector<double>>(nCells + 3, vector<double>(4)));
	
	std::vector<std::vector< std::vector<double>>> deltay0;
	deltay0.resize(nCells + 3,vector<vector<double>>(nCells + 3, vector<double>(4)));
	
	std::vector<std::vector< std::vector<double>>> deltay1;
	deltay1.resize(nCells + 3,vector<vector<double>>(nCells + 3, vector<double>(4)));
	
	std::vector<std::vector< std::vector<double>>> deltax;
	deltax.resize(nCells + 3,vector<vector<double>>(nCells + 3, vector<double>(4)));
	
	std::vector<std::vector< std::vector<double>>> deltay;
	deltay.resize(nCells + 3,vector<vector<double>>(nCells + 3, vector<double>(4)));
		
	for(int i = 2; i < nCells + 2; i++){
		for(int j = 2; j < nCells + 2; j++){
			double x = x0 + (i - 1) * dx;
			double y = y0 + (j - 1) * dy;
			if(x * x + y * y <= 0.4 * 0.4){
				u_pri[i][j][0] = 1;
				u_pri[i][j][1] = 0;
				u_pri[i][j][2] = 0;
				u_pri[i][j][3] = 1;
			}else if(x * x + y * y > 0.4 * 0.4){
				u_pri[i][j][0] = 0.125;
				u_pri[i][j][1] = 0;
				u_pri[i][j][2] = 0;
				u_pri[i][j][3] = 0.1;
			}
		}
	}

	do{
		for(int i = 0; i < nCells + 3; i++){
			u_pri[i][0][0] = u_pri[i][3][0];
			u_pri[i][1][0] = u_pri[i][2][0];
			u_pri[i][0][1] = -u_pri[i][3][1];
			u_pri[i][1][1] = -u_pri[i][2][1];
			u_pri[i][0][2] = u_pri[i][3][2];
			u_pri[i][1][2] = u_pri[i][2][2];
			u_pri[i][0][3] = u_pri[i][3][3];
			u_pri[i][1][3] = u_pri[i][2][3];
			u_pri[i][nCells + 2] = u_pri[i][nCells + 1];
		}
		
		for(int i = 0; i < nCells + 3; i++){
			u_pri[0][i][0] = u_pri[3][i][0];
			u_pri[1][i][0] = u_pri[2][i][0];
			u_pri[0][i][1] = -u_pri[3][i][1];
			u_pri[1][i][1] = -u_pri[2][i][1];
			u_pri[0][i][2] = u_pri[3][i][2];
			u_pri[1][i][2] = u_pri[2][i][2];
			u_pri[0][i][3] = u_pri[3][i][3];
			u_pri[1][i][3] = u_pri[2][i][3];
			u_pri[nCells + 2][i] = u_pri[nCells + 1][i];
		}

		for(int i = 0; i < nCells + 3; i++){
			for(int j = 0; j < nCells + 3; j++){
				u_con[i][j] = pri2con(u_pri[i][j], gama);
			}
		}
		
		double a_max;
		for(int i = 0; i < nCells + 3; i++){
			for(int j = 0; j < nCells + 3; j++){
				double v_ij = sqrt(u_pri[i][j][1] * u_pri[i][j][1] + u_pri[i][j][2] * u_pri[i][j][2]);
				if(i == 0 && j == 0){
					a_max = v_ij + calcs(u_pri[i][j], gama);
				}else if(a_max < v_ij + calcs(u_pri[i][j], gama)){
					a_max = v_ij + calcs(u_pri[i][j], gama);
				}
			}
		}
		
		double dt = c * dx / a_max;
		t = t + dt;
		
		for(int i = 2; i < nCells + 2; i++){
			for(int j = 2; j < nCells + 2; j++){
				for(int k = 0; k < 4; k++){
					deltax0[i][j][k] = u_con[i][j][k] - u_con[i][j - 1][k];
					deltax1[i][j][k] = u_con[i][j + 1][k] - u_con[i][j][k];
				}
			}
		}
		
		for(int i = 2; i < nCells + 2; i++){
			for(int j = 2; j < nCells + 2; j++){
				for(int k = 0; k < 4; k++){
					deltax[i][j][k] = 0.5 * (1 + wx[k]) * deltax0[i][j][k] + 0.5 * (1 - wx[k]) * deltax1[i][j][k];
				}
			}
		}
		
		for(int i = 2; i < nCells + 2; i++){
			deltax[i][0] = deltax[i][1];
			deltax[i][nCells + 1] = deltax[i][nCells];
		}
		
		std::vector<std::vector< std::vector<double>>> xL;
		xL.resize(nCells + 3,vector<vector<double>>(nCells + 3,vector<double>(4)));
		
		std::vector<std::vector< std::vector<double>>> xR;
		xR.resize(nCells + 3,vector<vector<double>>(nCells + 3,vector<double>(4)));
			
		for(int i = 2; i < nCells + 2; i++){
			for(int j = 2; j < nCells + 2; j++){
				std::vector<double> rx;
				rx.resize(4);
				for(int k = 0; k < 4; k++){
					rx[k] = (u_con[i][j][k] - u_con[i][j - 1][k])/(u_con[i][j + 1][k] - u_con[i][j][k]);
					if(u_con[i][j + 1][k] - u_con[i][j][k] == 0){
						rx[k] = 0;
					}
					xL[i][j][k] = u_con[i][j][k] - 0.5 * Minbee(rx[k]) * deltax[i][j][k];
					xR[i][j][k] = u_con[i][j][k] + 0.5 * Minbee(rx[k]) * deltax[i][j][k];
				}
			}
		}
		
		for(int i = 2; i < nCells + 2; i++){
			xL[i][0][0] = xL[i][3][0];
			xL[i][1][0] = xL[i][2][0];
			xL[i][0][1] = -xL[i][3][1];
			xL[i][1][1] = -xL[i][2][1];
			xL[i][0][2] = xL[i][3][2];
			xL[i][1][2] = xL[i][2][2];
			xL[i][0][3] = xL[i][3][3];
			xL[i][1][3] = xL[i][2][3];
			xL[i][nCells + 2] = xL[i][nCells + 1];
			
			xR[i][0][0] = xR[i][3][0];
			xR[i][1][0] = xR[i][2][0];
			xR[i][0][1] = -xR[i][3][1];
			xR[i][1][1] = -xR[i][2][1];
			xR[i][0][2] = xR[i][3][2];
			xR[i][1][2] = xR[i][2][2];
			xR[i][0][3] = xR[i][3][3];
			xR[i][1][3] = xR[i][2][3];
			xR[i][nCells + 2] = xR[i][nCells + 1];
		}
		
		std::vector<std::vector< std::vector<double>>> u_bar1;
		u_bar1.resize(nCells + 3,vector<vector<double>>(nCells + 3,vector<double>(4)));
		
		std::vector<std::vector< std::vector<double>>> halfxL;
		halfxL.resize(nCells + 3,vector<vector<double>>(nCells + 3,vector<double>(4)));
		
		std::vector<std::vector< std::vector<double>>> halfxR;
		halfxR.resize(nCells + 3,vector<vector<double>>(nCells + 3,vector<double>(4)));
		
		for(int i = 2; i < nCells + 2; i++){
			for(int j = 2; j < nCells + 2; j++){
				std::vector<double> fxL = f(xL[i][j], gama);
				std::vector<double> fxR = f(xR[i][j], gama);
				
				for(int k = 0; k < 4; k++){
					halfxL[i][j][k] = xL[i][j][k] - 0.5 * (dt/dx) * (fxR[k] - fxL[k]);
					halfxR[i][j][k] = xR[i][j][k] - 0.5 * (dt/dx) * (fxR[k] - fxL[k]);
				}
			}
		}
		
		for(int i = 2; i < nCells + 2; i++){
			halfxL[i][0][0] = halfxL[i][3][0];
			halfxL[i][1][0] = halfxL[i][2][0];
			halfxL[i][0][1] = -halfxL[i][3][1];
			halfxL[i][1][1] = -halfxL[i][2][1];
			halfxL[i][0][2] = halfxL[i][3][2];
			halfxL[i][1][2] = halfxL[i][2][2];
			halfxL[i][0][3] = halfxL[i][3][3];
			halfxL[i][1][3] = halfxL[i][2][3];
			halfxL[i][nCells + 2] = halfxL[i][nCells + 1];
			
			halfxR[i][0][0] = halfxR[i][3][0];
			halfxR[i][1][0] = halfxR[i][2][0];
			halfxR[i][0][1] = -halfxR[i][3][1];
			halfxR[i][1][1] = -halfxR[i][2][1];
			halfxR[i][0][2] = halfxR[i][3][2];
			halfxR[i][1][2] = halfxR[i][2][2];
			halfxR[i][0][3] = halfxR[i][3][3];
			halfxR[i][1][3] = halfxR[i][2][3];
			halfxR[i][nCells + 2] = halfxR[i][nCells + 1];
		}
		
		std::vector<std::vector< std::vector<double>>> flux_x;
		flux_x.resize(nCells + 3,vector<vector<double>>(nCells + 2,vector<double>(4)));
		
		for(int i = 2; i < nCells + 2; i++){
			for(int j = 0; j < nCells + 2; j++){
				flux_x[i][j] = FORCE_flux(halfxR[i][j], halfxL[i][j + 1], dt, dx, gama, 'x');
			}
		}

		for(int i = 2; i < nCells + 2; i++){
			for(int j = 2; j < nCells + 2; j++){
				for(int k = 0; k < 4; k++){
					u_bar1[i][j][k] = u_con[i][j][k] - (dt/dx) * (flux_x[i][j][k] - flux_x[i][j - 1][k]);
				}
			}
		}
		
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		for(int i = 0; i < nCells + 3; i++){
			u_bar1[i][0][0] = u_bar1[i][3][0];
			u_bar1[i][1][0] = u_bar1[i][2][0];
			u_bar1[i][0][1] = -u_bar1[i][3][1];
			u_bar1[i][1][1] = -u_bar1[i][2][1];
			u_bar1[i][0][2] = u_bar1[i][3][2];
			u_bar1[i][1][2] = u_bar1[i][2][2];
			u_bar1[i][0][3] = u_bar1[i][3][3];
			u_bar1[i][1][3] = u_bar1[i][2][3];
			u_bar1[i][nCells + 2] = u_bar1[i][nCells + 1];
		}
		
		for(int i = 0; i < nCells + 3; i++){
			u_bar1[0][i][0] = u_bar1[3][i][0];
			u_bar1[1][i][0] = u_bar1[2][i][0];
			u_bar1[0][i][1] = -u_bar1[3][i][1];
			u_bar1[1][i][1] = -u_bar1[2][i][1];
			u_bar1[0][i][2] = u_bar1[3][i][2];
			u_bar1[1][i][2] = u_bar1[2][i][2];
			u_bar1[0][i][3] = u_bar1[3][i][3];
			u_bar1[1][i][3] = u_bar1[2][i][3];
			u_bar1[nCells + 2][i] = u_bar1[nCells + 1][i];
		}
		
		for(int i = 2; i < nCells + 2; i++){
			for(int j = 2; j < nCells + 2; j++){
				for(int k = 0; k < 4; k++){
					deltay0[i][j][k] = u_bar1[i][j][k] - u_bar1[i - 1][j][k];
					deltay1[i][j][k] = u_bar1[i + 1][j][k] - u_bar1[i][j][k];
				}
			}
		}
		
		for(int i = 2; i < nCells + 2; i++){
			deltay0[0][i] = deltay0[1][i];
			deltay0[nCells + 1][i] = deltay0[nCells][i];
			deltay1[0][i] = deltay1[1][i];
			deltay1[nCells + 1][i] = deltay1[nCells][i];
		}
		
		for(int i = 2; i < nCells + 2; i++){
			for(int j = 2; j < nCells + 2; j++){
				for(int k = 0; k < 4; k++){
					deltay[i][j][k] = 0.5 * (1 + wy[k]) * deltay0[i][j][k] + 0.5 * (1 - wy[k]) * deltay1[i][j][k];
				}
			}
		}
		
		for(int i = 2; i < nCells + 2; i++){
			deltay[0][i] = deltay[1][i];
			deltay[nCells + 1][i] = deltay[nCells][i];
		}
		
		std::vector<std::vector< std::vector<double>>> yL;
		yL.resize(nCells + 3,vector<vector<double>>(nCells + 3,vector<double>(4)));
		
		std::vector<std::vector< std::vector<double>>> yR;
		yR.resize(nCells + 3,vector<vector<double>>(nCells + 3,vector<double>(4)));

		for(int i = 2; i < nCells + 2; i++){
			for(int j = 2; j < nCells + 2; j++){
				std::vector<double> ry;
				ry.resize(4);
				for(int k = 0; k < 4; k++){
					ry[k] = (u_bar1[i][j][k] - u_bar1[i - 1][j][k])/(u_bar1[i + 1][j][k] - u_bar1[i][j][k]);
					if(u_bar1[i + 1][j][k] - u_bar1[i][j][k] == 0){
						ry[k] = 0;
					}
					yL[i][j][k] = u_bar1[i][j][k] - 0.5 * Minbee(ry[k]) * deltay[i][j][k];
					yR[i][j][k] = u_bar1[i][j][k] + 0.5 * Minbee(ry[k]) * deltay[i][j][k];
				}
			}
		}
	
		for(int i = 2; i < nCells + 2; i++){
			yL[0][i][0] = yL[3][i][0];
			yL[1][i][0] = yL[2][i][0];
			yL[0][i][1] = -yL[3][i][1];
			yL[1][i][1] = -yL[2][i][1];
			yL[0][i][2] = yL[3][i][2];
			yL[1][i][2] = yL[2][i][2];
			yL[0][i][3] = yL[3][i][3];
			yL[1][i][3] = yL[2][i][3];
			yL[nCells + 2][i] = yL[nCells + 1][i];
			
			yR[0][i][0] = yR[3][i][0];
			yR[1][i][0] = yR[2][i][0];
			yR[0][i][1] = -yR[3][i][1];
			yR[1][i][1] = -yR[2][i][1];
			yR[0][i][2] = yR[3][i][2];
			yR[1][i][2] = yR[2][i][2];
			yR[0][i][3] = yR[3][i][3];
			yR[1][i][3] = yR[2][i][3];
			yR[nCells + 2][i] = yR[nCells + 1][i];
		}
				
		std::vector<std::vector< std::vector<double>>> u_bar2;
		u_bar2.resize(nCells + 3,vector<vector<double> >(nCells + 3,vector<double>(4)));
		
		std::vector<std::vector< std::vector<double>>> halfyL;
		halfyL.resize(nCells + 3,vector<vector<double>>(nCells + 3,vector<double>(4)));
		
		std::vector<std::vector< std::vector<double>>> halfyR;
		halfyR.resize(nCells + 3,vector<vector<double>>(nCells + 3,vector<double>(4)));

		for(int i = 2; i < nCells + 2; i++){
			for(int j = 2; j < nCells + 2; j++){
				std::vector<double> gyL = g(yL[i][j], gama);
				std::vector<double> gyR = g(yR[i][j], gama);
				
				for(int k = 0; k < 4; k++){
					halfyL[i][j][k] = yL[i][j][k] - 0.5 * (dt/dy) * (gyR[k] - gyL[k]);
					halfyR[i][j][k] = yR[i][j][k] - 0.5 * (dt/dy) * (gyR[k] - gyL[k]);
				}
			}
		}

		for(int i = 2; i < nCells + 2; i++){
			halfyL[0][i][0] = halfyL[3][i][0];
			halfyL[1][i][0] = halfyL[2][i][0];
			halfyL[0][i][1] = -halfyL[3][i][1];
			halfyL[1][i][1] = -halfyL[2][i][1];
			halfyL[0][i][2] = halfyL[3][i][2];
			halfyL[1][i][2] = halfyL[2][i][2];
			halfyL[0][i][3] = halfyL[3][i][3];
			halfyL[1][i][3] = halfyL[2][i][3];
			halfyL[nCells + 2][i] = halfyL[nCells + 1][i];
			
			halfyR[0][i][0] = halfyR[3][i][0];
			halfyR[1][i][0] = halfyR[2][i][0];
			halfyR[0][i][1] = -halfyR[3][i][1];
			halfyR[1][i][1] = -halfyR[2][i][1];
			halfyR[0][i][2] = halfyR[3][i][2];
			halfyR[1][i][2] = halfyR[2][i][2];
			halfyR[0][i][3] = halfyR[3][i][3];
			halfyR[1][i][3] = halfyR[2][i][3];
			halfyR[nCells + 2][i] = halfyR[nCells + 1][i];
		}
		
		std::vector<std::vector< std::vector<double>>> flux_y;
		flux_y.resize(nCells + 2,vector<vector<double>>(nCells + 3,vector<double>(4)));
		
		for(int j = 2; j < nCells + 2; j++){
			for(int i = 0; i < nCells + 2; i++){
				flux_y[i][j] = FORCE_flux(halfyR[i][j], halfyL[i + 1][j], dt, dy, gama, 'y');
			}
		}

		for(int j = 2; j < nCells + 2; j++){
			for(int i = 2; i < nCells + 2; i++){
				for(int k = 0; k < 4; k++){
					u_bar2[i][j][k] = u_bar1[i][j][k] - (dt/dy) * (flux_y[i][j][k] - flux_y[i - 1][j][k]);
				}
			}
		}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
		std::vector<std::vector<std::vector<double>>> K1;
		std::vector<std::vector<std::vector<double>>> K2;
		K1.resize(nCells + 3, vector<vector<double>>(nCells + 3,vector<double>(4)));
		K2.resize(nCells + 3, vector<vector<double>>(nCells + 3,vector<double>(4)));
		
		std::vector<std::vector< std::vector<double>>> u_con_plus;
		u_con_plus.resize(nCells + 3,vector<vector<double>>(nCells + 3,vector<double>(4)));
		
		for(int i = 2; i < nCells + 2; i++){
			for(int j = 2; j < nCells + 2; j++){
				double x = x0 + dx * (i - 0.5 - 1);
				std::vector<double> u;
				std::vector<double> s;
				s.resize(4);
				
				u = u_bar2[i][j];
				
				std::vector<double> fr_value = fr(u, gama);
				
				for(int k = 0; k < 4; k++){
					s[k] = (-x) * fr_value[k];
					K1[i][j][k] = dt * s[k];
				}
			}
		}
		
		for(int i = 2; i < nCells + 2; i++){
			for(int j = 2; j < nCells + 2; j++){
				double x = x0 + dx * (i - 0.5 - 1);
				std::vector<double> u;
				u.resize(4);
				std::vector<double> s;
				s.resize(4);
				
				for(int k = 0; k < 4; k++){
					u[k] = u_bar2[i][j][k] + 0.5 * K1[i][j][k];
				}
				
				std::vector<double> fr_value = fr(u, gama);
				
				for(int k = 0; k < 4; k++){
					s[k] = (-x) * fr_value[k];
					K2[i][j][k] = dt * s[k];
				}
			}
		}
		
		for(int i = 2; i < nCells + 2; i++){
			for(int j = 2; j < nCells + 2; j++){
				for(int k = 0; k < 4; k++){
					u_con_plus[i][j][k] = u_bar2[i][j][k] + (1.0/2.0) * (K1[i][j][k] + K2[i][j][k]);
				}
			}
		}
		
		u_con = u_con_plus;
		
		for(int i = 2; i < nCells + 2; i++){
			for(int j = 2; j < nCells + 2; j++){
				u_pri[i][j] = con2pri(u_con[i][j], gama);
			}
		}
	}while(t < time);
	
	ofstream outFile("/Users/chenzefeng/Desktop/rho.dat");
	for(int i = 2; i < nCells + 2; i++){
		for(int j = 2; j < nCells + 2; j++){
			outFile << x0 + dx * (i - 0.5 - 1) << " " << y0 + dy * (j - 0.5 - 1) << " " << u_pri[i][j][0] << endl;
		}
	}
	outFile.close();
}