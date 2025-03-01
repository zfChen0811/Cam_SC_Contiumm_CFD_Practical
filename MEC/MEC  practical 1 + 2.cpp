#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>
#include <fstream>

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
	double x1 = 2;
	double y0 = 0;
	double y1 = 2;
	double gama = 1.4;
	double c = 0.8;
	int nCells = 100;
	double nxPoints = nCells + 1;
	double nyPoints = nCells + 1;
	double time = 0.25;
	double t = 0;
	double dx = (x1 - x0)/(nxPoints-1);
	double dy = (y1 - y0)/(nyPoints-1);
	
	std::vector<std::vector< std::vector<double>>> u_pri;
	u_pri.resize(nCells + 2,vector<vector<double>>(nCells + 2,vector<double>(4)));

	std::vector<std::vector< std::vector<double>>> u_con;
	u_con.resize(nCells + 2,vector<vector<double>>(nCells + 2,vector<double>(4)));
	
	std::vector<double> wx;
	wx.resize(4,0);
	
	std::vector<double> wy;
	wy.resize(4,0);
	
	std::vector<std::vector< std::vector<double>>> deltax0;
	deltax0.resize(nCells + 2,vector<vector<double>>(nCells + 2,vector<double>(4)));
	
	std::vector<std::vector< std::vector<double>>> deltax1;
	deltax1.resize(nCells + 2,vector<vector<double>>(nCells + 2,vector<double>(4)));
	
	std::vector<std::vector< std::vector<double>>> deltay0;
	deltay0.resize(nCells + 2,vector<vector<double>>(nCells + 2,vector<double>(4)));
	
	std::vector<std::vector< std::vector<double>>> deltay1;
	deltay1.resize(nCells + 2,vector<vector<double>>(nCells + 2,vector<double>(4)));
	
	std::vector<std::vector< std::vector<double>>> deltax;
	deltax.resize(nCells + 2,vector<vector<double>>(nCells + 2,vector<double>(4)));
	
	std::vector<std::vector< std::vector<double>>> deltay;
	deltay.resize(nCells + 2,vector<vector<double>>(nCells + 2,vector<double>(4)));
		
	for(int i = 1; i < nCells + 1; i++){
		for(int j = 1; j < nCells + 1; j++){
			double x = x0 + (i - 1) * dx;
			double y = y0 + (j - 1) * dy;
			if((x - 1) * (x - 1) + (y - 1) * (y - 1) <= 0.4 * 0.4){
				u_pri[i][j][0] = 1;
				u_pri[i][j][1] = 0;
				u_pri[i][j][2] = 0;
				u_pri[i][j][3] = 1;
			}else if((x - 1) * (x - 1) + (y - 1) * (y - 1) > 0.4 * 0.4){
				u_pri[i][j][0] = 0.125;
				u_pri[i][j][1] = 0;
				u_pri[i][j][2] = 0;
				u_pri[i][j][3] = 0.1;
			}
//			if(x >= 0.5 && y >= 0.5){
//				u_pri[i][j][0] = 1.5;
//				u_pri[i][j][1] = 0;
//				u_pri[i][j][2] = 0;
//				u_pri[i][j][3] = 1.5;
//			}else if(x >= 0.5 && y < 0.5){
//				u_pri[i][j][0] = 0.5323;
//				u_pri[i][j][1] = 0;
//				u_pri[i][j][2] = 1.206;
//				u_pri[i][j][3] = 1.3;
//			}else if(x < 0.5 && y>= 0.5){
//				u_pri[i][j][0] = 0.5323;
//				u_pri[i][j][1] = 1.206;
//				u_pri[i][j][2] = 0;
//				u_pri[i][j][3] = 0.3;
//			}else{
//				u_pri[i][j][0] = 0.138;
//				u_pri[i][j][1] = 1.206;
//				u_pri[i][j][2] = 1.206;
//				u_pri[i][j][3] = 0.029;
//			}
		}
	}

	do{
		for(int i = 0; i < nCells + 2; i++){
			u_pri[i][0] = u_pri[i][1];
			u_pri[i][nCells + 1] = u_pri[i][nCells];
		}
		
		for(int i = 0; i < nCells + 2; i++){
			u_pri[0][i] = u_pri[1][i];
			u_pri[nCells + 1][i] = u_pri[nCells][i];
		}

		for(int i = 0; i < nCells + 2; i++){
			for(int j = 0; j < nCells + 2; j++){
				u_con[i][j] = pri2con(u_pri[i][j], gama);
			}
		}
		
		double a_max;
		for(int i = 0; i < nCells + 2; i++){
			for(int j = 0; j < nCells + 2; j++){
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
		
		for(int i = 1; i < nCells + 1; i++){
			for(int j = 1; j < nCells + 1; j++){
				for(int k = 0; k < 4; k++){
					deltax0[i][j][k] = u_con[i][j][k] - u_con[i - 1][j][k];
					deltax1[i][j][k] = u_con[i + 1][j][k] - u_con[i][j][k];
				}
			}
		}
		
		for(int i = 1; i < nCells + 1; i++){
			for(int j = 1; j < nCells + 1; j++){
				for(int k = 0; k < 4; k++){
					deltax[i][j][k] = 0.5 * (1 + wx[k]) * deltax0[i][j][k] + 0.5 * (1 - wx[k]) * deltax1[i][j][k];
				}
			}
		}
		
		for(int i = 1; i < nCells + 1; i++){
			deltax[i][0] = deltax[i][1];
			deltax[i][nCells + 1] = deltax[i][nCells];
		}
		
		std::vector<std::vector< std::vector<double>>> xL;
		xL.resize(nCells + 2,vector<vector<double>>(nCells + 2,vector<double>(4)));
		
		std::vector<std::vector< std::vector<double>>> xR;
		xR.resize(nCells + 2,vector<vector<double>>(nCells + 2,vector<double>(4)));
			
		for(int i = 1; i < nCells + 1; i++){
			for(int j = 1; j < nCells + 1; j++){
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
		
		for(int i = 1; i < nCells + 1; i++){
			xL[i][0] = xL[i][1];
			xL[i][nCells + 1] = xL[i][nCells];
			
			xR[i][0] = xR[i][1];
			xR[i][nCells + 1] = xR[i][nCells];
		}
		
		std::vector<std::vector< std::vector<double>>> u_bar;
		u_bar.resize(nCells + 2,vector<vector<double>>(nCells + 2,vector<double>(4)));
		
		std::vector<std::vector< std::vector<double>>> halfxL;
		halfxL.resize(nCells + 2,vector<vector<double>>(nCells + 2,vector<double>(4)));
		
		std::vector<std::vector< std::vector<double>>> halfxR;
		halfxR.resize(nCells + 2,vector<vector<double>>(nCells + 2,vector<double>(4)));
		
		for(int i = 1; i < nCells + 1; i++){
			for(int j = 1; j < nCells + 1; j++){
				std::vector<double> fxL = f(xL[i][j], gama);
				std::vector<double> fxR = f(xR[i][j], gama);
				
				for(int k = 0; k < 4; k++){
					halfxL[i][j][k] = xL[i][j][k] - 0.5 * (dt/dx) * (fxR[k] - fxL[k]);
					halfxR[i][j][k] = xR[i][j][k] - 0.5 * (dt/dx) * (fxR[k] - fxL[k]);
				}
			}
		}
		
		for(int i = 1; i < nCells + 1; i++){
			halfxL[i][0] = halfxL[i][1];
			halfxL[i][nCells + 1] = halfxL[i][nCells];
			
			halfxR[i][0] = halfxR[i][1];
			halfxR[i][nCells + 1] = halfxR[i][nCells];
		}
		
		std::vector<std::vector< std::vector<double>>> flux_x;
		flux_x.resize(nCells + 2,vector<vector<double>>(nCells + 1,vector<double>(4)));
		
		for(int i = 1; i < nCells + 1; i++){
			for(int j = 0; j < nCells + 1; j++){
				flux_x[i][j] = FORCE_flux(halfxR[i][j], halfxL[i][j + 1], dt, dx, gama, 'x');
			}
		}

		for(int i = 1; i < nCells + 1; i++){
			for(int j = 1; j < nCells + 1; j++){
				for(int k = 0; k < 4; k++){
					u_bar[i][j][k] = u_con[i][j][k] - (dt/dx) * (flux_x[i][j][k] - flux_x[i][j - 1][k]);
				}
			}
		}
		
		
		
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		for(int i = 0; i < nCells + 2; i++){
			u_bar[i][0] = u_bar[i][1];
			u_bar[i][nCells + 1] = u_bar[i][nCells];
		}
		
		for(int i = 0; i < nCells + 2; i++){
			u_bar[0][i] = u_bar[1][i];
			u_bar[nCells + 1][i] = u_bar[nCells][i];
		}
		
		for(int i = 1; i < nCells + 1; i++){
			for(int j = 1; j < nCells + 1; j++){
				for(int k = 0; k < 4; k++){
					deltay0[i][j][k] = u_bar[i][j][k] - u_bar[i][j - 1][k];
					deltay1[i][j][k] = u_bar[i][j + 1][k] - u_bar[i][j][k];
				}
			}
		}
		
		for(int i = 1; i < nCells + 1; i++){
			deltay0[0][i] = deltay0[1][i];
			deltay0[nCells + 1][i] = deltay0[nCells][i];
			deltay1[0][i] = deltay1[1][i];
			deltay1[nCells + 1][i] = deltay1[nCells][i];
		}
		
		for(int i = 1; i < nCells + 1; i++){
			for(int j = 1; j < nCells + 1; j++){
				for(int k = 0; k < 4; k++){
					deltay[i][j][k] = 0.5 * (1 + wy[k]) * deltay0[i][j][k] + 0.5 * (1 - wy[k]) * deltay1[i][j][k];
				}
			}
		}
		
		for(int i = 1; i < nCells + 1; i++){
			deltay[0][i] = deltay[1][i];
			deltay[nCells + 1][i] = deltay[nCells][i];
		}
		
		std::vector<std::vector< std::vector<double>>> yL;
		yL.resize(nCells + 2,vector<vector<double>>(nCells + 2,vector<double>(4)));
		
		std::vector<std::vector< std::vector<double>>> yR;
		yR.resize(nCells + 2,vector<vector<double>>(nCells + 2,vector<double>(4)));

		for(int i = 1; i < nCells + 1; i++){
			for(int j = 1; j < nCells + 1; j++){
				std::vector<double> ry;
				ry.resize(4);
				for(int k = 0; k < 4; k++){
					ry[k] = (u_bar[i][j][k] - u_bar[i - 1][j][k])/(u_bar[i + 1][j][k] - u_bar[i][j][k]);
					if(u_bar[i + 1][j][k] - u_bar[i][j][k] == 0){
						ry[k] = 0;
					}
					yL[i][j][k] = u_bar[i][j][k] - 0.5 * Minbee(ry[k]) * deltay[i][j][k];
					yR[i][j][k] = u_bar[i][j][k] + 0.5 * Minbee(ry[k]) * deltay[i][j][k];
				}
			}
		}
	
		for(int i = 1; i < nCells + 1; i++){
				yL[0][i] = yL[1][i];
				yL[nCells + 1][i] = yL[nCells][i];
				
				yR[0][i] = yR[1][i];
				yR[nCells + 1][i] = yR[nCells][i];
		}
				
		std::vector<std::vector< std::vector<double>>> u_con_plus;
		u_con_plus.resize(nCells + 2,vector<vector<double> >(nCells + 2,vector<double>(4)));
		
		std::vector<std::vector< std::vector<double>>> halfyL;
		halfyL.resize(nCells + 2,vector<vector<double>>(nCells + 2,vector<double>(4)));
		
		std::vector<std::vector< std::vector<double>>> halfyR;
		halfyR.resize(nCells + 2,vector<vector<double>>(nCells + 2,vector<double>(4)));

		for(int i = 1; i < nCells + 1; i++){
			for(int j = 1; j < nCells + 1; j++){
				std::vector<double> gyL = g(yL[i][j], gama);
				std::vector<double> gyR = g(yR[i][j], gama);
				
				for(int k = 0; k < 4; k++){
					halfyL[i][j][k] = yL[i][j][k] - 0.5 * (dt/dy) * (gyR[k] - gyL[k]);
					halfyR[i][j][k] = yR[i][j][k] - 0.5 * (dt/dy) * (gyR[k] - gyL[k]);
				}
			}
		}

		for(int i = 1; i < nCells + 1; i++){
			halfyL[0][i] = halfyL[1][i];
			halfyL[nCells + 1][i] = halfyL[nCells][i];
			
			halfyR[0][i] = halfyR[1][i];
			halfyR[nCells + 1][i] = halfyR[nCells][i];
		}
		
		std::vector<std::vector< std::vector<double>>> flux_y;
		flux_y.resize(nCells + 1,vector<vector<double>>(nCells + 2,vector<double>(4)));
		
		for(int j = 1; j < nCells + 1; j++){
			for(int i = 0; i < nCells + 1; i++){
				flux_y[i][j] = FORCE_flux(halfyR[i][j], halfyL[i + 1][j], dt, dy, gama, 'y');
			}
		}
		
		for(int j = 1; j < nCells + 1; j++){
			for(int i = 1; i < nCells + 1; i++){
				for(int k = 0; k < 4; k++){
					u_con_plus[i][j][k] = u_bar[i][j][k] - (dt/dy) * (flux_y[i][j][k] - flux_y[i - 1][j][k]);
				}
			}
		}
		
		u_con = u_con_plus;
		
		for(int i = 1; i < nCells + 1; i++){
			for(int j = 1; j < nCells + 1; j++){
				u_pri[i][j] = con2pri(u_con[i][j], gama);
			}
		}
		std::cout << t << endl;
	}while(t < time);
	
	ofstream outFilerho("/Users/chenzefeng/Desktop/rho.dat");
	for(int i = 1; i < nCells + 1; i++){
		for(int j = 1; j < nCells + 1; j++){
			outFilerho << x0 + dx * (i - 0.5) << " " << y0 + dy * (j - 0.5) << " " << u_pri[i][j][0] << endl;
		}
	}
	outFilerho.close();
	
	ofstream outFilevx("/Users/chenzefeng/Desktop/vx.dat");
	for(int i = 1; i < nCells + 1; i++){
		for(int j = 1; j < nCells + 1; j++){
			outFilevx << x0 + dx * (i - 0.5) << " " << y0 + dy * (j - 0.5) << " " << u_pri[i][j][1] << endl;
		}
	}
	outFilevx.close();
	
	ofstream outFilevy("/Users/chenzefeng/Desktop/vy.dat");
	for(int i = 1; i < nCells + 1; i++){
		for(int j = 1; j < nCells + 1; j++){
			outFilevy << x0 + dx * (i - 0.5) << " " << y0 + dy * (j - 0.5) << " " << u_pri[i][j][2] << endl;
		}
	}
	outFilevy.close();
	
	ofstream outFilep("/Users/chenzefeng/Desktop/p.dat");
	for(int i = 1; i < nCells + 1; i++){
		for(int j = 1; j < nCells + 1; j++){
			outFilep << x0 + dx * (i - 0.5) << " " << y0 + dy * (j - 0.5) << " " << u_pri[i][j][3] << endl;
		}
	}
	outFilep.close();
	
	ofstream outFileE("/Users/chenzefeng/Desktop/E.dat");
	for(int i = 1; i < nCells + 1; i++){
		for(int j = 1; j < nCells + 1; j++){
			outFileE << x0 + dx * (i - 0.5) << " " << y0 + dy * (j - 0.5) << " " << u_pri[i][j][3] << endl;
		}
	}
	outFileE.close();
}