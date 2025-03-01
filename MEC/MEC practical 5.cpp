#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>
#include <fstream>

double Minbee(double r);
std::vector<double> f(std::vector<double>con,double gama, double Bx);
std::vector<double> con2pri(std::vector<double>con, double gama, double Bx);
std::vector<double> pri2con(std::vector<double>pri, double gama, double Bx);
std::vector<double> Fri_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama, double Bx);
std::vector<double> RI_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama, double Bx);
std::vector<double> FORCE_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama, double Bx);
double calecon(std::vector<double> con, double Bx);
double calepri(std::vector<double> pri, double gama);
double calcf(std::vector<double> pri, double gama, double Bx);

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

std::vector<double> f(std::vector<double>con, double gama, double Bx){
	// calculate f
	// based on conservative value roh rohv E
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double vz = con[3] / con[0];
	double U = con[4];
	double By = con[5];
	double Bz = con[6];
	std::vector<double> pri = con2pri(con, gama, Bx);
	double p = pri[4];
	std::vector<double> f;
	f.resize(7);
	double B_square = Bx * Bx + By * By + Bz * Bz;
	double vB = vx * Bx + vy * By + vz * Bz;
	f[0] = rho * vx;
	f[1] = rho * vx * vx + p + 0.5 * B_square - Bx * Bx;
	f[2] = rho * vy * vx - By * Bx;
	f[3] = rho * vz * vx - Bz * Bx;
	f[4] = (U + p + 0.5 * B_square) * vx - vB * Bx;
	f[5] = By * vx - vy * Bx;
	f[6] = Bz * vx - vz * Bx;
	return f;
}

std::vector<double> con2pri(std::vector<double>con, double gama, double Bx){
	// convert conservative (rho rhov E) to primitive (rho v p)
	std::vector<double> pri;
	pri.resize(7);
	// rho v p
	double e = calecon(con, Bx);
	pri[0] = con[0];
	pri[1] = con[1] / con[0];
	pri[2] = con[2] / con[0];
	pri[3] = con[3] / con[0];
	pri[4] = (gama - 1) * con[0] * e;
	pri[5] = con[5];
	pri[6] = con[6];
	return pri;
}

std::vector<double> pri2con(std::vector<double>pri, double gama, double Bx){
	// convert primitive (rho v p) to conservative  (rho rhov E)
	std::vector<double> con;
	con.resize(7);
	double e = calepri(pri, gama);
	// rho rhov E
	double vx = pri[1];
	double vy = pri[2];
	double vz = pri[3];
	double By = pri[5];
	double Bz = pri[6];
	double B_square = Bx * Bx + By * By + Bz * Bz;
	double v_square = vx * vx + vy * vy + vz * vz;
	con[0] = pri[0];
	con[1] = pri[0] * pri[1];
	con[2] = pri[0] * pri[2];
	con[3] = pri[0] * pri[3];
	con[4] = pri[0] * e + 0.5 * pri[0] * v_square + 0.5 * B_square;
	con[5] = pri[5];
	con[6] = pri[6];
	return con;
}

std::vector<double> Fri_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama, double Bx){
	// calculate the flux for lax-friedrich
	std::vector<double> Fri;
	Fri.resize(7);
	std::vector<double> f0 = f(con0, gama, Bx);
	std::vector<double> f1 = f(con1, gama, Bx);
	for(int i = 0; i < 7; i++){
		Fri[i] = 0.5 * (dx/dt) * (con0[i] - con1[i]) + 0.5 * (f1[i] + f0[i]);
	}
	return Fri;
}

std::vector<double> RI_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama, double Bx){
	std::vector<double> RI;
	std::vector<double> u;
	u.resize(7);
	std::vector<double> f0 = f(con0, gama, Bx);
	std::vector<double> f1 = f(con1, gama, Bx);
	for(int i = 0; i < 7; i++){
		u[i] = 0.5 * (con0[i] + con1[i]) - 0.5 * (dt/dx) * (f1[i] - f0[i]);
	}
	RI = f(u, gama, Bx);
	return RI;
}

std::vector<double> FORCE_flux(std::vector<double> con0, std::vector<double> con1, double dt, double dx, double gama, double Bx){
	std::vector<double> flux;
	flux.resize(7);
	std::vector<double> RI = RI_flux(con0, con1, dt, dx, gama, Bx);
	std::vector<double> Fri = Fri_flux(con0, con1, dt, dx, gama, Bx);
	for(int i = 0; i < 7 ; i++){
		flux[i] = 0.5 * (Fri[i] + RI[i]);
	}
	return flux;
}

double calecon(std::vector<double> con, double Bx){
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double vz = con[3] / con[0];
	double U = con[4];
	double By = con[5];
	double Bz = con[6];
	double B_square = Bx * Bx + By * By + Bz * Bz;
	double v_square = vx * vx + vy * vy + vz * vz;
	double e = (U - 0.5 * rho * v_square - 0.5 * B_square)/rho;
	return e;
}

double calepri(std::vector<double> pri, double gama){
	double e = pri[4]/((gama - 1) * pri[0]);
	return e;
}

double calcf(std::vector<double> pri, double gama, double Bx){
	double cs = sqrt(gama * pri[4] / pri[0]);
	double ca = fabs(Bx) / sqrt(pri[0]);
	double cf = sqrt(0.5 * ((cs * cs + ca * ca) + sqrt((cs * cs + ca * ca) * (cs * cs + ca * ca) - ((4 * (cs * cs * Bx * Bx))/pri[0]))));
	return cf;
}

using namespace std;
int main(int argc, char *argv[]) {
	double x0 = -0.75;
	double x1 = 0.75;
	double gama = 5.0 / 3.0;
	double c = 0.8;
	int nCells = 2048;
	double nxPoints = nCells + 1;
	double time = 0.2;
	double t = 0;
	double dx = (x1 - x0)/(nxPoints-2);
	
	std::vector<std::vector<double>> u_pri;
	u_pri.resize(nCells + 4,vector<double>(7));

	std::vector<std::vector<double>> u_con;
	u_con.resize(nCells + 4,vector<double>(7));
	
	std::vector<double> wx;
	wx.resize(7,0);
	
	std::vector<double> Bx;
	Bx.resize(nCells + 4,0);
	
	std::vector<std::vector<double>> deltax0;
	deltax0.resize(nCells + 4,vector<double>(7));
	
	std::vector<std::vector<double>> deltax1;
	deltax1.resize(nCells + 4,vector<double>(7));
	
	std::vector<std::vector<double>> deltax;
	deltax.resize(nCells + 4, vector<double>(7));
	
//	for(int i = 2; i < nCells + 2; i++){
//		double x = x0 + (i - 1.5) * dx;
//		if(x <= 0.5){
//			u_pri[i][0] = 1.08;
//			u_pri[i][1] = 1.2;
//			u_pri[i][2] = 0.01;
//			u_pri[i][3] = 0.5;
//			u_pri[i][4] = 3.6 / sqrt(4 * M_PI);
//			Bx[i] = 2 / sqrt(4 * M_PI);
//			u_pri[i][5] = 2 / sqrt(4 * M_PI);
//			u_pri[i][6] = 0.95;
//		}else if(x > 0.5){
//			u_pri[i][0] = 1;
//			u_pri[i][1] = 0;
//			u_pri[i][2] = 0;
//			u_pri[i][3] = 0;
//			u_pri[i][4] = 4 / sqrt(4 * M_PI);
//			Bx[i] = 2 / sqrt(4 * M_PI);
//			u_pri[i][5] = 2 / sqrt(4 * M_PI);
//			u_pri[i][6] = 1;
//		}
//	}
	
	for(int i = 2; i < nCells + 2; i++){
		double x = x0 + (i - 1.5) * dx;
		if(x <= 0){
			u_pri[i][0] = 1;
			u_pri[i][1] = 0;
			u_pri[i][2] = 0;
			u_pri[i][3] = 0;
			u_pri[i][4] = 1;
			Bx[i] = 0.75;
			u_pri[i][5] = 1;
			u_pri[i][6] = 0;
		}else if(x > 0){
			u_pri[i][0] = 0.125;
			u_pri[i][1] = 0;
			u_pri[i][2] = 0;
			u_pri[i][3] = 0;
			u_pri[i][4] = 0;
			Bx[i] = 0.75;
			u_pri[i][5] = 0;
			u_pri[i][6] = 0;
		}
	}
//
	
//	for(int i = 2; i < nCells + 2; i++){
//		double x = x0 + (i - 1.5) * dx;
//		if(x <= 0.5){
//			u_pri[i][0] = 1;
//			u_pri[i][1] = 0;
//			u_pri[i][2] = 0;
//			u_pri[i][3] = 0;
//			u_pri[i][4] = 1;
//			Bx[i] = 0;
//			u_pri[i][5] = 0;
//			u_pri[i][6] = 0;
//		}else if(x > 0.5){
//			u_pri[i][0] = 0.125;
//			u_pri[i][1] = 0;
//			u_pri[i][2] = 0;
//			u_pri[i][3] = 0;
//			u_pri[i][4] = 0.1;
//			Bx[i] = 0;
//			u_pri[i][5] = 0;
//			u_pri[i][6] = 0;
//		}
//	}

	
	do{
		for(int k = 0; k < 7; k++){
			u_pri[1][k] = u_pri[2][k];
			u_pri[0][k] = u_pri[2][k];
			u_pri[nCells + 2][k] = u_pri[nCells + 1][k];
			u_pri[nCells + 3][k] = u_pri[nCells + 1][k];
		}
		
		Bx[1] = Bx[2];
		Bx[0] = Bx[2];
		Bx[nCells + 2] = Bx[nCells + 1];
		Bx[nCells + 3] = Bx[nCells + 1];
		
		for(int i = 0; i < nCells + 4; i++){
			u_con[i] = pri2con(u_pri[i], gama, Bx[i]);
		}
		
		double a_max = -1e8;
		
		for(int i = 0; i < nCells + 4; i++){
			double cf = calcf(u_pri[i], gama, Bx[i]);
			double v = sqrt(fabs(u_pri[i][1] * u_pri[i][1] + u_pri[i][2] * u_pri[i][2] + u_pri[i][3] * u_pri[i][3]));
			if(a_max < v + cf){
				a_max = v + cf;
			}
		}
		
		double dt = c * dx / a_max;
		t = t + dt;
				
		for(int i = 1; i < nCells + 3; i++){
			for(int k = 0; k < 7; k++){
				deltax0[i][k] = u_con[i][k] - u_con[i - 1][k];
				deltax1[i][k] = u_con[i + 1][k] - u_con[i][k];
			}
		}
		
		for(int i = 1; i < nCells + 3; i++){
			for(int k = 0; k < 7; k++){
				deltax[i][k] = 0.5 * (1 + wx[k]) * deltax0[i][k] + 0.5 * (1 - wx[k]) * deltax1[i][k];
			}
		}

		std::vector<std::vector<double>> xL;
		xL.resize(nCells + 4,vector<double>(7));
		
		std::vector<std::vector<double>> xR;
		xR.resize(nCells + 4,vector<double>(7));
			
		for(int i = 1; i < nCells + 3; i++){
			std::vector<double> rx;
			rx.resize(7);
			for(int k = 0; k < 7; k++){
				rx[k] = (u_con[i][k] - u_con[i - 1][k])/(u_con[i + 1][k] - u_con[i][k]);
				if(u_con[i + 1][k] - u_con[i][k] == 0){
					rx[k] = 0;
				}
				xL[i][k] = u_con[i][k] - 0.5 * Minbee(rx[k]) * deltax[i][k];
				xR[i][k] = u_con[i][k] + 0.5 * Minbee(rx[k]) * deltax[i][k];
			}
		}
		
		std::vector<std::vector<double>> u_bar;
		u_bar.resize(nCells + 4,vector<double>(7));
		
		std::vector<std::vector<double>> halfxL;
		halfxL.resize(nCells + 4,vector<double>(7));
		
		std::vector<std::vector<double>> halfxR;
		halfxR.resize(nCells + 4,vector<double>(7));
		
		for(int i = 1; i < nCells + 3; i++){
			std::vector<double> fxL = f(xL[i], gama, Bx[i]);
			std::vector<double> fxR = f(xR[i], gama, Bx[i]);
				
			for(int k = 0; k < 7; k++){
				halfxL[i][k] = xL[i][k] - 0.5 * (dt/dx) * (fxR[k] - fxL[k]);
				halfxR[i][k] = xR[i][k] - 0.5 * (dt/dx) * (fxR[k] - fxL[k]);
			}
		}
		
		std::vector<std::vector<double>> flux_x;
		flux_x.resize(nCells + 3,vector<double>(7));
		
		for(int i = 1; i < nCells + 2; i++){
			flux_x[i] = FORCE_flux(halfxR[i], halfxL[i + 1], dt, dx, gama, Bx[i]);
		}
		
		for(int i = 2; i < nCells + 2; i++){
			for(int k = 0; k < 7; k++){
				u_bar[i][k] = u_con[i][k] - (dt/dx) * (flux_x[i][k] - flux_x[i - 1][k]);
			}
		}
		
		u_con = u_bar;
		
		for(int i = 2; i < nCells + 2; i++){
			u_pri[i] = con2pri(u_con[i], gama, Bx[i]);
		}
//		std::cout << dt << std::endl;
	}while(t < time);
	
	ofstream outFilerho("/Users/chenzefeng/Desktop/rho.dat");
	for(int i = 2; i < nCells + 2; i++){
		outFilerho << x0 + dx * (i - 2) << " " << u_pri[i][0] << endl;
	}
	
	outFilerho.close();
	
	ofstream outFilevx("/Users/chenzefeng/Desktop/vx.dat");
	for(int i = 2; i < nCells + 2; i++){
		outFilevx << x0 + dx * (i - 2) << " " << u_pri[i][1] << endl;
	}
	
	outFilevx.close();
	
	ofstream outFilevy("/Users/chenzefeng/Desktop/vy.dat");
	for(int i = 2; i < nCells + 2; i++){
		outFilevy << x0 + dx * (i - 2) << " " << u_pri[i][2] << endl;
	}
	
	outFilevy.close();
	
	ofstream outFilevz("/Users/chenzefeng/Desktop/vz.dat");
	for(int i = 2; i < nCells + 2; i++){
		outFilevz << x0 + dx * (i - 2) << " " << u_pri[i][3] << endl;
	}
	
	outFilevz.close();
	
	ofstream outFilep("/Users/chenzefeng/Desktop/p.dat");
	for(int i = 2; i < nCells + 2; i++){
		outFilep << x0 + dx * (i - 2) << " " << u_pri[i][4] << endl;
	}
	
	outFilep.close();
	
	ofstream outFileBy("/Users/chenzefeng/Desktop/By.dat");
	for(int i = 2; i < nCells + 2; i++){
		outFileBy << x0 + dx * (i - 2) << " " << u_pri[i][5] << endl;
	}
	
	outFileBy.close();
	
	ofstream outFileBz("/Users/chenzefeng/Desktop/Bz.dat");
	for(int i = 2; i < nCells + 2; i++){
		outFileBz << x0 + dx * (i - 2) << " " << u_pri[i][6] << endl;
	}
	
	outFileBz.close();
}