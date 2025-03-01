#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>
#include <fstream>

std::vector<double> f(std::vector<double>con);
std::vector<double> con2pri(std::vector<double>con);
std::vector<double> pri2con(std::vector<double>pri);
double calecon(std::vector<double> con);
double calepri(std::vector<double> pri);
double calcs(std::vector<double> pri);

std::vector<double> f(std::vector<double>con){
	// calculate f
	// based on conservative value roh rohv E
	double gama = 1.4;
	double rho = con[0];
	double vx = con[1] / con[0];
	double E = con[2];
	std::vector<double> pri = con2pri(con);
	double p = pri[2];
	std::vector<double> f;
	f.resize(3);
	f[0] = rho * vx;
	f[1] = rho * vx * vx + p;
	f[2] = (E + p) * vx;
	return f;
}

std::vector<double> con2pri(std::vector<double>con){
	// convert conservative (rho rhov E) to primitive (rho v p)
	double gama = 1.4;
	std::vector<double> pri;
	pri.resize(3);
	// rho v p
	double e = calecon(con);
	pri[0] = con[0];
	pri[1] = con[1] / con[0];
	pri[2] = (gama - 1) * con[0] * e;
	return pri;
}

std::vector<double> pri2con(std::vector<double>pri){
	// convert primitive (rho v p) to conservative  (rho rhov E)
	std::vector<double> con;
	con.resize(3);
	double e = calepri(pri);
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

double calepri(std::vector<double> pri){
	double gama = 1.4;
	double e = pri[2]/((gama - 1) * pri[0]);
	return e;
}

double calcs(std::vector<double> pri){
	double gama = 1.4;
	double cs = sqrt(gama *  pri[2] / pri[0]);
	return cs;
}

using namespace std;
int main(int argc, char *argv[]) {
	double x0 = 0;
	double x1 = 1;
	double gama = 1.4;
	double c = 0.8;
	int nCells = 10;
	double nxPoints = nCells + 1;
	double time = 0.25;
	double t = 0;
	double dx = (x1 - x0)/(nxPoints-1);
	
	std::vector<std::vector<double>> u_pri;
	u_pri.resize(nCells + 2,vector<double>(3));
	
	std::vector<std::vector<double>> u_con;
	u_con.resize(nCells + 2,vector<double>(3));
	
	std::vector<std::vector<double>> flux;
	flux.resize(nCells + 1,vector<double>(3));
	
	std::vector<double> v_star;
	v_star.resize(nCells + 1,0);
	
	for(int i = 1; i < nCells + 1; i++){
		double x = x0 + (i - 0.5) * dx;
		if(x < 0.5){
			u_pri[i][0] = 1;
			u_pri[i][1] = 0;
			u_pri[i][2] = 1;
		}else if(x >= 0.5){
			u_pri[i][0] = 0.125;
			u_pri[i][1] = 0;
			u_pri[i][2] = 0.1;
		}
	}
	
	do{
		for(int k = 0; k < 3; k++){
			u_pri[0][k] = u_pri[1][k];
			u_pri[nCells + 1][k] = u_pri[nCells][k];
		}
		
		for(int i = 0; i < nCells + 2; i++){
			u_con[i] = pri2con(u_pri[i]);
		}
		
		double a_max = -1e8;
		
		for(int i = 0; i < nCells + 1; i++){
			double cs = calcs(u_pri[i]);
			if(a_max < fabs(u_pri[i][1]) + cs){
				a_max = fabs(u_pri[i][1]) + cs;
			}
		}
		
		double dt = c * dx / a_max;
		t = t + dt;

		for(int i = 0; i < nCells + 1; i++){
			std::vector<double> uL = u_con[i];
			std::vector<double> uR = u_con[i + 1];
			std::vector<double> priL = u_pri[i];
			std::vector<double> priR = u_pri[i + 1];
			std::vector<double> fL = f(uL);
			std::vector<double> fR = f(uR);
			double pL = priL[2];
			double pR = priR[2];
			double EL = uL[2];
			double ER = uR[2];
			double rhoL = uL[0];
			double rhoR = uR[0];
			double vL = priL[1];
			double vR = priR[1];
			double csL = calcs(priL);
			double csR = calcs(priR);
			
//			double S_plus = fmax(fabs(vL) + csL, fabs(vR) + csR);
//			double SL = -S_plus;
//			double SR = S_plus;
			double SL = fmin(vL - csL, vR - csR);
			double SR = fmax(vL + csL, vR + csR);
			std::vector<double> uL_HLLC;
			uL_HLLC.resize(3,0);
			std::vector<double> uR_HLLC;
			uR_HLLC.resize(3,0);
			
			double S_star = (pR - pL + rhoL * vL * (SL - vL) - rhoR * vR * (SR - vR))/(rhoL * (SL - vL) - rhoR * (SR - vR));
			v_star[i] = S_star;
			
			uL_HLLC[0] = rhoL * ((SL - vL)/(SL - S_star)) * 1;
			uL_HLLC[1] = rhoL * ((SL - vL)/(SL - S_star)) * S_star;
			uL_HLLC[2] = rhoL * ((SL - vL)/(SL - S_star)) * ((EL / rhoL) + (S_star - vL) * (S_star + pL / (rhoL * (SL - vL))));
			
			uR_HLLC[0] = rhoR * ((SR - vR)/(SR - S_star)) * 1;
			uR_HLLC[1] = rhoR * ((SR - vR)/(SR - S_star)) * S_star;
			uR_HLLC[2] = rhoR * ((SR - vR)/(SR - S_star)) * ((ER / rhoR) + (S_star - vR) * (S_star + pR / (rhoR * (SR - vR))));
			
			std::vector<double> uL_HLLC_pri;
			uL_HLLC_pri.resize(3,0);
			std::vector<double> uR_HLLC_pri;
			uR_HLLC_pri.resize(3,0);
			uL_HLLC_pri = con2pri(uL_HLLC);
			uR_HLLC_pri = con2pri(uR_HLLC);
			
			std::cout << uL_HLLC_pri[2] << " " << uR_HLLC_pri[2] << " " << pL + rhoL * (SL - vL) * (S_star - vL)<< " " << pR + rhoR * (SR - vR) * (S_star - vR)<< endl;
			
			if(SL >= 0){
				flux[i] = fL;
			}else if(SL < 0 && S_star >=0){
				for(int k = 0; k < 3; k++){
					flux[i][k] = fL[k] + SL * (uL_HLLC[k] - uL[k]);
				}
			}else if(S_star < 0 && SR >= 0){
				for(int k = 0; k < 3; k++){
					flux[i][k] = fR[k] + SR * (uR_HLLC[k] - uR[k]);
				}
			}else{
				flux[i] = fR;
			}
			
		}
		
		std::vector<std::vector<double>> u_bar;
		u_bar.resize(nCells + 2,vector<double>(3));
		
		for(int i = 1; i < nCells + 1; i++){
			for(int k = 0; k < 3; k++){
				u_bar[i][k] = u_con[i][k] - (dt/dx) * (flux[i][k] - flux[i - 1][k]);
			}
		}
		
		u_con = u_bar;
		
		for(int i = 1; i < nCells + 1; i++){
			u_pri[i] = con2pri(u_con[i]);
		}
		
	}while(t < time);
	
	ofstream outFile("/Users/chenzefeng/Desktop/velocity.dat");
	for(int i = 1; i < nCells + 1; i++){
		outFile << x0 + dx * (i - 0.5 - 1) << " " << u_pri[i][0] << endl;
	}
	outFile.close();
}