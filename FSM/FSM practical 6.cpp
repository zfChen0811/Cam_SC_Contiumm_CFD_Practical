#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>
#include <fstream>

std::vector<double> f(std::vector<double>con, double gama1, double gama2, double p_inf1, double p_inf2);
std::vector<double> con2pri(std::vector<double>con, double gama1, double gama2, double p_inf1, double p_inf2);
std::vector<double> pri2con(std::vector<double> pri, double gama1, double gama2, double p_inf1, double p_inf2);
double calecon(std::vector<double> con, double alpha1, double alpha2, double gama1, double gama2, double p_inf1, double p_inf2);
double calepri(std::vector<double> pri, double gama, double p_inf);
double calcs(std::vector<double> pri, double gama1, double gama2, double p_inf1, double p_inf2);
double calgama(double alpha1, double alpha2, double gama1, double gama2);
double calp_inf(double alpha1, double alpha2, double gama1, double gama2, double p_inf1, double p_inf2);
double calrho(std::vector<double> con);
double calv(std::vector<double> con);
double calp(std::vector<double> con, double alpha1, double alpha2, double gama1, double gama2, double p_inf1, double p_inf2);

std::vector<double> f(std::vector<double>con, double gama1, double gama2, double p_inf1, double p_inf2){
	// calculate f
	// based on conservative value roh rohv E
	double alpha1 = con[0];
	double alpha2 = 1 - alpha1;
	double rho = calrho(con);
	double v = con[3] / rho;
	double E = con[4];
	double p = calp(con, alpha1, alpha2, gama1, gama2, p_inf1, p_inf2);
	std::vector<double> f;
	f.resize(5);
	f[0] = alpha1 * v;
	f[1] = con[1] * v;
	f[2] = con[2] * v;
	f[3] = rho * v * v + p;
	f[4] = alpha1 * (E + p) * v;
	return f;
}

std::vector<double> con2pri(std::vector<double>con, double gama1, double gama2, double p_inf1, double p_inf2){
	// convert conservative (rho rhov E) to primitive (rho v p)
	std::vector<double> pri;
	pri.resize(5);
	// rho v p
	double alpha1 = con[0];
	double alpha2 = 1 - alpha1;
	double gama = calgama(alpha1, alpha2, gama1, gama2);
	double p_inf = calp_inf(alpha1, alpha2, gama1, gama2, p_inf1, p_inf2);
	double rho = calrho(con);
	double v = calv(con);
	double E = con[4];
	double e = calecon(con, alpha1, alpha2, gama1, gama2, p_inf1, p_inf2);
	double rho1 = con[1] / alpha1;
	double rho2 = con[2] / alpha2;
	pri[0] = alpha1;
	pri[1] = rho1;
	pri[2] = rho2;
	pri[3] = v;
	pri[4] = calp(con, alpha1, alpha2, gama1, gama2, p_inf1, p_inf2);
	return pri;
}

std::vector<double> pri2con(std::vector<double> pri, double gama1, double gama2, double p_inf1, double p_inf2){
	// convert primitive (rho v p) to conservative  (rho rhov E)
	std::vector<double> con;
	con.resize(5);
	double alpha1 = pri[0];
	double alpha2 = 1 - pri[0];
	double gama = calgama(alpha1, alpha2, gama1, gama2);
	double p_inf = calp_inf(alpha1, alpha2, gama1, gama2, p_inf1, p_inf2);
	double e = calepri(pri, gama, p_inf);
	double rho1 = pri[1];
	double rho2 = pri[2];
	double rho = alpha1 * rho1 + alpha2 * rho2;
	double v = pri[3];
	con[0] = alpha1;
	con[1] = alpha1 * rho1;
	con[2] = alpha2 * rho2;
	con[3] = rho * v;
	con[4] = rho * e + 0.5 * rho * (v * v);
	return con;
}

double calecon(std::vector<double> con, double alpha1, double alpha2, double gama1, double gama2, double p_inf1, double p_inf2){
	double rho = calrho(con);
	double v = calv(con);
	double E = con[4];
	double gama = calgama(alpha1, alpha2, gama1, gama2);
	double e = (E - 0.5 * rho * (v * v))/rho;
	return e;
}

double calepri(std::vector<double> pri, double gama, double p_inf){
	double alpha1 = pri[0];
	double alpha2 = 1 - alpha1;
	double rho1 = pri[1];
	double rho2 = pri[2];
	double rho = alpha1 * rho1 + alpha2 * rho2;
	double p = pri[4];
	double e = (p + gama * p_inf)/((gama - 1) * rho);
	return e;
}

double calcs(std::vector<double> pri, double gama1, double gama2, double p_inf1, double p_inf2){
	double alpha1 = pri[0];
	double alpha2 = 1 - alpha1;
	double rho1 = pri[1];
	double rho2 = pri[2];
	double gama = calgama(alpha1, alpha2, gama1, gama2);
	double p_inf = calp_inf(alpha1, alpha2, gama1, gama2, p_inf1, p_inf2);
	double rho = alpha1 * rho1 + alpha2 * rho2;
	std::vector<double> con = pri2con(pri, gama1, gama2, p_inf1, p_inf2);
	double p = calp(con, alpha1, alpha2, gama1, gama2, p_inf1, p_inf2);
	double cs = sqrt(gama *  (p + p_inf)/ rho);
	return cs;
}

double calxi(double alpha1, double alpha2, double gama1, double gama2){
	double xi = alpha1 / (gama1 - 1) + alpha2 / (gama2 - 1);
	return xi;
}


double calgama(double alpha1, double alpha2, double gama1, double gama2){
	double gama = 1 / calxi(alpha1, alpha2, gama1, gama2) + 1;
	return gama;
}

double calp_inf(double alpha1, double alpha2, double gama1, double gama2, double p_inf1, double p_inf2){
	double gama = calgama(alpha1, alpha2, gama1, gama2); 
	double gamap_inf = (alpha1 * gama1 * p_inf1 / (gama1 - 1) + alpha2 * gama2 * p_inf2 / (gama2 - 1)) * (gama - 1);
	double p_inf = gamap_inf / gama;
	return p_inf;
}

double calp(std::vector<double> con, double alpha1, double alpha2, double gama1, double gama2, double p_inf1, double p_inf2){
	double p_inf = calp_inf(alpha1, alpha2, gama1, gama2, p_inf1, p_inf2);
	double gama = calgama(alpha1, alpha2, gama1, gama2);
	double rho = calrho(con); 
	double e = calecon(con, alpha1, alpha2, gama1, gama2, p_inf1, p_inf2);
	double p = (gama - 1) * rho * e - gama * p_inf;
	return p;
}

double calrho(std::vector<double> con){
	double rho = con[1] + con[2];
	return rho;
}

double calv(std::vector<double> con){
	double rho = calrho(con);
	double v = con[3] / rho;
	return v;
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
	int nCells = 100;
	double nxPoints = nCells + 1;
	double time = 237.44 * pow(10, -6);
	double t = 0;
	double dx = (x1 - x0)/(nxPoints-1);
		
	std::vector<std::vector<double>> u_pri;
	u_pri.resize(nCells + 2,vector<double>(5));
	
	std::vector<std::vector<double>> u_con;
	u_con.resize(nCells + 2,vector<double>(5));
	
	std::vector<std::vector<double>> flux;
	flux.resize(nCells + 1,vector<double>(5));
	
	std::vector<double> S_star;
	S_star.resize(nCells + 1,0);
	
	std::vector<double> v_star;
	v_star.resize(nCells + 1,0);
	
	std::vector<double> alphaL1X;
	alphaL1X.resize(nCells + 1,0);
	
	std::vector<double> alphaL2X;
	alphaL2X.resize(nCells + 1,0);
	
	std::vector<double> alphaR1X;
	alphaR1X.resize(nCells + 1,0);
	
	std::vector<double> alphaR2X;
	alphaR2X.resize(nCells + 1,0);
	
	std::vector<double> ELX;
	ELX.resize(nCells + 1,0);
	
	std::vector<double> ERX;
	ERX.resize(nCells + 1,0);
	
	std::vector<double> rhoLX;
	rhoLX.resize(nCells + 1,0);
	
	std::vector<double> rhoRX;
	rhoRX.resize(nCells + 1,0);
	
	std::vector<double> rhoL1X;
	rhoL1X.resize(nCells + 1,0);
	
	std::vector<double> rhoL2X;
	rhoL2X.resize(nCells + 1,0);
	
	std::vector<double> vLX;
	vLX.resize(nCells + 1,0);
	
	std::vector<double> pLX;
	pLX.resize(nCells + 1,0);
	
	std::vector<double> rhoR1X;
	rhoR1X.resize(nCells + 1,0);
	
	std::vector<double> rhoR2X;
	rhoR2X.resize(nCells + 1,0);
	
	std::vector<double> vRX;
	vRX.resize(nCells + 1,0);
	
	std::vector<double> pRX;
	pRX.resize(nCells + 1,0);
	
	std::vector<double> csLX;
	csLX.resize(nCells + 1,0);
	
	std::vector<double> csRX;
	csRX.resize(nCells + 1,0);
	
	std::vector<double> S_plusX;
	S_plusX.resize(nCells + 1,0);
	
	std::vector<std::vector<double>> fLX;
	fLX.resize(nCells + 1,vector<double>(5));
	
	std::vector<std::vector<double>> fRX;
	fRX.resize(nCells + 1,vector<double>(5));
	
	for(int i = 1; i < nCells + 1; i++){
		double x = x0 + (i - 0.5) * dx;
		if(x <= 0.7){
			u_pri[i][0] = 1 - 10e-12;
			u_pri[i][1] = 1000;
			u_pri[i][2] = 50;
			u_pri[i][3] = 0;
			u_pri[i][4] = 1 * pow(10, 9);
		}
		if(x >= 0.7){
			u_pri[i][0] = 10e-12;
			u_pri[i][1] = 1000;
			u_pri[i][2] = 50;
			u_pri[i][3] = 0;
			u_pri[i][4] = 1 * pow(10, 5);
		}
	}

	do{
		for(int k = 0; k < 5; k++){
			u_pri[0][k] = u_pri[1][k];
			u_pri[nCells + 1][k] = u_pri[nCells][k];
		}
				
		for(int i = 0; i < nCells + 2; i++){
			u_con[i] = pri2con(u_pri[i], gama1, gama2, p_inf1, p_inf2);
		}
		
		double a_max = -1e8;
		double cs;
		for(int i = 0; i < nCells + 2; i++){
			cs = calcs(u_pri[i], gama1, gama2, p_inf1, p_inf2);
			if(a_max < fabs(u_pri[i][3]) + cs){
				a_max = fabs(u_pri[i][3]) + cs;
			}
		}
		
		double dt = c * dx / a_max;
		t = t + dt;
		
		for(int i = 0; i < nCells + 1; i++){
			std::vector<double> uL = u_con[i];
			std::vector<double> uR = u_con[i + 1];
			std::vector<double> priL = u_pri[i];
			std::vector<double> priR = u_pri[i + 1];
			double alphaL1 = uL[0];
			double alphaL2 = 1 - alphaL1;
			double alphaR1 = uR[0];
			double alphaR2 = 1 - alphaR1;
			double EL = uL[4];
			double ER = uR[4];
			double rhoL = calrho(uL);
			double rhoR = calrho(uR);
			double rhoL1 = priL[1];
			double rhoL2 = priL[2];
			double vL = priL[3];
			double pL = priL[4];
			double rhoR1 = priR[1];
			double rhoR2 = priR[2];
			double vR = priR[3];
			double pR = priR[4];
			double csL = calcs(priL, gama1, gama2, p_inf1, p_inf2);
			double csR = calcs(priR, gama1, gama2, p_inf1, p_inf2);
			std::vector<double> fL = f(uL, gama1, gama2, p_inf1, p_inf2);
			std::vector<double> fR = f(uR, gama1, gama2, p_inf1, p_inf2);
			
			double S_plus = fmax(fabs(vL) + csL, fabs(vR) + csR);
			double SL = -S_plus;
			double SR = S_plus;
			
			
			
			
			alphaL1X[i] = alphaL1;
			alphaL2X[i] = alphaL2;
			alphaR1X[i] = alphaR1;
			alphaR2X[i] = alphaR2;
			ELX[i] = EL;
			ERX[i] = ER;
			rhoLX[i] = rhoL;
			rhoRX[i] = rhoR;
			rhoL1X[i] = rhoL1;
			rhoL2X[i] = rhoL2;
			vLX[i] = vL;
			pLX[i] = pL;
			rhoR1X[i] = rhoR1;
			rhoR2X[i] = rhoR2;
			vRX[i] = vR;
			pRX[i] = pR;
			csLX[i] = csL;
			csRX[i] = csR;
			S_plusX[i] = S_plus;
			fLX[i] = fL;
			fRX[i] = fR;
			
			
			
			
			S_star[i] = (pR - pL + rhoL * vL * (SL - vL) - rhoR * vR * (SR - vR))/(rhoL * (SL - vL) - rhoR * (SR - vR));
			
			std::vector<double> uL_HLLC;
			uL_HLLC.resize(5,0);
			std::vector<double> uR_HLLC;
			uR_HLLC.resize(5,0);
			
			uL_HLLC[0] = alphaL1;
			uL_HLLC[1] = alphaL1 * rhoL1 * ((SL - vL) / (SL - S_star[i]));
			uL_HLLC[2] = alphaL2 * rhoL2 * ((SL - vL) / (SL - S_star[i]));
			uL_HLLC[3] = rhoL * ((SL - vL)/(SL - S_star[i])) * S_star[i];
			uL_HLLC[4] = rhoL * ((SL - vL)/(SL - S_star[i])) * ((EL / rhoL) + (S_star[i] - vL) * (S_star[i] + pL / (rhoL * (SL - vL))));
			std::vector<double> priL_HLLC = con2pri(uL_HLLC, gama1, gama2, p_inf1, p_inf2);
			
			uR_HLLC[0] = alphaR1;
			uR_HLLC[1] = alphaR1 * rhoR1 * ((SR - vR) / (SR - S_star[i]));
			uR_HLLC[2] = alphaR2 * rhoR2 * ((SR - vR) / (SR - S_star[i]));
			uR_HLLC[3] = rhoR * ((SR - vR)/(SR - S_star[i])) * S_star[i];
			uR_HLLC[4] = rhoR * ((SR - vR)/(SR - S_star[i])) * ((ER / rhoR) + (S_star[i] - vR) * (S_star[i] + pR / (rhoR * (SR - vR))));
			std::vector<double> priR_HLLC = con2pri(uR_HLLC, gama1, gama2, p_inf1, p_inf2);
			
			if(SL >= 0){
				flux[i] = fL;
				v_star[i] = priL[3];
			}else if(SL < 0 && S_star[i] >= 0){
				for(int k = 0; k < 5; k++){
					flux[i][k] = fL[k] + SL * (uL_HLLC[k] - uL[k]);
				}
				v_star[i] = priL_HLLC[3];
			}else if(S_star[i] < 0 && SR >= 0){
				for(int k = 0; k < 5; k++){
					flux[i][k] = fR[k] + SR * (uR_HLLC[k] - uR[k]);
				}
				v_star[i] = priR_HLLC[3];
			}else{
				flux[i] = fR;
				v_star[i] = priR[3];
			}
		}
		
		std::vector<std::vector<double>> u_bar;
		u_bar.resize(nCells + 2,vector<double>(5));
		
		for(int i = 1; i < nCells + 1; i++){
			for(int k = 0; k < 5; k++){
				u_bar[i][k] = u_con[i][k] - (dt/dx) * (flux[i][k] - flux[i - 1][k]);
			}
		}
		
		for(int i = 1; i < nCells + 1; i++){
			u_con[i][0] = u_bar[i][0] + (dt/dx) * u_con[i][0] * (S_star[i] - S_star[i - 1]);
			for(int k = 1; k < 5; k++){
				u_con[i][k] = u_bar[i][k];
			}
		}
		
		for(int i = 1; i < nCells + 1; i++){
			u_pri[i] = con2pri(u_con[i], gama1, gama2, p_inf1, p_inf2);
		}
	}while(t < time);

	ofstream outFile("/Users/chenzefeng/Desktop/velocity.dat");
	for(int i = 1; i < nCells + 1; i++){
		outFile << x0 + dx * (i - 2) << " " << u_pri[i][1] << endl;
//		std::cout << calrho(u_con[i]) << endl;
//		std::cout << i << endl;
	}
	outFile.close();
}