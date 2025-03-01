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
double fK(double p_star, std::vector<double>uK, double gama, double p_inf);
double fK_d(double p_star, std::vector<double>uK, double gama, double p_inf);
double f_star(double p_star, std::vector<double>uL, std::vector<double>uR, double gama1, double gama2, double p_inf_1, double p_inf_2);
double f_star_d(double p_star, std::vector<double>uL, std::vector<double>uR, double gama1, double gama2, double p_inf_1, double p_inf_2);



double fK(double p_star, std::vector<double>uK, double gama, double p_inf){
	std::vector<double> pri = con2pri(uK, gama, p_inf);
	double rhoK = pri[0];
	double pK = pri[2];
	double f_value;
	double AK = 2 / ((gama + 1) * rhoK);
	double BK = ((gama - 1) / (gama + 1)) * pK + 2 * gama * p_inf / (gama + 1);
	double csK = calcs(pri, gama, p_inf);
	if(p_star > pK){
		f_value = (p_star - pK) * sqrt(AK/(p_star + BK));
	}else{
		f_value = ((2 * csK) / (gama - 1)) * (pow(((p_star + p_inf) / (pK + p_inf)), (gama - 1)/(2 * gama)) - 1);
	}
	return f_value;
}

double fK_d(double p_star, std::vector<double>uK, double gama, double p_inf){
	std::vector<double> pri = con2pri(uK, gama, p_inf);
	double rhoK = pri[0];
	double pK = pri[2];
	double f_d;
	double AK = 2 / ((gama + 1) * rhoK);
	double BK = ((gama - 1) / (gama + 1)) * pK + 2 * gama * p_inf / (gama + 1);
	double csK = calcs(pri, gama, p_inf);
	if(p_star > pK){
		f_d = sqrt(AK / (BK + p_star)) * (1 - ((p_star - pK) / (2 * (BK + p_star))));
	}else{
		f_d = (1 / (rhoK * csK)) * pow((p_star + p_inf) / (pK + p_inf), -(gama + 1) / (2 * gama)) ;
	}
	return f_d;
}

double f_star(double p_star, std::vector<double>uL, std::vector<double>uR, double gama1, double gama2, double p_inf_1, double p_inf_2){
	double vL = uL[1] / uL[0];
	double vR = uR[1] / uR[0];
	double dv = vR - vL;
	double f_value;
	f_value = fK(p_star, uR, gama2, p_inf_2) + fK(p_star, uL, gama1, p_inf_1) + dv;
	return f_value;
}

double f_star_d(double p_star, std::vector<double>uL, std::vector<double>uR, double gama1, double gama2, double p_inf_1, double p_inf_2){
	double f_d_value;
	f_d_value = fK_d(p_star, uR, gama2, p_inf_2) + fK_d(p_star, uL, gama1, p_inf_1);
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
	double cs = sqrt(gama * (pri[2] + p_inf)/ pri[0]);
	return cs;
}

using namespace std;
int main(int argc, char *argv[]) {
	double x0 = 0;
	double x1 = 1;
	double c = 1;
	int nCells = 10;
	double nxPoints = nCells + 1;
	double time = 0.0007;
//	double time = 0.25;
	double t = 0;
	double dx = (x1 - x0)/(nxPoints-1);
	double p_old, p_new;
	double eplison = 10e-8;
	double interface = 0.5;
	double gama1 = 1.4;
	double gama2 = 1.67;
	double p_inf_1 = 0;
	double p_inf_2 = 0;
	
	std::vector<std::vector<double>> u_pri;
	u_pri.resize(nCells + 2,vector<double>(3));
	
	std::vector<std::vector<double>> u_con;
	u_con.resize(nCells + 2,vector<double>(3));
	
	std::vector<std::vector<double>> flux;
	flux.resize(nCells + 1,vector<double>(3));
	
//	for(int i = 1; i < nCells + 1; i++){
//		double x = x0 + (i - 0.5) * dx;
//		if(x < interface){
//			u_pri[i][0] = 1000;
//			u_pri[i][1] = 0;
//			u_pri[i][2] = pow(10,9);
//		}else if(x >= interface){
//			u_pri[i][0] = 50;
//			u_pri[i][1] = 0;
//			u_pri[i][2] = pow(10,5);
//		}
//	}
	for(int i = 1; i < nCells + 1; i++){
		double x = x0 + (i - 0.5) * dx;
		if(x < 0.5){
			u_pri[i][0] = 1;
			u_pri[i][1] = 0;
			u_pri[i][2] = 1 * pow(10,5);
		}else if(x >= 0.5){
			u_pri[i][0] = 0.1379;
			u_pri[i][1] = 0;
			u_pri[i][2] = 1 * pow(10,5);
		}
	}

	do{
		for(int k = 0; k < 3; k++){
			u_pri[0][k] = u_pri[1][k];
			u_pri[nCells + 1][k] = u_pri[nCells][k];
		}
		
		for(int i = 0; i < nCells + 2; i++){
			double x = x0 + (i - 0.5) * dx;
			if(x < interface){
				u_con[i] = pri2con(u_pri[i], gama1, p_inf_1);
			}else if(x >= interface){
				u_con[i] = pri2con(u_pri[i], gama2, p_inf_2);
			}
		}
		
		double a_max = -1e8;
		double cs;
		
		for(int i = 1; i < nCells + 1; i++){
			double x = x0 + (i - 0.5) * dx;
			if(x < interface){
				cs = calcs(u_pri[i], gama1, p_inf_1);
			}else if(x >= interface){
				cs = calcs(u_pri[i], gama2, p_inf_2);
			}
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
			double pL = u_pri[i][2];
			double pR = u_pri[i + 1][2];
			double rhoL = uL[0];
			double rhoR = uR[0];
			double vL = priL[1];
			double vR = priR[1];
			double x = x0 + (i - 0.5) * dx;
			double csL, csR;
			double gamaL, gamaR;
			double p_inf_L, p_inf_R;
			std::vector<double> u_half;
			
			if(x < interface){
				gamaL = gama1;
				p_inf_L = p_inf_1;
			}else{
				gamaL = gama2;
				p_inf_L = p_inf_2;
			}
			
			if(x + dx< interface){
				gamaR = gama1;
				p_inf_R = p_inf_1;
			}else{
				gamaR = gama2;
				p_inf_R = p_inf_2;
			}
			
			csL = calcs(priL, gamaL, p_inf_L);
			csR = calcs(priR, gamaR, p_inf_R);
			
			std::vector<double> u_pri_half;
			u_pri_half.resize(3,0);
			
			p_new = pow((csL + csR - 0.5 * (0.5 * (gamaL + gamaR) - 1) * (vR - vL)) / (csL/pow(pL,(gamaL - 1) / (2 * gamaL)) + (csR/pow(pR, (gamaR - 1) / (2 * gamaR)))), (2 * 0.5 * (gamaL + gamaR)) / (0.5 * (gamaL + gamaR) - 1));
			
			do{
				p_old = p_new;
				p_new = p_old - f_star(p_old, uL, uR, gamaL, gamaR, p_inf_L, p_inf_R) / f_star_d(p_old, uL, uR, gamaL, gamaR, p_inf_L, p_inf_R);
			}while(fabs(p_new - p_old) / fabs(p_old) > eplison);
			
			double p_star = p_new;
			double v_star = 0.5 * (vL + vR) + 0.5 * (fK(p_star, uR, gamaR, p_inf_R) - fK(p_star, uL, gamaL, p_inf_L));
			
			double csL_star = csL * pow((p_star/pL), ((gamaL - 1)/ (2 * gamaL)));
			double S_HL = vL - csL;
			double S_TL = v_star - csL_star;
			
			double csR_star = csR * pow((p_star/pR), ((gamaR - 1)/ (2 * gamaR)));
			double S_HR = vR + csR;
			double S_TR = v_star + csR_star;
			
			double sL = vL - csL * sqrt(((gamaL + 1) / (2 * gamaL)) * (p_star / pL) + (gamaL - 1) / (2 * gamaL));
			double sR = vR + csR * sqrt(((gamaR + 1) / (2 * gamaR)) * (p_star / pR) + (gamaR - 1) / (2 * gamaR));
			
			std::vector<double> priL_star;
			priL_star.resize(3,0);
			
			std::vector<double> priR_star;
			priR_star.resize(3,0);
			
			if(p_star > pL){
				// left shock
				double rhoL_star_shock = rhoL * (((p_star / pL) +  ((gamaL - 1) / (gamaL + 1))) / (((gamaL - 1) / (gamaL + 1)) * (p_star / pL) + 1));
				priL_star[0] = rhoL_star_shock;
				priL_star[1] = v_star;
				priL_star[2] = p_star;
			}else{
				// left rarefaction
				double rhoL_star_rare = rhoL * pow((p_star/pL), (1/gamaL));
				priL_star[0] = rhoL_star_rare;
				priL_star[1] = v_star;
				priL_star[2] = p_star;
			}
			
			if(p_star > pR){
				// right shock
				double rhoR_star_shock = rhoR * (((p_star / pR) +  ((gamaR - 1) / (gamaR + 1))) / (((gamaR - 1) / (gamaR + 1)) * (p_star / pR) + 1));
				priR_star[0] = rhoR_star_shock;
				priR_star[1] = v_star;
				priR_star[2] = p_star;
			}else{
				// right rarefaction
				double rhoR_star_rare = rhoR * pow((p_star/pR), (1/gamaR));
				priR_star[0] = rhoR_star_rare;
				priR_star[1] = v_star;
				priR_star[2] = p_star;
			}
			
			if(v_star > 0){
				if(p_star > pL){
					// left shock
					if(sL > 0){
						u_pri_half = uL;
					}else{
						u_pri_half = priL_star;
					}
					
				}else{
					// left fan
					if(S_HL > 0){
						u_pri_half = uL;
					}else if(S_TL < 0){
						u_pri_half = priL_star;
					}else{
						std::vector<double> uL_rare;
						uL_rare.resize(3,0);
						uL_rare[0] = rhoL * pow((2 / (gamaL + 1) + ((gamaL - 1) / ((gamaL + 1) * csL)) * (vL)), (2 / (gamaL - 1)));
						uL_rare[1] = (2 / (gamaL + 1)) * (csL + ((gamaL - 1) / 2) * vL);
						uL_rare[2] = pL * pow(((2 / (gamaL + 1)) + ((gamaL - 1) / ((gamaL + 1) * csL)) * (vL)), (2 * gamaL) / (gamaL - 1));
						u_pri_half = uL_rare;
					}
				}
			}else if(v_star <= 0){
				if(p_star > pR){
					// right shock
					if(sR < 0){
						u_pri_half = uR;
					}else{
						u_pri_half = priR_star;
					}
				}else{
					// right fan
					if(S_HR < 0){
						u_pri_half = uR;
					}else if(S_TR > 0){
						u_pri_half = priR_star;
					}else{
						std::vector<double> uR_rare;
						uR_rare.resize(3,0);
						uR_rare[0] = rhoR * pow((2 / (gamaR + 1) - ((gamaR - 1) / ((gamaR + 1) * csR)) * (vR)), (2 / (gamaR - 1)));
						uR_rare[1] = (2 / (gamaR + 1)) * (-csR + ((gamaR - 1) / 2) * vR);
						uR_rare[2] = pR * pow(((2 / (gamaR + 1)) - ((gamaR - 1) / ((gamaR + 1) * csR)) * (vR)), (2 * gamaR) / (gamaR - 1));
						u_pri_half = uR_rare;
					}
				}
			}
			double gama_half, p_inf_half;
			if(x + 0.5 * dx < interface){
				gama_half = gama1;
				p_inf_half = p_inf_1;
			}else{
				gama_half = gama2;
				p_inf_half = p_inf_2;
			}
			u_half = pri2con(u_pri_half, gama_half, p_inf_half);
			flux[i] = f(u_half, gama_half, p_inf_half);
			
			//std::cout<< "stop" << std::endl;
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
			double x = x0 + (i - 0.5) * dx;
			double gama, p_inf;
			if(x + 0.5 * dx < interface){
				gama = gama1;
				p_inf = p_inf_1;
			}else{
				gama = gama2;
				p_inf = p_inf_2;
			}
			u_pri[i] = con2pri(u_con[i], gama, p_inf);
		}
		
	}while(t < time);
	
	ofstream outFile("/Users/chenzefeng/Desktop/velocity.dat");
	for(int i = 1; i < nCells + 1; i++){
		std::cout << u_pri[i][0] << std::endl;
		outFile << x0 + dx * (i - 0.5 - 1) << " " << u_pri[i][0] << endl;
	}
	outFile.close();
}