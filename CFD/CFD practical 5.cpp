#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>
#include <fstream>

double Minbee(double r);
std::vector<double> f(double rho, double v, double p, double E);
std::vector<double> con2pri(double rho, double rhov, double E, double gama);
std::vector<double> pri2con(double rho, double v, double p, double gama);
double FORCE(double dx, double dt, double rho0, double rho1, double rhov0, double rhov1, double E0, double E1);
std::vector<double> Fri_flux(std::vector<double> con0, std::vector<double> con1, double p0, double p1, double dt, double dx);
std::vector<double> RI_flux(std::vector<double> con0, std::vector<double> con1,double rho0, double rho1, double v0, double v1, double p0, double p1, double E0, double E1, double dt, double dx, double gama);
double calecon(double rho, double v, double E);
double calepri(double p, double rho, double gama);
double calcs(double p, double rho, double gama);

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



std::vector<double> f(double rho, double v, double p, double E){
	// calculate f
	// based on conservative value roh rohv E
	std::vector<double> f;
	f.resize(3);
	f[0] = rho * v;
	f[1] = rho * v * v + p;
	f[2] = (E + p) * v;
	return f;
}

std::vector<double> con2pri(double rho, double rhov, double E, double gama){
	// convert conservative (rho rhov E) to primitive (rho v p)
	std::vector<double> pri;
	pri.resize(3);
	// rho v p
	double v = rhov / rho;
	double e = calecon(rho, v, E);
	pri[0] = rho;
	pri[1] = v;
	pri[2] = (gama - 1) * rho * e;
	return pri;
}

std::vector<double> pri2con(double rho, double v, double p, double gama){
	// convert primitive (rho v p) to conservative  (rho rhov E)
	std::vector<double> con;
	con.resize(3);
	double e = calepri(p, rho, gama);
	// rho rhov E
	con[0] = rho;
	con[1] = rho * v;
	con[2] = rho * e + 0.5 * rho * v * v;
	return con;
}

std::vector<double> Fri_flux(std::vector<double> con0, std::vector<double> con1, double p0, double p1, double dt, double dx){
	// calculate the flux for lax-friedrich
	double rho0 = con0[0];
	double rhov0 = con0[1];
	double E0 = con0[2];
	double v0 = rhov0 / rho0;
	double rho1 = con1[0];
	double rhov1 = con1[1];
	double E1 = con1[2];
	double v1 = rhov1 / rho1;
	std::vector<double> Fri;
	Fri.resize(3);
	std::vector<double> f0;
	f0.resize(3);
	std::vector<double> f1;
	f1.resize(3);
	f0 = f(rho0, v0, p0, E0);
	f1 = f(rho1, v1, p1, E1);
	for(int i = 0; i < 3; i++){
		Fri[i] = 0.5 * (dx/dt) * (con0[i] - con1[i]) + 0.5 * (f1[i] + f0[i]);
	}
	return Fri;
}

std::vector<double> RI_flux(std::vector<double> con0, std::vector<double> con1, double p0, double p1, double dt, double dx, double gama){
	double rho0 = con0[0];
	double rhov0 = con0[1];
	double E0 = con0[2];
	double v0 = rhov0 / rho0;
	double rho1 = con1[0];
	double rhov1 = con1[1];
	double E1 = con1[2];
	double v1 = rhov1 / rho1;
	std::vector<double> u;
	u.resize(3);
	std::vector<double> f0;
	f0.resize(3);
	std::vector<double> f1;
	f1.resize(3);
	f0 = f(rho0, v0, p0, E0);
	f1 = f(rho1, v1, p1, E1);
	for(int i = 0; i < 3; i++){
		u[i] = 0.5 * (con0[i] + con1[i]) - 0.5 * (dt/dx) * (f1[i] - f0[i]);
	}
	double rho = u[0];
	double rhov = u[1];
	double E = u[2];
	double v = rhov / rho;
	
	std::vector<double> pri = con2pri(rho, rhov, E, gama);
	double p = pri[2];
	std::vector<double> RI = f(rho, v, p, E);
	return RI;
}

std::vector<double> FORCE_flux(std::vector<double> con0, std::vector<double> con1, double p0, double p1, double dt, double dx, double gama){
	std::vector<double> flux;
	flux.resize(3);
	std::vector<double> RI = RI_flux(con0, con1, p0, p1, dt, dx, gama);
	std::vector<double> Fri = Fri_flux(con0, con1, p0, p1, dt, dx);
	for(int i = 0; i < 3 ; i++){
		flux[i] = 0.5 * (Fri[i] + RI[i]);
	}
	return flux;
}


double calecon(double rho, double v, double E){
	double e = (E - 0.5 * rho * v * v)/rho;
	return e;
}

double calepri(double p, double rho, double gama){
	double e = p/((gama - 1) * rho);
	return e;
}

double calcs(double p, double rho, double gama){
	double cs = sqrt(gama * p / rho);
	return cs;
}

using namespace std;
int main(int argc, char *argv[]){
	double x0 = 0;
	double x1 = 1;
	double gama = 1.4;
	double c = 0.8;
	double nCells = 200;
	double nPoints = nCells + 1;
	double time = 0.15;
	double t = 0;
	double dx = (x1 - x0)/(nPoints-1);
	std::vector<double> rho;
	std::vector<double> v;
	std::vector<double> p;
	rho.resize(nCells + 2);
	v.resize(nCells + 2);
	p.resize(nCells + 2);
	
	std::vector<double> rhodelta;
	rhodelta.resize(nCells + 2);
	std::vector<double> rhodelta0;
	rhodelta0.resize(nCells + 2);
	std::vector<double> rhodelta1;
	rhodelta1.resize(nCells + 2);
	std::vector<double> rhouL; 
	rhouL.resize(nCells + 2);
	std::vector<double> rhouR; 
	rhouR.resize(nCells + 2);
	double rhow = 0;
	
	std::vector<double> rhovdelta;
	rhovdelta.resize(nCells + 2);
	std::vector<double> rhovdelta0;
	rhovdelta0.resize(nCells + 2);
	std::vector<double> rhovdelta1;
	rhovdelta1.resize(nCells + 2);
	std::vector<double> rhovuL; 
	rhovuL.resize(nCells + 2);
	std::vector<double> rhovuR; 
	rhovuR.resize(nCells + 2);
	double rhovw = 0;
	
	std::vector<double> Edelta;
	Edelta.resize(nCells + 2);
	std::vector<double> Edelta0;
	Edelta0.resize(nCells + 2);
	std::vector<double> Edelta1;
	Edelta1.resize(nCells + 2);
	std::vector<double> EuL; 
	EuL.resize(nCells + 2);
	std::vector<double> EuR; 
	EuR.resize(nCells + 2);
	double Ew = 0;
	
	for(int i = 1; i < nCells + 1; i++){
		double x = x0 + (i - 1) * dx;
		if(x < 0.5){
			rho[i] = 1;
			v[i] = -2;
			p[i] = 0.4;
		}else if(x >= 0.5){
			rho[i] = 1.0;
			v[i] = 2.0;
			p[i] = 0.4;
		}
	}
	
	std::vector<double> rhov;
	rhov.resize(nCells + 2);
	std::vector<double> E;
	E.resize(nCells + 2);
	std::vector<double> inE;
	inE.resize(nCells + 2);
	
	do{
		rho[0] = rho[1];
		rho[nCells + 1] = rho[nCells];
		v[0] = v[1];
		v[nCells + 1] = v[nCells];
		p[0] = p[1];
		p[nCells + 1] = p[nCells];
		inE[0] = inE[1];
		inE[nCells + 1] = inE[nCells];
		
		for(int i = 0; i < nCells + 2; i++){
			std::vector <double> conser = pri2con(rho[i], v[i], p[i], gama);
			rhov[i] = conser[1];
			E[i] = conser[2];
		}
		
		for(int i = 1; i < nCells + 1; i++){
			rhodelta0[i] = rho[i] - rho[i - 1];
			rhodelta1[i] = rho[i + 1] - rho[i];
			rhovdelta0[i] = rhov[i] - rhov[i - 1];
			rhovdelta1[i] = rhov[i + 1] - rhov[i];
			Edelta0[i] = E[i] - E[i - 1];
			Edelta1[i] = E[i + 1] - E[i];
		}
		
		rhodelta0[0] = rhodelta0[1];
		rhodelta0[nCells + 1] = rhodelta0[nCells];
		rhodelta1[0] = rhodelta1[1];
		rhodelta1[nCells + 1] = rhodelta1[nCells];
		rhovdelta0[0] = rhovdelta0[1];
		rhovdelta0[nCells + 1] = rhovdelta0[nCells];
		rhovdelta1[0] = rhovdelta1[1];
		rhovdelta1[nCells + 1] = rhovdelta1[nCells];
		Edelta0[0] = Edelta0[1];
		Edelta0[nCells + 1] = Edelta0[nCells];
		Edelta1[0] = Edelta1[1];
		Edelta1[nCells + 1] = Edelta1[nCells];

		for(int i = 1; i < nCells + 1; i++){
			rhodelta[i] = 0.5 * (1 + rhow) * rhodelta0[i] + 0.5 * (1 - rhow) * rhodelta1[i];
			rhovdelta[i] = 0.5 * (1 + rhovw) * rhovdelta0[i] + 0.5 * (1 - rhovw) * rhovdelta1[i];
			Edelta[i] = 0.5 * (1 + Ew) * Edelta0[i] + 0.5 * (1 - Ew) * Edelta1[i];
		}
		
		rhodelta[0] = rhodelta[1];
		rhodelta[nCells + 1] = rhodelta[nCells];
		rhovdelta[0] = rhovdelta[1];
		rhovdelta[nCells + 1] = rhovdelta[nCells];
		Edelta[0] = Edelta[1];
		Edelta[nCells + 1] = Edelta[nCells];
		
		std::vector<double> rhoL, rhoR;
		rhoL.resize(nCells + 2); 
		rhoR.resize(nCells + 2);
		std::vector<double> rhovL, rhovR;
		rhovL.resize(nCells + 2); 
		rhovR.resize(nCells + 2);
		std::vector<double> EL, ER;
		EL.resize(nCells + 2); 
		ER.resize(nCells + 2);
		
		for(int i = 1; i < nCells + 1; i++){
			double rhor = (rho[i] - rho[i - 1])/(rho[i + 1] - rho[i]);
			if(rho[i + 1] - rho[i] == 0){
				rhor = 0;
			}
			rhoL[i] = rho[i] - 0.5 * Minbee(rhor) * rhodelta[i];
			rhoR[i] = rho[i] + 0.5 * Minbee(rhor) * rhodelta[i];
			
			
			double rhovr = (rhov[i] - rhov[i - 1])/(rhov[i + 1] - rhov[i]);
			if(rhov[i + 1] - rhov[i] == 0){
				rhovr = 0;
			}
			
			rhovL[i] = rhov[i] - 0.5 * Minbee(rhovr) * rhovdelta[i];
			rhovR[i] = rhov[i] + 0.5 * Minbee(rhovr) * rhovdelta[i];
			
			double Er = (E[i] - E[i - 1])/(E[i + 1] - E[i]);
			if(E[i + 1] - E[i] == 0){
				Er = 0;
			}
			EL[i] = E[i] - 0.5 * Minbee(Er) * Edelta[i];
			ER[i] = E[i] + 0.5 * Minbee(Er) * Edelta[i];
		}
		
		rhoL[0] = rhoL[1];
		rhoL[nCells + 1] = rhoL[nCells];
		rhovL[0] = rhovL[1];
		rhovL[nCells + 1] = rhovL[nCells];
		EL[0] = EL[1];
		EL[nCells + 1] = EL[nCells];
		
		rhoR[0] = rhoR[1];
		rhoR[nCells + 1] = rhoR[nCells];
		rhovR[0] = rhovR[1];
		rhovR[nCells + 1] = rhovR[nCells];
		ER[0] = ER[1];
		ER[nCells + 1] = ER[nCells];
		
		double a_max;
		for(int i = 0; i < nCells + 2; i++){
			if(i == 0){
				a_max = fabs(v[i] + calcs(p[i], rho[i], gama));
			}else if (fabs(a_max) < fabs(v[i] + calcs(p[i], rho[i], gama))){
				a_max = fabs(v[i] + calcs(p[i], rho[i], gama));
			}
		}
		
		double dt = c * dx / a_max;
		t = t + dt;
		
		std::vector<double> rhoPlus;
		rhoPlus.resize(nCells + 2);
		std::vector<double> rhovPlus;
		rhovPlus.resize(nCells + 2);
		std::vector<double> EPlus;
		EPlus.resize(nCells + 2);
		std::vector<double> rhoflux; 
		rhoflux.resize(nCells + 1);
		std::vector<double> rhovflux; 
		rhovflux.resize(nCells + 1);
		std::vector<double> Eflux; 
		Eflux.resize(nCells + 1);
		
		std::vector<double> rhoRhalf;
		rhoRhalf.resize(nCells + 2);
		std::vector<double> rhoLhalf;
		rhoLhalf.resize(nCells + 2);
		std::vector<double> rhovRhalf;
		rhovRhalf.resize(nCells + 2);
		std::vector<double> rhovLhalf;
		rhovLhalf.resize(nCells + 2);
		std::vector<double> ERhalf;
		ERhalf.resize(nCells + 2);
		std::vector<double> ELhalf;
		ELhalf.resize(nCells + 2);

		for(int i = 1; i < nCells + 1; i++){
			double vL = rhovL[i]/rhoL[i];
			double vR = rhovR[i]/rhoR[i];
			std::vector<double> priL= con2pri(rhoL[i], rhovL[i], EL[i], gama);
			std::vector<double> priR= con2pri(rhoR[i], rhovR[i], ER[i], gama);
			double pL = priL[2];
			double pR = priR[2];
			std::vector<double> fL = f(rhoL[i],vL,pL,EL[i]);
			std::vector<double> fR = f(rhoR[i],vR,pR,ER[i]);
			
			rhoLhalf[i] = rhoL[i] - 0.5 * (dt/dx) * (fR[0] - fL[0]);
			rhovLhalf[i] = rhovL[i] - 0.5 * (dt/dx) * (fR[1] - fL[1]);
			ELhalf[i] = EL[i] - 0.5 * (dt/dx) * (fR[2] - fL[2]);
			
			rhoRhalf[i] = rhoR[i] - 0.5 * (dt/dx) * (fR[0] - fL[0]);
			rhovRhalf[i] = rhovR[i] - 0.5 * (dt/dx) * (fR[1] - fL[1]);
			ERhalf[i] = ER[i] - 0.5 * (dt/dx) * (fR[2] - fL[2]);
		}
		
		rhoLhalf[0] = rhoLhalf[1];
		rhoLhalf[nCells + 1] = rhoLhalf[nCells];
		rhoRhalf[0] = rhoRhalf[1];
		rhoRhalf[nCells + 1] = rhoRhalf[nCells];
		rhovLhalf[0] = rhovLhalf[1];
		rhovLhalf[nCells + 1] = rhovLhalf[nCells];
		rhovRhalf[0] = rhovRhalf[1];
		rhovRhalf[nCells + 1] = rhovRhalf[nCells];
		ELhalf[0] = ELhalf[1];
		ELhalf[nCells + 1] = ELhalf[nCells];
		ERhalf[0] = ERhalf[1];
		ERhalf[nCells + 1] = ERhalf[nCells];
		
		for(int i = 0; i < nCells + 1; i++){
			double rho0 = rhoRhalf[i];
			double rhov0 = rhovRhalf[i];
			double E0 = ERhalf[i];
			std::vector<double> pri0 = con2pri(rho0, rhov0, E0, gama);
			double p0 = pri0[2];
			
			std::vector<double> con0;
			con0.resize(3);
			con0[0] = rho0;
			con0[1] = rhov0;
			con0[2] = E0;
			
			double rho1 = rhoLhalf[i + 1];
			double rhov1 = rhovLhalf[i + 1];
			double E1 = ELhalf[i + 1];
			std::vector<double> pri1 = con2pri(rho1, rhov1, E1, gama);
			double p1 = pri1[2];
			std::vector<double> con1;
			con1.resize(3);
			con1[0] = rho1;
			con1[1] = rhov1;
			con1[2] = E1;

			std::vector<double> flux = FORCE_flux(con0, con1, p0, p1, dt, dx, gama);
			
			rhoflux[i] = flux[0];
			rhovflux[i] = flux[1];
			Eflux[i] = flux[2];
		}
		
		for(int i = 1; i < nCells + 1; i++){
			rhoPlus[i] = rho[i] - (dt/dx) * (rhoflux[i] - rhoflux[i - 1]);
			rhovPlus[i] = rhov[i] - (dt/dx) * (rhovflux[i] - rhovflux[i - 1]);
			EPlus[i] = E[i] - (dt/dx) * (Eflux[i] - Eflux[i - 1]);
		}
		rho = rhoPlus;
		rhov = rhovPlus;
		E = EPlus;
		
		for(int i = 1; i < nCells + 1; i++){
			std::vector<double> prim = con2pri(rho[i], rhov[i], E[i], gama);
			v[i] = prim[1];
			p[i] = prim[2];
			inE[i] = calecon(rho[i], v[i], E[i]);
		}
	}while(t < time);
	
	double select;
	std::cin >> select;
	if(select == 1){
		ofstream outFile("/Users/chenzefeng/Desktop/SLIC-rho.dat");
		for(int i = 1; i < nCells + 1; i++){
			outFile << x0 + dx * (i - 0.5) << " " << rho[i] << endl;
		}
		outFile.close();
	}else if(select == 2){
		ofstream outFile("/Users/chenzefeng/Desktop/SLIC-v.dat");
		for(int i = 1; i < nCells + 1; i++){
			outFile << x0 + dx * (i - 0.5) << " " << v[i] << endl;
		}
		outFile.close();
	}else if(select == 3){
		ofstream outFile("/Users/chenzefeng/Desktop/SLIC-E.dat");
		for(int i = 1; i < nCells + 1; i++){
			outFile << x0 + dx * (i - 0.5) << " " << inE[i] << endl;
		}
		outFile.close();
	}else if(select == 4){
		ofstream outFile("/Users/chenzefeng/Desktop/SLIC-P.dat");
		for(int i = 1; i < nCells + 1; i++){
			outFile << x0 + dx * (i - 0.5) << " " << p[i] << endl;
		}
		outFile.close();
	}
	
}