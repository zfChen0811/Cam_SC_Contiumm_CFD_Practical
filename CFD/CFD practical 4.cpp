#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>
#include <fstream>

std::vector<double> f(double rho, double v, double p, double E);
std::vector<double> con2pri(double rho, double rhov, double E, double gama);
std::vector<double> pri2con(double rho, double v, double p, double gama);
double FORCE(double dx, double dt, double rho0, double rho1, double rhov0, double rhov1, double E0, double E1);
std::vector<double> Fri_flux(std::vector<double> con0, std::vector<double> con1, double p0, double p1, double dt, double dx);
std::vector<double> RI_flux(std::vector<double> con0, std::vector<double> con1, double p0, double p1, double dt, double dx, double gama);
std::vector<double> FORCE_flux(std::vector<double> con0, std::vector<double> con1, double p0, double p1, double dt, double dx, double gama);
double calecon(double rho, double v, double E);
double calepri(double p, double rho, double gama);
double calcs(double p, double rho, double gama);


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
int main(int argc, char *argv[]) {
	double x0 = 0;
	double x1 = 1;
	double gama = 1.4;
	double c = 0.9;
	double nCells = 100;
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
	
	for(int i = 1; i < nCells + 1; i++){
		double x = x0 + (i - 1) * dx;
		if(x < 0.5){
			rho[i] = 1;
			v[i] = -2.0;
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
		
		for(int i = 1; i < nCells + 1; i++){
			std::vector <double> conser = pri2con(rho[i], v[i], p[i], gama);
			rhov[i] = conser[1];
			E[i] = conser[2];
		}
		rhov[0] = rhov[1];
		rhov[nCells + 1] = rhov[nCells];
		E[0] = E[1];
		E[nCells + 1] = E[nCells];
		
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

		for(int i = 0; i < nCells + 1; i++){
			std::vector<double> con0;
			con0.resize(3);
			con0[0] = rho[i];
			con0[1] = rhov[i];
			con0[2] = E[i];
			std::vector<double> con1;
			con1.resize(3);
			con1[0] = rho[i + 1];
			con1[1] = rhov[i + 1];
			con1[2] = E[i + 1];
			double p0 = p[i];
			double p1 = p[i + 1];
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
	
	ofstream outFile("/Users/chenzefeng/Desktop/u.dat");
	for(int i = 1; i < nCells + 1; i++){
		outFile << x0 + dx * (i - 0.5) << " " << p[i] << endl;
	}
	
	outFile.close();
}