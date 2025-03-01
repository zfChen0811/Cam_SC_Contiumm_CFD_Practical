#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;

double f(double x){
	return (0.5 * x * x);
}

double f_Fri(double dx, double dt, double u_1, double u_0){
	double value;
	value = 0.5 * (dx/dt) * (u_0 - u_1) + 0.5 * (f(u_1) + f(u_0));
	return value;
}

double f_RI(double u){
	return f(u);
}

double halftime(double dx, double dt, double u_1, double u_0){
	double value;
	value = 0.5 * (u_0 + u_1) - 0.5 * (dt/dx) * (f(u_1) - f(u_0));
	return value;
}

double f_FORCE(double dx, double dt, double u_1, double u_0){
	double value;
	value = 0.5 * (f_Fri(dx, dt, u_1, u_0) + f_RI(halftime(dx, dt, u_1, u_0)));
	return value;
}

double f_GOD(double u_1, double u_0){
	double value, s;
	if(u_0 > u_1){
		s = 0.5 * (u_0 + u_1);
		if(s > 0){
			value = u_0;
		}else{
			value = u_1;
		}
	}else{
		if(u_0 > 0){
			value = u_0;
		}else if(u_0 <= 0 && u_1 >= 0){
			value = 0;
		}else{
			value = u_1;
		}
	}
	return f(value);
}

double calculate_dt(double C, double a_max, double dx){
	double dt;
	dt = C * dx / a_max;
	return dt;
}

int main(int argc, char *argv[]) {
	double nPoints, dx, method; 
	std::cin >> method;
	double x_0 = 0;
	double x_1 = 1.5;
	double nCells = 75;
	double tStop = 0.5;
	double C = 0.8;
	std::vector<double> u_Fri;
	std::vector<double> u_FOR;
	std::vector<double> u_God;
	std::vector<double> flux_Fri;
	std::vector<double> flux_FOR;
	std::vector<double> flux_God;
	
	nPoints = nCells + 1;
	
	u_Fri.resize(nCells + 2);
	u_FOR.resize(nCells + 2);
	u_God.resize(nCells + 2);
	
	flux_Fri.resize(nCells + 1);
	flux_FOR.resize(nCells + 1);
	flux_God.resize(nCells + 1);
	
	dx = (x_1 - x_0)/(nPoints - 1);
		
	for(int i = 1; i < nCells + 1; i++){
		double x = x_0 + (i - 0.5) * dx;
		if(x <= 0.5){
			u_Fri[i] = -0.5;
			u_FOR[i] = -0.5;
			u_God[i] = -0.5;
		}else if(x <= 1){
			u_Fri[i] = 1;
			u_FOR[i] = 1;
			u_God[i] = 1;
		}else{
			u_Fri[i] = 0;
			u_FOR[i] = 0;
			u_God[i] = 0;
		}
	}
//	for(int i = 1; i < nCells + 1; i++){
//		double x = x_0 + (i - 0.5) * dx;
//		if(x <= 0.5){
//			u_Fri[i] = 2;
//			u_FOR[i] = 2;
//			u_God[i] = 2;
//		}else if(x <= 1){
//			u_Fri[i] = 1;
//			u_FOR[i] = 1;
//			u_God[i] = 1;
//		}
//	}
//	
//	for(int i = 1; i < nCells + 1; i++){
//		double x = x_0 + (i - 0.5) * dx;
//		if(x <= 0.5){
//			u_Fri[i] = 1;
//			u_FOR[i] = 1;
//			u_God[i] = 1;
//		}else if(x <= 1){
//			u_Fri[i] = 2;
//			u_FOR[i] = 2;
//			u_God[i] = 2;
//		}
//	}

	

	double t = 0; 
	
	if(method == 1){
		do{
			u_Fri[0] = u_Fri[1];
			u_Fri[nCells + 1] = u_Fri[nCells];
			
			double a_max;
			for(int i = 0; i < nCells + 2; i++){
				if(i == 0){
					a_max = fabs(u_Fri[0]);
				}else if (fabs(a_max) < fabs(u_Fri[i])) {
					a_max = fabs(u_Fri[i]);
				}
			}
			
			double dt = calculate_dt(C, a_max, dx);
			
			t = t + dt;
			
			std::vector<double> u_FriPlus;
			u_FriPlus.resize(nCells + 2);
			
			for(int i = 0; i < nCells + 1; i++){
				flux_Fri[i] = f_Fri(dx, dt, u_Fri[i + 1], u_Fri[i]);
			}
			
			for(int i = 1; i < nCells + 1; i++){
				u_FriPlus[i] = u_Fri[i] - (dt/dx) * (flux_Fri[i] - flux_Fri[i - 1]);
			}
			u_Fri = u_FriPlus;
		}while(t < tStop);
	}
	
	if(method == 2){
		do{
			u_FOR[0] = u_FOR[1];
			u_FOR[nCells + 1] = u_FOR[nCells];
			double a_max;
			
			for(int i = 0; i < nCells + 2; i++){
				if(i == 0){
					a_max = fabs(u_FOR[0]);
				}else if (fabs(a_max) < fabs(u_FOR[i])) {
					a_max = fabs(u_FOR[i]);
				}
			}
			
			double dt = calculate_dt(C, a_max, dx);
			t = t + dt;
			
			std::vector<double> u_FORPlus;
			u_FORPlus.resize(nCells + 2);
			
			for(int i = 0; i < nCells + 1; i++){
				flux_FOR[i] = f_FORCE(dx, dt, u_FOR[i + 1], u_FOR[i]);
			}
			
			for(int i = 1; i < nCells + 1; i++){
				u_FORPlus[i] = u_FOR[i] - (dt/dx) * (flux_FOR[i] - flux_FOR[i - 1]);
			}
			u_FOR = u_FORPlus;
		}while(t < tStop);
	}
	
	if(method == 3){
		do{
			u_God[0] = u_God[1];
			u_God[nCells + 1] = u_God[nCells];
			
			double a_max;
			
			for(int i = 0; i < nCells + 2; i++){
				if(i == 0){
					a_max = fabs(u_God[0]);
				}else if (fabs(a_max) < fabs(u_God[i])){
					a_max = fabs(u_God[i]);
				}
			}
			
			double dt = calculate_dt(C, a_max, dx);
			t = t + dt;
			
			std::vector<double> u_GodPlus;
			u_GodPlus.resize(nCells + 2);
			
			
			for(int i = 0; i < nCells + 1; i++){
				flux_God[i] = f_GOD(u_God[i + 1], u_God[i]);
			}
			
			for(int i = 1; i < nCells + 1; i++){
				u_GodPlus[i] = u_God[i] - (dt/dx) * (flux_God[i] - flux_God[i - 1]);
			}
			u_God = u_GodPlus;
			
		}while(t < tStop);
	}
	
	std::vector<double> u;
	u.resize(nCells + 2);
	if(method == 1){
		u_Fri[0] = u_Fri[1];
		u_Fri[nCells + 1] = u_Fri[nCells];
		u = u_Fri;
	}else if(method == 2){
		u_FOR[0] = u_FOR[1];
		u_FOR[nCells + 1] = u_FOR[nCells];
		u = u_FOR;
	}else if(method == 3){
		u_God[0] = u_God[1];
		u_God[nCells + 1] = u_God[nCells];
		u = u_God;
	}
	
	ofstream outFile("/Users/chenzefeng/Desktop/u.dat");
	for(int i = 1; i < nCells + 1; i++){
		outFile << x_0 + (i - 0.5) * dx << " " << u[i] << endl;
	}
	outFile.close();
	
}