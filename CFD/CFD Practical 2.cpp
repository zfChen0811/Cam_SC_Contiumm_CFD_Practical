#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
// upwind method
// Lax-Friedrichs method
// Lax-Wendroff method
// Warming-Beam method

using namespace std;
int main(int argc, char *argv[]) {
	int nPoints, a = 1, t_stop = 1, initial_data, method;
	nPoints = 1000;
	std::cin >> method;
	double x0 = -1.0; // x minimum value
	double x1 = 1.0; // x max value
	double t_0 = 0;
	double dx;
	double dt;
	double c =  0.8;
	dx = (x1 - x0)/(nPoints - 1);
	dt = (dx * c) / a;
	
	std::vector<double> x;
	std::vector<double> t;
	std::vector<double> up_u;
	std::vector<double> fri_u;
	std::vector<double> wen_u;
	std::vector<double> beam_u;
	up_u.resize(nPoints + 2);
	fri_u.resize(nPoints + 2);
	wen_u.resize(nPoints + 2);
	beam_u.resize(nPoints + 4);
	std::vector<double> u;
	
	for(int i = 1; i < nPoints + 1 ;i++){
		double x = x0 + (i-1) * dx;
		up_u[i] = exp(-8 * x * x);
		fri_u[i] = exp(-8 * x * x);
		wen_u[i] = exp(-8 * x * x);
		beam_u[i + 1] = exp(-8 * x * x);
	}
	
	for(double t = t_0; t < t_stop; t = t + dt){
		up_u[0] = up_u[nPoints];
		up_u[nPoints + 1] = up_u[1];
		fri_u[0] = fri_u[nPoints];
		fri_u[nPoints + 1] = fri_u[1];
		wen_u[0] = wen_u[nPoints];
		wen_u[nPoints + 1] = wen_u[1];
		beam_u[0] = beam_u[nPoints];
		beam_u[1] = beam_u[nPoints + 1];
		beam_u[nPoints + 2] = beam_u[2];
		beam_u[nPoints + 3] = beam_u[3];
		
		std::vector<double> up_u_new;
		std::vector<double> fri_u_new;
		std::vector<double> wen_u_new;
		std::vector<double> beam_u_new;
		
		up_u_new.resize(nPoints + 2);
		fri_u_new.resize(nPoints + 2);
		wen_u_new.resize(nPoints + 2);
		beam_u_new.resize(nPoints + 4);
		
		if(method == 1){
			if(a >= 0){
				for(int i = 1; i < nPoints + 1; i++){
					up_u_new[i] = up_u[i] - a * (dt/dx) * (up_u[i] - up_u[i - 1]);
				}
			}else{
				for(int i = 1; i < nPoints+1; i++){
					up_u_new[i] = up_u[i] - a * (dt/dx) * (up_u[i + 1] - up_u[i]);
				}
			}
			up_u = up_u_new;
			u = up_u;
		}else if(method == 2){
			for(int i = 1; i < nPoints + 1; i++){
				fri_u_new[i] = 0.5 * (1 + c) * fri_u[i - 1] + 0.5 * (1 - c) * fri_u[i + 1];
			}
			fri_u = fri_u_new;
			u = fri_u;
		}else if(method == 3){
			for(int i = 1; i < nPoints + 1; i++){
				wen_u_new[i] = 0.5 * c * (1 + c) * wen_u[i - 1] + (1 - c * c) * wen_u[i] - 0.5 * c * (1 - c) * wen_u[i + 1];
			}
			wen_u = wen_u_new;
			u = wen_u;
		}else if(method == 4){
			for(int i = 2; i < nPoints + 2; i++){
				beam_u_new[i] = -0.5 * c * (1 - c) * beam_u[i - 2] + c * (2 - c) * beam_u[i - 1] + 0.5 * (c - 1) * (c - 2) * beam_u[i];
			}
			beam_u = beam_u_new;
			u = beam_u;
		}else{
			std::cout << "Enter an error method";
			break;
		}
	}
	
	ofstream outFile("/Users/chenzefeng/Desktop/u.dat");
	for(int i = 1; i < nPoints+1; i++){
		outFile << x0 + dx * (i - 1) << " " << u[i] << endl;
	}
	outFile.close();
}