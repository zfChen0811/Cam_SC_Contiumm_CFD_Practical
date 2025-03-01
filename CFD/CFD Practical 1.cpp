#include <iostream>
#include <vector> 
#include <cmath>
#include <fstream>

using namespace std;
int main(int argc, char *argv[]) {
	int nPoints = 100; // Chosen such number of points
	double x0 = 0.0; // x minimum value
	double x1 = 1.0; // x max value
	double tStart = 0.0; // t start
	double tStop = 1.0; // t end
	double a = 1.0; // linear advection
	double dx = (x1 - x0) / nPoints; // calculate delta x
	// See problem sheet for more options here:
	double dt = dx; // dt = dx
	//define a vector with enough space to store discretised data
	// at each point
//	std::vector<double> u; 
//	u.resize(nPoints);   // initialize a vector
//	// calculate the	 boundary condition i.e. t = 0
//	for(int i = 0; i < u.size(); i++) {
//		double x = x0 + i * dx; 
//		u[i] = sin(x); // u0(x) = sin(x) 
//	}
//	double t = tStart;
//	do {t = t + dt;
//		//You may want to manually reduce dt if this
//		// would overshoot tStop
//		// Update the data
//		// Insert some code here...
//		std::vector<double> uPlus1; 
//		u.resize(nPoints);
//		//Update the data
//		for(int i = 0; i < u.size(); i++) {
//			// Use the forward difference defined earlier
//			uPlus1[i] = u[i] - a * (dt/dx) * (u[i+1] - u[i]);
//		}
//		// Now replace u with the updated data for the next time step
//		u = uPlus1;
//		//define a vector with enough space to store discretised data
//		// at each point
//	} while (t < tStop);
//	
//
//	
	std::vector<double> u1;
	std::vector<double> u2;
	std::vector<double> u3;
	// Allow for one additional point at each end of the domain
	u1.resize(nPoints+2);
	u2.resize(nPoints+2);
	u3.resize(nPoints+2);
	for(int i = 0; i < u1.size(); i++) {
		// x 0 is at point i=1
		double x = x0 + (i-1) * dx; 
		u1[i] = sin(x); // u0(x) = sin(x) 
		u2[i] = sin(x);
		u3[i] = sin(x);
	}
	// 注意：u[0]此时为x0 - deltax。u[nPoints]此时为x1 + deltax。都超出范围。
	double t = tStart;
	do {
		t = t + dt;
		//You may want to manually reduce dt if this
		// would overshoot tStop
		//Periodic boundary conditions
		u1[0] = u1[nPoints];
		u1[nPoints+1] = u1[1];
		u2[0] = u2[nPoints];
		u2[nPoints+1] = u2[1];
		u3[0] = u3[nPoints];
		u3[nPoints+1] = u3[1];
		// 第一个和最后一个都是不在范围内的，用周期重新赋值
		std::vector<double> uPlus1; 
		uPlus1.resize(nPoints+2);
		std::vector<double> uPlus2; 
		uPlus2.resize(nPoints+2);
		std::vector<double> uPlus3; 
		uPlus3.resize(nPoints+2);
		//Update the data
		//Make sure the limits on this loop correspond only to the
		// true domain
		for(int i = 1; i < nPoints+1; i++){
			// Use the forward difference defined earlier
			uPlus1[i] = u1[i] - a * (dt/dx) * (u1[i+1] - u1[i]);
		}
		// backward
		for(int i = 1; i < nPoints+1; i++){
			uPlus2[i] = u2[i] - a * (dt/dx) * (u2[i] - u2[i-1]);
		}
		// central difference
		for(int i = 1; i < nPoints+1; i++){
			uPlus3[i] = u3[i] - a * (dt/(2*dx)) * (u3[i+1] - u3[i-1]);
		}
		// Now replace u with the updated data for the next time step
		u1 = uPlus1;
		u2 = uPlus2;
		u3 = uPlus3;
	} while (t < tStop);
	
//	ofstream outFile("/Users/chenzefeng/Desktop/advectionsResults_1.dat");
//	for(int i = 1; i < nPoints+1;i++){
//		outFile << x0 + dx * (i - 1) << " " << u1[i] << endl;
//	}
//	outFile.close();

//	ofstream outFile("/Users/chenzefeng/Desktop/advectionsResults_2.dat");
//	for(int i = 1; i < nPoints+1;i++){
//		outFile << x0 + dx * (i - 1) << " " << u2[i] << endl;
//	}
//	outFile.close();
	ofstream outFile("/Users/chenzefeng/Desktop/advectionsResults_3.dat");
	for(int i = 1; i < nPoints+1;i++){
		outFile << x0 + dx * (i - 1) << " " << u3[i] << endl;
	}
	outFile.close();
}