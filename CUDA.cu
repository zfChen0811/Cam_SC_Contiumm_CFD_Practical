#include <cuda.h>
#include <thrust/device_vector.h>
#include <thrust/inner_product.h>
#include <iostream>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>
#include <fstream>
using namespace std;

#define GAMMA 1.4
#define C 0.8
#define NUM_VARS 4
int num_var = 4;
double t = 0, stop_time = 0.3, dt;
// double x_0 = 0,x_1 = 1,y_0 = 0,y_1 = 1;
double x_0 = 0,x_1 = 225,y_0 = -44.5,y_1 = 44.5;
int nxCells = 100, nyCells = 40, boundary_cells = 2;
int nxPoints = nxCells + 1, nyPoints = nyCells + 1;
double dx = (x_1 - x_0)/(nxPoints-1);
double dy = (y_1 - y_0)/(nyPoints-1);

enum Processor {CPU, GPU};
struct Grid
{
	double* data;
	int xCells;
	int yCells;
	int BCells;
	double xMin;
	double yMin;
	double xMax;
	double yMax;
	int num_vars;
	Processor proc; 
	// Access to data as l−value
	
	Grid(Processor p){
		xCells = nxCells;
		yCells = nyCells;
		BCells = boundary_cells;
		xMin = 0;
		xMax = 1;
		yMin = 0;
		yMax = 1;
		proc = p;
		num_vars = 4;
		
		switch(proc){
			case CPU:
				data = new double[(xCells + 2 * BCells) * (yCells + 2 * BCells) * NUM_VARS];
				break;
			case GPU:
				cudaMalloc((void **)&data, (xCells + 2 * BCells) * (yCells + 2 * BCells) * NUM_VARS * sizeof(double));
				break; 
		}
	}
	
	Grid(int x_Cells, int y_Cells, Processor p){
		xCells = x_Cells;
		yCells = y_Cells;
		BCells = boundary_cells;
		xMin = 0;
		xMax = 1;
		yMin = 0;
		yMax = 1;
		proc = p;
		
		switch(proc){
			case CPU:
				data = new double[(xCells + 2 * BCells) * (yCells + 2 * BCells) * NUM_VARS];
				break;
			case GPU:
				cudaMalloc((void **)&data, (xCells + 2 * BCells) * (yCells + 2 * BCells) * NUM_VARS * sizeof(double));
				
				break; 
		}
	}
	__device__ __host__
	double operator()(int i, int j, int k)const
	{
		//return data[NUM_VARS * (i + j * (xCells + 2 * BCells)) + k];
		return data[i + j * (xCells + 2 * BCells) + k *  (xCells + 2 * BCells) * (yCells + 2 * BCells)];
	}
	__device__ __host__
	double& operator()(int i, int j, int k){
		//return data[NUM_VARS * (i + j * (xCells + 2 * BCells)) + k];
		return data[i + j * (xCells + 2 * BCells) + k *  (xCells + 2 * BCells) * (yCells + 2 * BCells)];
	}
};

#define CUDA_CHECK {cudaDeviceSynchronize();	\
  cudaError_t err = cudaGetLastError();\
  if(err){\
    std::cout << "Error: " << cudaGetErrorString(err) << " line " << __LINE__ << std::endl; \
    exit(1);\
  }}


__global__ void GPU_calamax(Grid u_pri, double *res);
double Minbee(double r);
void f_con(double *f_value, double *con);
void g_con(double *g_value, double *con);
void con2pri(Grid &pri, Grid &con, const int &i, const int&j);
void pri2con(Grid &con, Grid &pri, const int &i, const int &j);
template<int coord> void Fri_flux(double *Fri, double *con0, double *con1);
template<int coord> void RI_flux(double *RI, double *con0, double *con1);
template<int coord> void FORCE_flux(double *flux,double *con0, double *con1);
double calecon(Grid &con, const int& i, const int& j);
double calepri(Grid &pri, const int &i, const int& j);
void tran_bound(Grid &u);
void upri2ucon(Grid &u_con,Grid &u_pri);
void ucon2upri(Grid &u_pri,Grid &u_con);
void data_recon(Grid &L, Grid &R, Grid &u_con);
template<int coord>void half_time(Grid &halfL, Grid& halfR, Grid &u_con);
template<int coord>void SLIC(Grid &flux, Grid &u_con);
template<int coord>void update(Grid &u_con,Grid &flux);
double caldt(Grid u_pri, int blocks, dim3 dimGrid, dim3 dimBlock);



__global__ void GPU_calamax(Grid u_pri, double *res){
	__shared__ double a[32][32];
	int boundcells = u_pri.BCells;
	int xCells = u_pri.xCells;
	int yCells = u_pri.yCells;
  	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
  	int locI = threadIdx.x;
	int locJ = threadIdx.y;
	int blockSizeX = blockDim.x;
	int blockSizeY = blockDim.y;
	double v_ij = sqrt(u_pri(i, j, 1) * u_pri(i, j, 1) + u_pri(i, j, 2) * u_pri(i, j, 2));
	
	if(i < xCells + 2 * boundcells && j < yCells + 2 * boundcells){
		a[locI][locJ] = 1.0*(v_ij + sqrt(GAMMA *  u_pri(i,j,3)/ u_pri(i,j,0)));
	}else{
		a[locI][locJ] = 0.0;
	}
	
	__syncthreads();

	if( blockSizeX >= 32 && blockSizeY >= 32 && locI < 16 && locJ < 16){
    	a[locI][locJ] = fmax(a[locI][locJ], a[locI][locJ + 16]);
		a[locI][locJ] = fmax(a[locI][locJ], a[locI + 16][locJ]);
		a[locI][locJ] = fmax(a[locI][locJ], a[locI + 16][locJ + 16]);
		__syncthreads();
  	}

	if( blockSizeX >= 16 && blockSizeY >= 16 && locI < 8 && locJ < 8){
    	a[locI][locJ] = fmax(a[locI][locJ], a[locI][locJ + 8]);
		a[locI][locJ] = fmax(a[locI][locJ], a[locI + 8][locJ]);
		a[locI][locJ] = fmax(a[locI][locJ], a[locI + 8][locJ + 8]);
		__syncthreads();
  	}

	if( blockSizeX >= 8 && blockSizeY >= 8 && locI < 4 && locJ < 4){
    	a[locI][locJ] = fmax(a[locI][locJ], a[locI][locJ + 4]);
		a[locI][locJ] = fmax(a[locI][locJ], a[locI + 4][locJ]);
		a[locI][locJ] = fmax(a[locI][locJ], a[locI + 4][locJ + 4]);
		__syncthreads();
  	}

	if( blockSizeX >= 4 && blockSizeY >= 4 && locI < 2 && locJ < 2){
    	a[locI][locJ] = fmax(a[locI][locJ], a[locI][locJ + 2]);
		a[locI][locJ] = fmax(a[locI][locJ], a[locI + 2][locJ]);
		a[locI][locJ] = fmax(a[locI][locJ], a[locI + 2][locJ + 2]);
		__syncthreads();
  	}

	if( blockSizeX >= 2 && blockSizeY >= 2 && locI < 1 && locJ < 1){
    	a[locI][locJ] = fmax(a[locI][locJ], a[locI][locJ + 1]);
		a[locI][locJ] = fmax(a[locI][locJ], a[locI + 1][locJ]);
		a[locI][locJ] = fmax(a[locI][locJ], a[locI + 1][locJ + 1]);
		__syncthreads();
  	}

	if(locI == 0 && locJ == 0){
    	res[blockIdx.y * gridDim.x + blockIdx.x] = a[locI][locJ];
	}
}

__global__ void GPU_un_calamax(Grid u_pri, double *res){
	__shared__ double a[32][32];
	int boundcells = u_pri.BCells;
	int xCells = u_pri.xCells;
	int yCells = u_pri.yCells;
  	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
  	int locI = threadIdx.x;
	int locJ = threadIdx.y;
	double v_ij = sqrt(u_pri(i, j, 1) * u_pri(i, j, 1) + u_pri(i, j, 2) * u_pri(i, j, 2));
	
	if(i < xCells + 2 * boundcells && j < yCells + 2 * boundcells){
		a[locI][locJ] = 1.0*(v_ij + sqrt(GAMMA *  u_pri(i,j,3)/ u_pri(i,j,0)));
	}else{
		a[locI][locJ] = 0.0;
	}
	
	__syncthreads();
    res[j * blockDim.y + i] = a[locI][locJ];
}

double un_caldt(Grid u_pri, int blocks, dim3 dimGrid, dim3 dimBlock){
  	double *res;
  	cudaMalloc(&res, dimBlock.x * dimBlock.y * dimGrid.x * dimGrid.y * sizeof(double));
  	CUDA_CHECK;
  	GPU_un_calamax<<<dimGrid, dimBlock>>>(u_pri, res);
  	CUDA_CHECK;
  	double* resHost = new double[dimBlock.x * dimBlock.y * dimGrid.x * dimGrid.y];
  	cudaMemcpy(resHost, res, sizeof(double) * dimBlock.x * dimBlock.y * dimGrid.x * dimGrid.y, cudaMemcpyDeviceToHost);
  	CUDA_CHECK;
  	double a_max= -1e6;
  	for(int i = 0 ; i < dimBlock.x * dimBlock.y * dimGrid.x * dimGrid.y ; i++){
   		a_max = fmax(resHost[i],a_max);
  	}
  	delete[] resHost;
  	cudaFree(res);
  	CUDA_CHECK;
	double dt = C * (fmin(dx, dy) / a_max);
  	return dt;
}



__device__ double GPU_Minbee(double r){
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
template<int coord>
__device__ void GPU_f_con(double *f_value, double *con){
	// calculate f
	// based on conservative value roh rohv E

	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double E = con[3];
	double e = (E - 0.5 * rho * (vx * vx + vy * vy))/rho;
	double p = (GAMMA - 1) * con[0] * e;
	f_value[0] =  (coord == 0) ? (rho * vx) : (rho * vy);
	f_value[1] =  (coord == 0) ? (rho * vx * vx + p) : (rho * vx * vy);
	f_value[2] =  (coord == 0) ? (rho * vx * vy) : (rho * vy * vy + p);
	f_value[3] =  (coord == 0) ? ((E + p) * vx) : ((E + p) * vy);
}

__device__ void GPU_g_con(double *g_value, double *con){
	// calculate f
	// based on conservative value roh rohv E
	double rho = con[0];
	double vx = con[1] / con[0];
	double vy = con[2] / con[0];
	double E = con[3];
	double e = (E - 0.5 * rho * (vx * vx + vy * vy))/rho;
	double p = (GAMMA - 1) * con[0] * e;
	g_value[0] = rho * vy;
	g_value[1] = rho * vx * vy;
	g_value[2] = rho * vy * vy + p;
	g_value[3] = (E + p) * vy;
}

__device__ void GPU_con2pri(Grid pri, Grid con, const int &i, const int&j){
	double rho = con(i,j,0);
	double vx = con(i,j,1) / con(i,j,0);
	double vy = con(i,j,2) / con(i,j,0);
	double E = con(i,j,3);
	double e = (E - 0.5 * rho * (vx * vx + vy * vy))/rho;
	pri(i, j, 0) = rho;
	pri(i, j, 1) = vx;
	pri(i, j, 2) = vy;
	pri(i, j, 3) = (GAMMA - 1) * rho * e;
}

__device__ void GPU_pri2con(Grid con, Grid pri, const int &i, const int&j){
	double rho = pri(i,j,0);
	double vx = pri(i,j,1);
	double vy = pri(i,j,2);
	double p = pri(i,j,3);
	double e = p/((GAMMA - 1) * rho);
	con(i, j, 0) = rho;
	con(i, j, 1) = rho * vx;
	con(i, j, 2) = rho * vy;
	con(i, j, 3) = rho * e + 0.5 * rho * (vx * vx + vy * vy);
}
template<int coord>
__device__ void GPU_Fri_flux(double *Fri, double *con0, double *con1, double dx, double dy, double dt){
	// calculate the flux for lax-friedrich
	double f0[NUM_VARS], f1[NUM_VARS];
	GPU_f_con<coord>(f0, con0);
	GPU_f_con<coord>(f1, con1);
	double dz = (coord == 0) ? dx : dy;
	for(int i = 0; i < NUM_VARS; i++){
		Fri[i] = 0.5 * (dz/dt) * (con0[i] - con1[i]) + 0.5 * (f1[i] + f0[i]);
	}
}
template<int coord>
__device__ void GPU_RI_flux(double *RI, double *con0, double *con1,double dx, double dy, double dt){
	double u[NUM_VARS];
	double f0[NUM_VARS],f1[NUM_VARS];
	GPU_f_con<coord>(f0, con0);
	GPU_f_con<coord>(f1, con1);
	double dz = (coord == 0) ? dx : dy;
	for(int i = 0; i < NUM_VARS; i++){
		u[i] = 0.5 * (con0[i] + con1[i]) - 0.5 * (dt/dz) * (f1[i] - f0[i]);
	}
	GPU_f_con<coord>(RI, u);
}
template<int coord>
__device__ void GPU_FORCE_flux(double *flux,double *con0, double *con1, double dx, double dy, double dt){
	double RI[NUM_VARS], Fri[NUM_VARS];
	GPU_RI_flux<coord>(RI, con0, con1, dx,dy,dt);
	GPU_Fri_flux<coord>(Fri, con0, con1,dx,dy,dt);
	for(int i = 0; i < NUM_VARS; i++){
		flux[i] = 0.5 * (Fri[i] + RI[i]);
	}
}

__global__ void GPU_tran_bound_x(Grid u){
	int boundcells = u.BCells;
	int xCells = u.xCells;
	int yCells = u.yCells;
  	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int num_var = u.num_vars;
	if(i >= 0 && i < xCells + boundcells * 2 && j >= 0 && j < boundcells){
		for(int k = 0; k < num_var; k++){
			u(i, j, k) = u(i, boundcells, k);
			u(i,yCells + boundcells + j,k) = u(i, yCells + boundcells - 1, k);
		}
		__syncthreads();
	}
}
__global__ void GPU_tran_bound_y(Grid u){
	int boundcells = u.BCells;
	int xCells = u.xCells;
	int yCells = u.yCells;
  	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int num_var = u.num_vars;
	if(j >= 0 && j < yCells + boundcells * 2 &&  i >= 0 && i < boundcells){
		for(int k = 0; k < num_var; k++){
			u(i, j, k) = u(boundcells, j, k);
			u(xCells + boundcells + i, j, k) = u(xCells + boundcells - 1, j, k);
		}
		__syncthreads();
	}
}

__global__ void GPU_upri2ucon(Grid u_con,Grid u_pri){
	int boundcells = u_con.BCells;
	int xCells = u_con.xCells;
	int yCells = u_con.yCells;
  	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	if(i >= 0 && i < xCells + 2 * boundcells && j >=0 &&  j < yCells + 2 * boundcells){
		GPU_pri2con(u_con, u_pri, i, j);
	}
}

__global__ void GPU_ucon2upri(Grid u_pri,Grid u_con){
	int boundcells = u_con.BCells;
	int xCells = u_con.xCells;
	int yCells = u_con.yCells;
  	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	if(i >= 0 && i < xCells + 2 * boundcells && j >=0 &&  j < yCells + 2 * boundcells){
		GPU_con2pri(u_pri, u_con, i, j);
	}
}


__device__ void GPU_data_recon(double *L, double *R, double *u_conL, double *u_con,double *u_conR){
	double delta0[NUM_VARS], delta1[NUM_VARS],delta[NUM_VARS];
	double w[NUM_VARS];
	double r[NUM_VARS];
	for(int k = 0; k < 4; k++){
		w[k] = 0;
		delta0[k] = u_con[k] - u_conL[k];
		delta1[k] = u_conR[k] - u_con[k];
		delta[k] = 0.5 * (1 + w[k]) * delta0[k] + 0.5 * (1 - w[k]) * delta1[k];
		
		r[k] = (u_con[k] - u_conL[k])/(u_conR[k] - u_con[k]);
		if(u_conR[k] - u_con[k] == 0){
			r[k] = 0;
		}
		L[k] = u_con[k] - 0.5 * GPU_Minbee(r[k]) * delta[k];
		R[k] = u_con[k] + 0.5 * GPU_Minbee(r[k]) * delta[k];
	}
	__syncthreads();
}
template<int coord>
__device__ void GPU_half_time(double *halfL, double *halfR, double *L, double *R, double dx, double dy, double dt){
	double fxL[NUM_VARS],fxR[NUM_VARS];
	GPU_f_con<coord>(fxL, L);
	GPU_f_con<coord>(fxR, R);
	double dz = (coord == 0) ? dx : dy;
	for(int k = 0; k < NUM_VARS; k++){
		halfL[k] = L[k] - 0.5 * (dt/dz) * (fxR[k] - fxL[k]);
		halfR[k] = R[k] - 0.5 * (dt/dz) * (fxR[k] - fxL[k]);
	}
}

template<int coord>
__global__ void GPU_half_time_update(Grid u_con, Grid GPU_halfL, Grid GPU_halfR,double dx, double dy, double dt){
	int boundcells = u_con.BCells;
	int xCells = u_con.xCells;
	int yCells = u_con.yCells;
  	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	double ini_halfL[NUM_VARS], ini_halfR[NUM_VARS];
	double uL[NUM_VARS], uR[NUM_VARS], u[NUM_VARS];
	double L[NUM_VARS], R[NUM_VARS];
	int dimi_0 = boundcells -  ((coord == 0) ? 1 : 0);
	int dimi_1 = xCells + boundcells + ((coord == 0) ? 1 : 0);
	int dimj_0 = boundcells -  ((coord == 1) ? 1 : 0);
	int dimj_1 = yCells + boundcells + ((coord == 1) ? 1 : 0);
    int dimx = ((coord == 0) ? 1 : 0);
	int dimy = ((coord == 1) ? 1 : 0);
	if(i >= dimi_0 && i < dimi_1 && j >= dimj_0 && j < dimj_1){
		for(int k = 0; k < NUM_VARS; k++){
			uL[k] = u_con(i - dimx,j - dimy,k);
			u[k] = u_con(i,j,k);
			uR[k] = u_con(i + dimx,j + dimy,k);
		}
		__syncthreads();
		GPU_data_recon(L,R,uL,u,uR);
		GPU_half_time<coord>(ini_halfL, ini_halfR, L,R,dx,dy,dt);
		for(int k = 0; k < NUM_VARS; k++){
			GPU_halfL(i,j,k) = ini_halfL[k];
			GPU_halfR(i,j,k) = ini_halfR[k];
		}
		__syncthreads();
	}
}

template<int coord>
__global__ void GPU_SLIC(Grid flux, Grid u_con, Grid GPU_halfL, Grid GPU_halfR, double dx, double dy, double dt){
	int boundcells = u_con.BCells;
	int xCells = u_con.xCells;
	int yCells = u_con.yCells;
  	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	double dimi_0 = boundcells - ((coord == 0) ? 1 : 0);
	double dimi_1 = boundcells + xCells;
	double dimj_0 = boundcells - ((coord == 1) ? 1 : 0);
	double dimj_1 = yCells + boundcells;
	double ini_flux[NUM_VARS];
	double ini_halfR[NUM_VARS];
	double ini_halfL[NUM_VARS];
    int dimx = i + ((coord == 0) ? 1 : 0);
	int dimy = j + ((coord == 1) ? 1 : 0);
	if(j >= dimj_0 && j < dimj_1 && i>=dimi_0 && i < dimi_1){
		for(int k = 0;k<NUM_VARS;k++){
			ini_halfR[k] = GPU_halfR(i,j,k);
			ini_halfL[k] = GPU_halfL(dimx,dimy,k);
		}
		__syncthreads();
		GPU_FORCE_flux<coord>(ini_flux, ini_halfR, ini_halfL, dx,dy,dt);
		for(int k = 0; k< NUM_VARS; k++){
			flux(i,j,k) = ini_flux[k];
		}
		__syncthreads();
	}
}

template<int coord>
__global__ void GPU_update(Grid u_con, Grid flux, double dx, double dy, double dt){
	int boundcells = u_con.BCells;
	int xCells = u_con.xCells;
	int yCells = u_con.yCells;
  	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	int dimx = i - ((coord == 0) ? 1 : 0);
	int dimy = j - ((coord == 1) ? 1 : 0);
	double dz = (coord == 0) ? dx : dy;
	if(j >= boundcells && j < yCells + boundcells && i >= boundcells && i < xCells + boundcells){
		for(int k = 0; k < 4; k++){
			u_con(i, j, k) = u_con(i, j, k) - (dt/dz) * (flux(i, j, k) - flux(dimx, dimy, k));
		}
		__syncthreads();
	}
}
template<int coord>
__global__ void Overlap_GPU_SLIC(Grid u_con, Grid GPU_halfL, Grid GPU_halfR,double dx, double dy, double dt){
	__shared__ double u[6][6][NUM_VARS];
	__shared__ double uL[6][6][NUM_VARS];
	__shared__ double uR[6][6][NUM_VARS];
	__shared__ double flux[6][6][NUM_VARS];
	int boundcells = u_con.BCells;
	int xCells = u_con.xCells;
	int yCells = u_con.yCells;
	int i = blockIdx.x * blockDim.x + threadIdx.x - 2 * ((coord == 0) ? 1 : 0) * blockIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y - 2 * ((coord == 1) ? 1 : 0) * blockIdx.y;
	int locI = threadIdx.x;
	int locJ = threadIdx.y;
	int dimi_0_1 =  boundcells - ((coord == 0) ? 1 : 0);
	int dimi_1_1 =  xCells + boundcells + ((coord == 0) ? 1 : 0);
	int dimj_0_1 = boundcells - ((coord == 1) ? 1 : 0);
	int dimj_1_1 = yCells + boundcells + ((coord == 1) ? 1 : 0);
	int dimi_0_2 =  boundcells - ((coord == 0) ? 1 : 0);
	int dimi_1_2 =  xCells + boundcells;
	int dimj_0_2 = boundcells - ((coord == 1) ? 1 : 0);
	int dimj_1_2 = yCells + boundcells;
	int kI = ((coord == 0) ? 1 : 0);
	int kJ = ((coord == 1) ? 1 : 0);
	double dz = ((coord == 0) ? dx : dy);
	if(j >= dimj_0_1 && j < dimj_1_1 && i >= dimi_0_1 && i < dimi_1_1){
		for(int k = 0; k < NUM_VARS;k++){
			u[locI][locJ][k] = u_con(i,j,k);
		}
		__syncthreads();
		for(int k = 0; k < NUM_VARS;k++){
			if((coord == 0 && locI == 0) || (coord == 1 && locJ == 0)){
				uL[locI][locJ][k] = u_con(i - kI,j - kJ,k);
			}else{
				uL[locI][locJ][k] = u[locI - kI][locJ - kJ][k];
			}
			if((coord == 0 && locI == blockDim.x - 1) || (coord == 1 && locJ == blockDim.y - 1)){
				uR[locI][locJ][k] = u_con(i + kI,j + kJ,k);
			}else{
				uR[locI][locJ][k] = u[locI + kI][locJ + kJ][k];
			}
		}
		GPU_data_recon(uL[locI][locJ], uR[locI][locJ],uL[locI][locJ],u[locI][locJ],uR[locI][locJ]);
		GPU_half_time<coord>(uL[locI][locJ], uR[locI][locJ],uL[locI][locJ], uR[locI][locJ], dx,dy,dt);
		__syncthreads();
	}
	if(j >= dimj_0_2 && j < dimj_1_2 && i >= dimi_0_2 && i < dimi_1_2){
		int locK = ((coord == 0) ? locI : locJ);
		if(locK < 5){
			GPU_FORCE_flux<coord>(flux[locI][locJ], uR[locI][locJ], uL[locI + kI][locJ + kJ],dx,dy,dt);
		}
		__syncthreads();
	}
	if(j >= boundcells && j < yCells + boundcells && i >= boundcells && i < xCells + boundcells){
		int locK = ((coord == 0) ? locI : locJ);
		if(locK > 0 && locK < 5){
			for(int k = 0; k< NUM_VARS; k++){
				u_con(i,j,k) = u_con(i, j, k) - (dt/dz) * (flux[locI][locJ][k] - flux[locI - kI][locJ - kJ][k]);
			}
		}
		__syncthreads();
	}

}

double CPU_calcs(Grid &pri, const int &i, const int &j){
	double cs = sqrt(GAMMA *  pri(i,j,3)/ pri(i,j,0));
	return cs;
}

double CPU_calamax(Grid &u_pri){
	double a_max = 1e-16;
	double v_ij;
	for(int i = 0; i < nxCells + 2 * boundary_cells; i++){
		for(int j = 0; j < nyCells + 2 * boundary_cells; j++){
			v_ij = sqrt(u_pri(i, j, 1) * u_pri(i, j, 1) + u_pri(i, j, 2) * u_pri(i, j, 2));
			if(a_max < v_ij + CPU_calcs(u_pri, i, j)){
				a_max = v_ij + CPU_calcs(u_pri, i, j);
			}
		}
	}
	return a_max;
}


double CPU_caldt(Grid &u_pri){
	double a_max = CPU_calamax(u_pri);
	double dt = C * fmin(dx, dy) / a_max;
	return dt;
}


double caldt(Grid u_pri, int blocks, dim3 dimGrid, dim3 dimBlock){
  	double *res;
  	cudaMalloc(&res, blocks * sizeof(double));
  	CUDA_CHECK;
  	GPU_calamax<<<dimGrid, dimBlock>>>(u_pri, res);
  	CUDA_CHECK;
  	double* resHost = new double[blocks];
  	cudaMemcpy(resHost, res, sizeof(double) * blocks, cudaMemcpyDeviceToHost);
  	CUDA_CHECK;
  	double a_max= -1e6;
  	for(int i = 0 ; i < blocks ; i++){
   		a_max = fmax(resHost[i],a_max);
  	}
  	delete[] resHost;
  	cudaFree(res);
  	CUDA_CHECK;
	double dt = C * (fmin(dx, dy) / a_max);
  	return dt;
}

__global__ void GPU_initial_u_pri(Grid u_pri,double x_0,double y_0,double dx,double dy){
	int boundcells = u_pri.BCells;
	int xCells = u_pri.xCells;
	int yCells = u_pri.yCells;
  	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
    double x_center = 35.0;
    double y_center = 0.0;
    double r = 25;
	double x, y;
	// if(i >= boundcells && i < xCells + boundcells && j >= boundcells && j < yCells + boundcells){
	// 	x = x_0 + (i - boundcells + 0.5) * dx;
	// 	y = y_0 + (j - boundcells + 0.5) * dy;
	// 	if(x < 0.5 && y >= 0.5){
	// 		u_pri(i,j,0) = 0.5323;
	// 		u_pri(i,j,1) = 1.206;
	// 		u_pri(i,j,2) = 0;
	// 		u_pri(i,j,3) = 0.3;
	// 	}else if(x < 0.5 && y < 0.5){
	// 		u_pri(i,j,0) = 0.138;
	// 		u_pri(i,j,1) = 1.206;
	// 		u_pri(i,j,2) = 1.206;
	// 		u_pri(i,j,3) = 0.029;
	// 	}else if(x >= 0.5 && y >= 0.5){
	// 		u_pri(i,j,0) = 1.5;
	// 		u_pri(i,j,1) = 0;
	// 		u_pri(i,j,2) = 0;
	// 		u_pri(i,j,3) = 1.5;
	// 	}else if(x >= 0.5 && y < 0.5){
	// 		u_pri(i,j,0) = 0.5323;
	// 		u_pri(i,j,1) = 0;
	// 		u_pri(i,j,2) = 1.206;
	// 		u_pri(i,j,3) = 0.3;
	// 	}
	// }
    if(i >= boundcells && i < xCells + boundcells && j >= boundcells && j < yCells + boundcells){
		x = x_0 + (i - boundcells + 0.5) * dx;
		y = y_0 + (j - boundcells + 0.5) * dy;
		if(x < 5){
			u_pri(i,j,0) = 1.7755;
			u_pri(i,j,1) = 110.63;
			u_pri(i,j,2) = 0.0;
			u_pri(i,j,3) = 159060.0;
		}else if(sqrt((x - x_center) * (x - x_center) + (y - y_center) * (y - y_center)) > r){
			u_pri(i,j,0) = 1.29;
			u_pri(i,j,1) = 0.0;
			u_pri(i,j,2) = 0.0;
			u_pri(i,j,3) = 101325.0;
		}else{
            u_pri(i,j,0) = 0.214;
			u_pri(i,j,1) = 0.0;
			u_pri(i,j,2) = 0.0;
			u_pri(i,j,3) = 101325.0;
        } 
	}
}

void initial_u_pri(Grid u_pri){
	for(int i = boundary_cells;  i < nxCells + boundary_cells; i++){
		for(int j = boundary_cells; j < nyCells + boundary_cells; j++){
			double x = x_0 + (i - boundary_cells + 0.5) * dx;
			double y = y_0 + (j - boundary_cells + 0.5) * dy;
			if(x < 0.5 && y >= 0.5){
				u_pri(i,j,0) = 0.5323;
				u_pri(i,j,1) = 1.206;
				u_pri(i,j,2) = 0;
				u_pri(i,j,3) = 0.3;
			}else if(x < 0.5 && y < 0.5){
				u_pri(i,j,0) = 0.138;
				u_pri(i,j,1) = 1.206;
				u_pri(i,j,2) = 1.206;
				u_pri(i,j,3) = 0.029;
			}else if(x >= 0.5 && y >= 0.5){
				u_pri(i,j,0) = 1.5;
				u_pri(i,j,1) = 0;
				u_pri(i,j,2) = 0;
				u_pri(i,j,3) = 1.5;
			}else if(x >= 0.5 && y < 0.5){
				u_pri(i,j,0) = 0.5323;
				u_pri(i,j,1) = 0;
				u_pri(i,j,2) = 1.206;
				u_pri(i,j,3) = 0.3;
			}
			
		}
	}
	double x_center = 35.0;
	double y_center = 0.0;
	double r = 25;
	for(int i = boundary_cells;  i < nxCells + boundary_cells; i++){
		for(int j = boundary_cells; j < nyCells + boundary_cells; j++){
			double x = x_0 + (i - boundary_cells + 0.5) * dx;
			double y = y_0 + (j - boundary_cells + 0.5) * dy;
			if(x < 5){
				u_pri(i,j,0) = 1.7755;
				u_pri(i,j,1) = 110.63;
				u_pri(i,j,2) = 0.0;
				u_pri(i,j,3) = 159060.0;
			}else if(sqrt((x - x_center) * (x - x_center) + (y - y_center) * (y - y_center)) > r){
				u_pri(i,j,0) = 1.29;
				u_pri(i,j,1) = 0.0;
				u_pri(i,j,2) = 0.0;
				u_pri(i,j,3) = 101325.0;
			}else{
				u_pri(i,j,0) = 0.214;
				u_pri(i,j,1) = 0.0;
				u_pri(i,j,2) = 0.0;
				u_pri(i,j,3) = 101325.0;
			} 
		}
	}
}


int main(int argc, char *argv[]){
	
	Grid u_con_host(CPU);
	Grid u_pri_host(CPU);
	Grid u_con(GPU);
	Grid u_pri(GPU);
	Grid flux_x_host(nxCells - 1,nyCells,CPU);
	Grid flux_y_host(nxCells,nyCells - 1,CPU);
	Grid flux_x(nxCells - 1,nyCells,GPU);
	Grid flux_y(nxCells,nyCells - 1,GPU);
	Grid GPU_halfL(GPU), GPU_halfR(GPU);
	int blocks = (int)ceil((nxCells + 2 * boundary_cells)/ 32.0) * (int)ceil((nyCells + 2 * boundary_cells)/ 32.0);
	int xblock, yblock;
	
	dim3 dimBlock(32,32,1);
  	dim3 dimGrid((int)ceil((nxCells + 2 * boundary_cells)/ 32.0), (int)ceil((nyCells + 2 * boundary_cells)/ 32.0), 1);
	int xblockx = 6;
	int yblocky = 6;
	
	xblock = (int)ceil(1 + (nxCells + 2 * boundary_cells - xblockx)/(xblockx - 2 * 1) + 1);
	yblock = (int)ceil(1 + (nyCells + 2 * boundary_cells - yblocky)/(yblocky - 2 * 0) + 1);
	dim3 dimBlockx(xblockx,yblocky,1);
	dim3 dimGridx(xblock, yblock, 1);

	xblock = (int)ceil(1 + (nxCells + 2 * boundary_cells - xblockx)/(xblockx - 2 * 0) + 1);
	yblock = (int)ceil(1 + (nyCells + 2 * boundary_cells - yblocky)/(yblocky - 2 * 1) + 1);
	dim3 dimBlocky(xblockx,yblocky,1);
	dim3 dimGridy(xblock, yblock, 1);
	double elapsdt = 0;
	double elapsx = 0;
	double elapsy = 0;
	clock_t startx;
	clock_t endx;
	clock_t starty;
	clock_t endy;
	clock_t startdt;
	clock_t enddt;
	clock_t startini;
	clock_t endini;
	clock_t start;
	clock_t end;
	start = clock();
	GPU_initial_u_pri<<<dimGrid, dimBlock>>>(u_pri,x_0,y_0,dx,dy);
	do{
		GPU_tran_bound_x<<<dimGrid, dimBlock>>>(u_pri);
		CUDA_CHECK;

		GPU_tran_bound_y<<<dimGrid, dimBlock>>>(u_pri);
		CUDA_CHECK;

		GPU_upri2ucon<<<dimGrid, dimBlock>>>(u_con, u_pri);
		CUDA_CHECK;

		dt = caldt(u_pri, blocks, dimGrid, dimBlock);
		CUDA_CHECK;
		
		t = t + dt;
		Overlap_GPU_SLIC<0><<<dimGridx, dimBlockx>>>(u_con, GPU_halfL, GPU_halfR, dx, dy, dt);
		CUDA_CHECK;

		GPU_tran_bound_x<<<dimGrid, dimBlock>>>(u_con);
		CUDA_CHECK;

		GPU_tran_bound_y<<<dimGrid, dimBlock>>>(u_con);
		CUDA_CHECK;
		
		Overlap_GPU_SLIC<1><<<dimGridy, dimBlocky>>>(u_con, GPU_halfL, GPU_halfR, dx, dy, dt);
		CUDA_CHECK;

		GPU_ucon2upri<<<dimGrid, dimBlock>>>(u_pri,u_con);
		CUDA_CHECK;

	}while(t < stop_time);
	cudaFree(u_con.data);
	cudaFree(u_pri.data);
	cudaFree(flux_x.data);
	cudaFree(flux_y.data);
	cudaFree(GPU_halfL.data);
	cudaFree(GPU_halfR.data);
	delete[] u_con_host.data;
	delete[] u_pri_host.data;
	delete[] flux_x_host.data;
	delete[] flux_y_host.data;
	return 0;
}

