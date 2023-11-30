
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <chrono>
#include <random>

#include <cuComplexBinOp.cuh>
#include <cudaErr.h>
#include <QActFlowDef.cuh>

#include <Mesh.cuh>
#include <Qfunctions.cuh>
#include <QFldfuncs.cuh>
// #include <Qderiv.cuh>

// #define M_PI 3.141592653589

using namespace std;
// struct to faciliate
// typedef struct Field{
//     Qreal *phys;
//     Qcomp *spec;
//     Field(Qreal *pphys, Qcomp* pspec):phys(pphys), spec(pspec){}
// } Field;
// inline void waveinit(Qreal* kx, Qreal* ky, Qreal* k_squared, int Nxh, int Ny, Qreal Lx, Qreal Ly);
// funtions
// __global__ 
// void coeff(Qreal *f, int Nx, int Ny, int BSZ){
//     int i = blockIdx.x * BSZ + threadIdx.x;
//     int j = blockIdx.y * BSZ + threadIdx.y;
//     int index = j*Nx + i;
//     if(i<Nx && j<Ny){
//         f[index] = f[index]/(Nx*Ny);
//     }
// }
// void FwdTrans(Qcomp* ft, Qreal* f, Mesh* mesh){
//     cudaMemcpy(mesh->phys, f, mesh->physsize, cudaMemcpyDeviceToDevice);
//     cufft_error_func( cufftExecD2Z(mesh->transf, mesh->phys, ft));
// }

// void BwdTrans( Qreal* f, Qcomp* ft, Mesh* mesh){
//     cudaMemcpy(mesh->spec, ft, mesh->specsize, cudaMemcpyDeviceToDevice);
//     cufft_error_func( cufftExecZ2D(mesh->inv_transf, mesh->spec, f));
//     coeff<<<mesh->dimGridp,mesh->dimBlockp>>>(f, mesh->Nx, mesh->Ny, mesh->BSZ);
// }

// inline void field_visual(Qreal *f, string name, int Nx, int Ny){
//     ofstream fval;
//     string fname = name;
//     fval.open(fname);
//     for (int j=0; j<Ny; j++){
//         for (int i=0; i<Nx; i++){
//             int index = j*Nx + i;
//             fval << f[index] << ",";
//         }
//         fval << endl;
//     }
//     fval.close();
// }

// inline void waveinit(Qreal* kx, Qreal* ky, Qreal* k_squared, int Nxh, int Ny, Qreal Lx, Qreal Ly){
//     for (int i = 0; i < Nxh; i++){
//         for (int j = 0; j < Ny; j++){
//             int index = i + j*Nxh;
            
//             if (j<Ny/2+1){
//                 ky[index] = 2*M_PI/Ly * j;
//             }
//             else{
//                 ky[index] = 2*M_PI/Ly * (j-Ny);
//             }
//             kx[index] = 2*M_PI/Lx * i;
//             k_squared[index] = kx[index]*kx[index] + ky[index]*ky[index]; 
//         }
//     }
// }

// inline void winit(Qreal* w, int Nx, int Ny, Qreal dx, Qreal dy){
//     for (int j = 0; j < Ny; j++){
//         for (int i = 0; i < Nx; i++){
//             int index = i + j*Nx;
//             Qreal x = i*dx;
//             Qreal y = j*dy;
//             w[index] = -1*sin(x+y);
//         }
//     }
// }

// inline void r1init(Qreal* w, int Nx, int Ny, Qreal dx, Qreal dy){
//     for (int j = 0; j < Ny; j++){
//         for (int i = 0; i < Nx; i++){
//             int index = i + j*Nx;
//             Qreal x = i*dx;
//             Qreal y = j*dy;
//             w[index] = -25*sin(3*x+4*y);
//         }
//     }
// }

// inline void r2init(Qreal* w, int Nx, int Ny, Qreal dx, Qreal dy){
//     for (int j = 0; j < Ny; j++){
//         for (int i = 0; i < Nx; i++){
//             int index = i + j*Nx;
//             Qreal x = i*dx;
//             Qreal y = j*dy;
//             w[index] = -25*sin(3*x+4*y);
//         }
//     }
// }

// inline void print_spec(Qcomp* f, int Nxh, int Ny){
//     for(int j = 0; j < Ny; j++){
//         for (int i = 0; i < Nxh; i++){
//             int index = i + j*Nxh;
//             printf("(%.2f, %.2f)  ", f[index].x, f[index].y);
//         }
//         cout << endl;
//     }
//     cout << endl;
// }

// inline void print_spec(Qreal* f, int Nxh, int Ny){
//     for(int j = 0; j < Ny; j++){
//         for (int i = 0; i < Nxh; i++){
//             int index = i + j*Nxh;
//             printf("%.2f  ", f[index]);
//         }
//         cout << endl;
//     }
//     cout << endl;
// }

// __global__
// void vel_funcD(Qcomp *w_c, Qcomp *u_c, Qcomp *v_c, 
// Qreal* k_squared, Qreal* kx, Qreal*ky, int Nxh, int Ny, int BSZ){
//     int i = blockIdx.x * BSZ + threadIdx.x;
//     int j = blockIdx.y * BSZ + threadIdx.y;
//     int index = j*Nxh + i;
//     if (i<Nxh && j<Ny){
//         if (i==0 && j==0)
//         {
//             u_c[index] = make_cuDoubleComplex(0.0,0.0);
//             v_c[index] = make_cuDoubleComplex(0.0,0.0);
//         }
//         else{
//             //u = -D_y(\phi) -> u_spec = -1 * i* ky* w_spec/(-1* (kx^2+ky^2) )
//             u_c[index] = ky[index]*im()*w_c[index]/(k_squared[index]);
//             //v = D_x(\phi) -> v_spec = i* kx* w_spec/(-1* (kx^2+ky^2) )
//             v_c[index] = -1.0*kx[index]*im()*w_c[index]/(k_squared[index]);
//         }
//     }
// }
// inline void vel_func(Qcomp *w_c, Qcomp *u_c, Qcomp *v_c, Mesh *mesh){
//     vel_funcD<<<mesh->dimGridsp, mesh->dimBlocksp>>>(w_c, u_c, v_c, mesh->k_squared, mesh->kx, mesh->ky, mesh->Nxh, mesh->Ny, mesh->BSZ);
// }

// inline void coord(Qreal dx, Qreal dy, int Nx, int Ny){
//     ofstream xcoord("x.csv");
//     ofstream ycoord("y.csv");
//     for (int j=0; j<Ny; j++){
//         for ( int i=0; i< Nx; i++){
//             float x = dx*i;
//             float y = dy*j;
//             xcoord << x << ",";
//             ycoord << y << ",";
//         }
//         xcoord << endl;
//         ycoord << endl;
//     }
//     xcoord.close();
//     ycoord.close();
// }

int main(){
    int Nx = 8;
    int Ny = Nx;
    int BSZ = 16;
    int Nxh = Nx/2+1;
    int specsize = Nxh*Ny*sizeof(Qcomp);
    int physize = Nx*Ny*sizeof(Qreal);
    int wavesize = Nxh*Ny*sizeof(Qreal);
    Qreal Lx = 2*M_PI;
    Qreal Ly = Lx;
    Qreal dx = Lx/Nx;
    Qreal dy = Ly/Ny;

    cufftHandle transf;
    cufftHandle inv_transf;
    cufft_error_func( cufftPlan2d( &(transf), Ny, Nx, CUFFT_D2Z ) );
    cufft_error_func( cufftPlan2d( &(inv_transf), Ny, Nx, CUFFT_Z2D ) );

    dim3 dimGridp = dim3(int((Nx-0.5)/BSZ) + 1, int((Ny-0.5)/BSZ) + 1);
    dim3 dimBlockp = dim3(BSZ, BSZ);

    dim3 dimGridsp = dim3(int((Nxh-0.5)/BSZ) + 1, int((Ny-0.5)/BSZ) + 1);
    dim3 dimBlocksp = dim3(BSZ, BSZ);

    Mesh *mesh = new Mesh(Nx, Ny, Lx, Ly, BSZ);
    coord(dx, dy, Nx, Ny);

    Qreal *w;
    Qcomp *w_c, *u_c, *v_c;
    Qreal *u, *v;
    Qreal *kx, *ky, *k_squared;

    Qreal *p1, *p2, *p3;
    Qcomp *sp1, *sp2, *sp3;
    cudaMalloc((void**)&w, sizeof(Qreal)*Nx*Ny);
    cudaMalloc((void**)&u, sizeof(Qreal)*Nx*Ny);
    cudaMalloc((void**)&v, sizeof(Qreal)*Nx*Ny);

    cudaMalloc((void**)&w_c, sizeof(Qcomp)*Nxh*Ny);
    cudaMalloc((void**)&u_c, sizeof(Qcomp)*Nxh*Ny);
    cudaMalloc((void**)&v_c, sizeof(Qcomp)*Nxh*Ny);

    cudaMalloc((void**)&kx, sizeof(Qreal)*Nxh*Ny);
    cudaMalloc((void**)&ky, sizeof(Qreal)*Nxh*Ny);
    cudaMalloc((void**)&k_squared, sizeof(Qreal)*Nxh*Ny);

    p1 = (Qreal*)malloc(sizeof(Qreal)*Nx*Ny);
    p2 = (Qreal*)malloc(sizeof(Qreal)*Nx*Ny);
    p3 = (Qreal*)malloc(sizeof(Qreal)*Nx*Ny);
    sp1 = (Qcomp*)malloc(sizeof(Qcomp)*Nxh*Ny); 
    sp2 = (Qcomp*)malloc(sizeof(Qcomp)*Nxh*Ny);
    sp3 = (Qcomp*)malloc(sizeof(Qcomp)*Nxh*Ny);

    Qreal *kxh = (Qreal*)malloc(sizeof(Qreal)*Nxh*Ny);
    Qreal *kyh = (Qreal*)malloc(sizeof(Qreal)*Nxh*Ny);
    Qreal *ksh = (Qreal*)malloc(sizeof(Qreal)*Nxh*Ny);
    
    
    winit(p1, Nx, Ny, dx, dy);
    cudaMemcpy(w, p1, sizeof(Qreal)*Nx*Ny, cudaMemcpyHostToDevice);
    FwdTrans(w_c, w, mesh);
    // cufft_error_func( cufftExecD2Z(transf, w, w_c));
    cudaMemcpy(p2, w, sizeof(Qreal)*Nx*Ny, cudaMemcpyDeviceToHost);
    field_visual(p2, "w0.csv", Nx, Ny);

    cudaMemcpy(sp2, w_c, sizeof(Qcomp)*Nxh*Ny, cudaMemcpyDeviceToHost);
    cuda_error_func(cudaDeviceSynchronize());
    cout << "w_c" << endl;
    print_spec(sp2, Nxh, Ny);
    cuda_error_func(cudaDeviceSynchronize());
    
    // cuda_error_func(cudaDeviceSynchronize());

    // BwdTrans(w, w_c, Nx, Ny, BSZ, dimGridp, dimBlockp, inv_transf);
    BwdTrans(w, w_c, mesh);
    // cufft_error_func( cufftExecZ2D(inv_transf, w_c, w));
    // cudaMemcpy(sp2, w_c, sizeof(Qcomp)*Nxh*Ny, cudaMemcpyDeviceToHost);
    // cuda_error_func(cudaDeviceSynchronize());
    // cout << "w_c 1" << endl;
    // print_spec(sp2, Nxh, Ny);
    // cuda_error_func(cudaDeviceSynchronize());
    // coeff<<<dimGridp, dimBlockp>>>(w, Nx, Ny, BSZ);
    // cudaMemcpy(sp2, w_c, sizeof(Qcomp)*Nxh*Ny, cudaMemcpyDeviceToHost);
    // cuda_error_func(cudaDeviceSynchronize());
    // cout << "w_c 2" << endl;
    // print_spec(sp2, Nxh, Ny);
    cudaMemcpy(p2, w, sizeof(Qreal)*Nx*Ny, cudaMemcpyDeviceToHost);
    field_visual(p2, "w.csv", Nx, Ny);

    cufft_error_func( cufftExecD2Z(transf, w, w_c));
    cudaMemcpy(sp2, w_c, sizeof(Qcomp)*Nxh*Ny, cudaMemcpyDeviceToHost);
    cuda_error_func(cudaDeviceSynchronize());
    cout << "w_c 3" << endl;
    print_spec(sp2, Nxh, Ny);
    cuda_error_func(cudaDeviceSynchronize());

    waveinit(kxh, kyh, ksh, Nxh, Ny, Lx, Ly);
    cudaMemcpy(kx, kxh, sizeof(Qreal)*Nxh*Ny, cudaMemcpyHostToDevice);
    cudaMemcpy(ky, kyh, sizeof(Qreal)*Nxh*Ny, cudaMemcpyHostToDevice);
    cudaMemcpy(k_squared, ksh, sizeof(Qreal)*Nxh*Ny, cudaMemcpyHostToDevice);

    cudaMemcpy(kxh, mesh->kx, sizeof(Qreal)*Nxh*Ny, cudaMemcpyDeviceToHost);
    cout << "kx" << endl;
    print_spec(kxh,Nxh,Ny);
    cudaMemcpy(kyh, mesh->ky, sizeof(Qreal)*Nxh*Ny, cudaMemcpyDeviceToHost);
    cout << "ky" << endl;
    print_spec(kyh,Nxh,Ny);
    cudaMemcpy(ksh, mesh->k_squared, sizeof(Qreal)*Nxh*Ny, cudaMemcpyDeviceToHost);
    cout << "k_squared" << endl;
    print_spec(ksh,Nxh,Ny);

    cuda_error_func(cudaDeviceSynchronize());
    // vel_func<<<dimGridsp, dimBlocksp>>>(w_c, u_c, v_c, mesh->k_squared, mesh->kx, mesh->ky, mesh->Nxh, mesh->Ny, mesh->BSZ);
    vel_func(w_c, u_c, v_c, mesh);
    cudaMemcpy(sp2, u_c, sizeof(Qcomp)*Nxh*Ny, cudaMemcpyDeviceToHost);
    cuda_error_func(cudaDeviceSynchronize());
    cout << "u_c" << endl;
    print_spec(sp2, Nxh, Ny);
    cudaMemcpy(sp2, v_c, sizeof(Qcomp)*Nxh*Ny, cudaMemcpyDeviceToHost);
    cuda_error_func(cudaDeviceSynchronize());
    cout << "v_c" << endl;
    print_spec(sp2,Nxh, Ny);
    cuda_error_func(cudaDeviceSynchronize());
    // BwdTrans(u, u_c, Nx, Ny, BSZ, dimGridp, dimBlockp, inv_transf);
    BwdTrans(u,u_c,mesh);
    // cufft_error_func( cufftExecZ2D(inv_transf, u_c, u));
    // coeff<<<dimGridp, dimBlockp>>>(u, Nx, Ny, BSZ);
    // BwdTrans(v, v_c, Nx, Ny, BSZ, dimGridp, dimBlockp, inv_transf);
    BwdTrans(v,v_c,mesh);
    // cufft_error_func( cufftExecZ2D(inv_transf, v_c, v));
    // coeff<<<dimGridp, dimBlockp>>>(v, Nx, Ny, BSZ);
    
    cudaMemcpy(p1, u, sizeof(Qreal)*Nx*Ny, cudaMemcpyDeviceToHost);
    field_visual(p1, "u.csv", Nx, Ny);
    cudaMemcpy(p2, v, sizeof(Qreal)*Nx*Ny, cudaMemcpyDeviceToHost);
    field_visual(p2, "v.csv", Nx, Ny);

    delete mesh;
    return 0;

}
// int main(){
//     int Nx = 8;
//     int Ny = Nx;
//     int BSZ = 16;
//     int Nxh = Nx/2+1;
//     int specsize = Nxh*Ny*sizeof(Qcomp);
//     int physize = Nx*Ny*sizeof(Qreal);
//     int wavesize = Nxh*Ny*sizeof(Qreal);
//     Qreal Lx = 2*M_PI;
//     Qreal Ly = Lx;
//     Qreal dx = Lx/Nx;
//     Qreal dy = Ly/Ny;

//     cufftHandle transf;
//     cufftHandle inv_transf;
//     cufft_error_func( cufftPlan2d( &(transf), Ny, Nx, CUFFT_D2Z ) );
//     cufft_error_func( cufftPlan2d( &(inv_transf), Ny, Nx, CUFFT_Z2D ) );

//     dim3 dimGridp = dim3(int((Nx-0.5)/BSZ) + 1, int((Ny-0.5)/BSZ) + 1);
//     dim3 dimBlockp = dim3(BSZ, BSZ);

//     dim3 dimGridsp = dim3(int((Nxh-0.5)/BSZ) + 1, int((Ny-0.5)/BSZ) + 1);
//     dim3 dimBlocksp = dim3(BSZ, BSZ);


//     // define parameters
//     Qreal Rf = 0.0000075;
//     Qreal lambda = 0.1;
//     Qreal Re = 0.1;
//     Qreal Er = 0.1;
//     coord(dx, dy, Nx, Ny);

//     Qreal *w;
//     Qcomp *w_c, *u_c, *v_c;
//     Qreal *u, *v;
//     Qreal *kx, *ky, *k_squared;

//     Qreal *p1, *p2, *p3;
//     Qcomp *sp1, *sp2, *sp3;
//     cudaMalloc((void**)&w, sizeof(Qreal)*Nx*Ny);
//     cudaMalloc((void**)&u, sizeof(Qreal)*Nx*Ny);
//     cudaMalloc((void**)&v, sizeof(Qreal)*Nx*Ny);

//     cudaMalloc((void**)&w_c, sizeof(Qcomp)*Nxh*Ny);
//     cudaMalloc((void**)&u_c, sizeof(Qcomp)*Nxh*Ny);
//     cudaMalloc((void**)&v_c, sizeof(Qcomp)*Nxh*Ny);

//     cudaMalloc((void**)&kx, sizeof(Qreal)*Nxh*Ny);
//     cudaMalloc((void**)&ky, sizeof(Qreal)*Nxh*Ny);
//     cudaMalloc((void**)&k_squared, sizeof(Qreal)*Nxh*Ny);

//     p1 = (Qreal*)malloc(sizeof(Qreal)*Nx*Ny);
//     p2 = (Qreal*)malloc(sizeof(Qreal)*Nx*Ny);
//     p3 = (Qreal*)malloc(sizeof(Qreal)*Nx*Ny);
//     sp1 = (Qcomp*)malloc(sizeof(Qcomp)*Nxh*Ny); 
//     sp2 = (Qcomp*)malloc(sizeof(Qcomp)*Nxh*Ny);
//     sp3 = (Qcomp*)malloc(sizeof(Qcomp)*Nxh*Ny);

//     Qreal *kxh = (Qreal*)malloc(sizeof(Qreal)*Nxh*Ny);
//     Qreal *kyh = (Qreal*)malloc(sizeof(Qreal)*Nxh*Ny);
//     Qreal *ksh = (Qreal*)malloc(sizeof(Qreal)*Nxh*Ny);
    
    
//     winit(p1, Nx, Ny, dx, dy);
//     cudaMemcpy(w, p1, sizeof(Qreal)*Nx*Ny, cudaMemcpyHostToDevice);
//     // FwdTrans(w_c, w, transf);
//     cufft_error_func( cufftExecD2Z(transf, w, w_c));
//     cudaMemcpy(p2, w, sizeof(Qreal)*Nx*Ny, cudaMemcpyDeviceToHost);
//     field_visual(p2, "w0.csv", Nx, Ny);
//     cudaMemcpy(sp2, w_c, sizeof(Qcomp)*Nxh*Ny, cudaMemcpyDeviceToHost);
//     cuda_error_func(cudaDeviceSynchronize());
//     cout << "w_c" << endl;
//     print_spec(sp2, Nxh, Ny);
//     cuda_error_func(cudaDeviceSynchronize());
    
//     // cuda_error_func(cudaDeviceSynchronize());

//     // BwdTrans(w, w_c, Nx, Ny, BSZ, dimGridp, dimBlockp, inv_transf);
//     cufft_error_func( cufftExecZ2D(inv_transf, w_c, w));
//     // cudaMemcpy(sp2, w_c, sizeof(Qcomp)*Nxh*Ny, cudaMemcpyDeviceToHost);
//     // cuda_error_func(cudaDeviceSynchronize());
//     // cout << "w_c 1" << endl;
//     // print_spec(sp2, Nxh, Ny);
//     // cuda_error_func(cudaDeviceSynchronize());
//     // coeff<<<dimGridp, dimBlockp>>>(w, Nx, Ny, BSZ);
//     // cudaMemcpy(sp2, w_c, sizeof(Qcomp)*Nxh*Ny, cudaMemcpyDeviceToHost);
//     // cuda_error_func(cudaDeviceSynchronize());
//     // cout << "w_c 2" << endl;
//     // print_spec(sp2, Nxh, Ny);
//     cudaMemcpy(p2, w, sizeof(Qreal)*Nx*Ny, cudaMemcpyDeviceToHost);
//     field_visual(p2, "w.csv", Nx, Ny);

//     cufft_error_func( cufftExecD2Z(transf, w, w_c));
//     cudaMemcpy(sp2, w_c, sizeof(Qcomp)*Nxh*Ny, cudaMemcpyDeviceToHost);
//     cuda_error_func(cudaDeviceSynchronize());
//     cout << "w_c 3" << endl;
//     print_spec(sp2, Nxh, Ny);
//     cuda_error_func(cudaDeviceSynchronize());

//     waveinit(kxh, kyh, ksh, Nxh, Ny, Lx, Ly);
//     cudaMemcpy(kx, kxh, sizeof(Qreal)*Nxh*Ny, cudaMemcpyHostToDevice);
//     cudaMemcpy(ky, kyh, sizeof(Qreal)*Nxh*Ny, cudaMemcpyHostToDevice);
//     cudaMemcpy(k_squared, ksh, sizeof(Qreal)*Nxh*Ny, cudaMemcpyHostToDevice);

//     cudaMemcpy(kxh, kx, sizeof(Qreal)*Nxh*Ny, cudaMemcpyDeviceToHost);
//     cout << "kxh" << endl;
//     print_spec(kxh,Nxh,Ny);
//     cudaMemcpy(kyh, ky, sizeof(Qreal)*Nxh*Ny, cudaMemcpyDeviceToHost);
//     cout << "kyh" << endl;
//     print_spec(kyh,Nxh,Ny);
//     cudaMemcpy(ksh, k_squared, sizeof(Qreal)*Nxh*Ny, cudaMemcpyDeviceToHost);
//     cout << "ksh" << endl;
//     print_spec(ksh,Nxh,Ny);

//     cuda_error_func(cudaDeviceSynchronize());
//     vel_func<<<dimGridsp, dimBlocksp>>>(w_c, u_c, v_c, k_squared, kx, ky, Nxh, Ny, BSZ);
//     cudaMemcpy(sp2, u_c, sizeof(Qcomp)*Nxh*Ny, cudaMemcpyDeviceToHost);
//     cuda_error_func(cudaDeviceSynchronize());
//     cout << "u_c" << endl;
//     print_spec(sp2, Nxh, Ny);
//     cudaMemcpy(sp2, v_c, sizeof(Qcomp)*Nxh*Ny, cudaMemcpyDeviceToHost);
//     cuda_error_func(cudaDeviceSynchronize());
//     cout << "v_c" << endl;
//     print_spec(sp2,Nxh, Ny);
//     cuda_error_func(cudaDeviceSynchronize());
//     // BwdTrans(u, u_c, Nx, Ny, BSZ, dimGridp, dimBlockp, inv_transf);
//     cufft_error_func( cufftExecZ2D(inv_transf, u_c, u));
//     coeff<<<dimGridp, dimBlockp>>>(u, Nx, Ny, BSZ);
//     // BwdTrans(v, v_c, Nx, Ny, BSZ, dimGridp, dimBlockp, inv_transf);
//     cufft_error_func( cufftExecZ2D(inv_transf, v_c, v));
//     coeff<<<dimGridp, dimBlockp>>>(v, Nx, Ny, BSZ);
    
//     cudaMemcpy(p1, u, sizeof(Qreal)*Nx*Ny, cudaMemcpyDeviceToHost);
//     field_visual(p1, "u.csv", Nx, Ny);
//     cudaMemcpy(p2, v, sizeof(Qreal)*Nx*Ny, cudaMemcpyDeviceToHost);
//     field_visual(p2, "v.csv", Nx, Ny);

//     return 0;
// }