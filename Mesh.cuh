#ifndef __MESH_CUH
#define __MESH_CUH
#include <QActFlowDef.cuh>
#include <cudaErr.h>

inline void waveinit(Qreal* kx, Qreal* ky, Qreal* k_squared, int Nxh, int Ny, Qreal Lx, Qreal Ly){
    for (int i = 0; i < Nxh; i++){
        for (int j = 0; j < Ny; j++){
            int index = i + j*Nxh;
            
            if (j<Ny/2+1){
                ky[index] = 2*M_PI/Ly * j;
            }
            else{
                ky[index] = 2*M_PI/Ly * (j-Ny);
            }
            kx[index] = 2*M_PI/Lx * i;
            k_squared[index] = kx[index]*kx[index] + ky[index]*ky[index]; 
        }
    }
}

typedef struct Mesh{
    int Nx, Nxh, Ny;
    Qreal Lx, Ly, dx, dy;
    int specsize, physsize;

    int BSZ;
    dim3 dimGridp, dimBlockp, dimGridsp, dimBlocksp;
    
    cufftHandle transf, inv_transf;

    Qcomp* spec; Qreal* phys;
    Qreal *kx, *ky, *k_squared;
    Mesh(int pNx, int pNy, Qreal pLx, Qreal pLy, int pBSZ):
    Nx(pNx),Ny(pNy),Nxh(Nx/2+1),BSZ(pBSZ), Lx(pLx), Ly(pLy)
    {
        
        specsize = Nxh*Ny*sizeof(Qcomp);
        physsize = Nx*Ny*sizeof(Qreal);
        cudaMalloc((void**)&spec, specsize);
        cudaMalloc((void**)&phys, physsize);

        
        cufft_error_func( cufftPlan2d( &(transf), Ny, Nx, CUFFT_D2Z ) );
        cufft_error_func( cufftPlan2d( &(inv_transf), Ny, Nx, CUFFT_Z2D ) );

        // thread information for physical space
        dimGridp = dim3(int((Nx-0.5)/BSZ) + 1, int((Ny-0.5)/BSZ) + 1);
        dimBlockp = dim3(BSZ, BSZ);

        // thread information for spectral space
        dimGridsp = dim3(int((Nxh-0.5)/BSZ) + 1, int((Ny-0.5)/BSZ) + 1);
        dimBlocksp = dim3(BSZ, BSZ);


        cudaMalloc((void**)&kx, sizeof(Qreal)*Nxh*Ny);
        cudaMalloc((void**)&ky, sizeof(Qreal)*Nxh*Ny);
        cudaMalloc((void**)&k_squared, sizeof(Qreal)*Nxh*Ny);

        Qreal *kxh = (Qreal*)malloc(sizeof(Qreal)*Nxh*Ny);
        Qreal *kyh = (Qreal*)malloc(sizeof(Qreal)*Nxh*Ny);
        Qreal *ksh = (Qreal*)malloc(sizeof(Qreal)*Nxh*Ny);
        waveinit(kxh, kyh, ksh, Nxh, Ny, Lx, Ly);

        cudaMemcpy(kx, kxh, sizeof(Qreal)*Nxh*Ny, cudaMemcpyHostToDevice);
        cudaMemcpy(ky, kyh, sizeof(Qreal)*Nxh*Ny, cudaMemcpyHostToDevice);
        cudaMemcpy(k_squared, ksh, sizeof(Qreal)*Nxh*Ny, cudaMemcpyHostToDevice);

        free(kxh);
        free(kyh);
        free(ksh);

    }

    ~Mesh(){
        cudaFree(spec);
        cudaFree(phys);
    }
} Mesh;

#endif