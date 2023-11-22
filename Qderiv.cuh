#ifndef QDERIV_CUH
#define QDERIV_CUH

#include <Mesh.cuh>
#include <cuComplexBinOp.cuh>

__global__ 
void coeff(Qreal *f, int Nx, int Ny, int BSZ){
    int i = blockIdx.x * BSZ + threadIdx.x;
    int j = blockIdx.y * BSZ + threadIdx.y;
    int index = j*Nx + i;
    if(i<Nx && j<Ny){
        f[index] = f[index]/(Nx*Ny);
    }
}

void FwdTrans(Qcomp* ft, Qreal* f, Mesh* mesh){
    cudaMemcpy(mesh->phys, f, mesh->physsize, cudaMemcpyDeviceToDevice);
    cufft_error_func( cufftExecD2Z(mesh->transf, mesh->phys, ft));
}

void BwdTrans( Qreal* f, Qcomp* ft, Mesh* mesh){
    cudaMemcpy(mesh->spec, ft, mesh->specsize, cudaMemcpyDeviceToDevice);
    cufft_error_func( cufftExecZ2D(mesh->inv_transf, mesh->spec, f));
    coeff<<<mesh->dimGridp,mesh->dimBlockp>>>(f, mesh->Nx, mesh->Ny, mesh->BSZ);
}

__global__
void vel_funcD(Qcomp *w_c, Qcomp *u_c, Qcomp *v_c, 
Qreal* k_squared, Qreal* kx, Qreal*ky, int Nxh, int Ny, int BSZ){
    int i = blockIdx.x * BSZ + threadIdx.x;
    int j = blockIdx.y * BSZ + threadIdx.y;
    int index = j*Nxh + i;
    if (i<Nxh && j<Ny){
        if (i==0 && j==0)
        {
            u_c[index] = make_cuDoubleComplex(0.0,0.0);
            v_c[index] = make_cuDoubleComplex(0.0,0.0);
        }
        else{
            //u = -D_y(\phi) -> u_spec = -1 * i* ky* w_spec/(-1* (kx^2+ky^2) )
            u_c[index] = ky[index]*im()*w_c[index]/(k_squared[index]);
            //v = D_x(\phi) -> v_spec = i* kx* w_spec/(-1* (kx^2+ky^2) )
            v_c[index] = -1.0*kx[index]*im()*w_c[index]/(k_squared[index]);
        }
    }
}
inline void vel_func(Qcomp *w_c, Qcomp *u_c, Qcomp *v_c, Mesh *mesh){
    vel_funcD<<<mesh->dimGridsp, mesh->dimBlocksp>>>(w_c, u_c, v_c, mesh->k_squared, mesh->kx, mesh->ky, mesh->Nxh, mesh->Ny, mesh->BSZ);
}


#endif