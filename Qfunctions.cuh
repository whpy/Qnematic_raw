#ifndef QFUNCTIONS_CUH
#define QFUNCTIONS_CUH

#include <iostream>
#include <fstream>
#include <math.h>
#include <QActFlowDef.cuh>

using std::string;
using std::endl;
using std::ofstream;
using std::cout;

inline void winit(Qreal* w, int Nx, int Ny, Qreal dx, Qreal dy){
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            int index = i + j*Nx;
            Qreal x = i*dx;
            Qreal y = j*dy;
            w[index] = -1*sin(x+y);
        }
    }
}

inline void r1init(Qreal* w, int Nx, int Ny, Qreal dx, Qreal dy){
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            int index = i + j*Nx;
            Qreal x = i*dx;
            Qreal y = j*dy;
            w[index] = -25*sin(3*x+4*y);
        }
    }
}

inline void r2init(Qreal* w, int Nx, int Ny, Qreal dx, Qreal dy){
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            int index = i + j*Nx;
            Qreal x = i*dx;
            Qreal y = j*dy;
            w[index] = -25*sin(3*x+4*y);
        }
    }
}

inline void print_spec(Qcomp* f, int Nxh, int Ny){
    for(int j = 0; j < Ny; j++){
        for (int i = 0; i < Nxh; i++){
            int index = i + j*Nxh;
            printf("(%.2f, %.2f)  ", f[index].x, f[index].y);
        }
        cout << endl;
    }
    cout << endl;
}

inline void print_spec(Qreal* f, int Nxh, int Ny){
    for(int j = 0; j < Ny; j++){
        for (int i = 0; i < Nxh; i++){
            int index = i + j*Nxh;
            printf("%.2f  ", f[index]);
        }
        cout << endl;
    }
    cout << endl;
}

inline void field_visual(Qreal *f, string name, int Nx, int Ny){
    ofstream fval;
    string fname = name;
    fval.open(fname);
    for (int j=0; j<Ny; j++){
        for (int i=0; i<Nx; i++){
            int index = j*Nx + i;
            fval << f[index] << ",";
        }
        fval << endl;
    }
    fval.close();
}

inline void coord(Qreal dx, Qreal dy, int Nx, int Ny){
    ofstream xcoord("x.csv");
    ofstream ycoord("y.csv");
    for (int j=0; j<Ny; j++){
        for ( int i=0; i< Nx; i++){
            float x = dx*i;
            float y = dy*j;
            xcoord << x << ",";
            ycoord << y << ",";
        }
        xcoord << endl;
        ycoord << endl;
    }
    xcoord.close();
    ycoord.close();
}

#endif