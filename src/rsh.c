//
//  rsh.c
//  libmsym
//
//  Created by Marcus Johansson on 21/10/15.
//  Copyright (c) 2015 Marcus Johansson.
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

// Ref:
//
// J. Ivanic and K. Ruedenberg, "Rotation Matrices for Real Spherical Harmonics. Direct Determination by Recursion",
// J. Phys. Chem., vol. 100, no. 15, pp. 6342-6347, 1996. http://pubs.acs.org/doi/pdf/10.1021/jp953350u

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include "rsh.h"
#include "linalg.h"
#include "symop.h"

#define SQR(x) ((x)*(x))

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

void rshSymmetryOperationRepresentation(msym_symmetry_operation_t *sops, int index, int l, rsh_representations_t *lts);
void rshCalculateUVWCoefficients(int l, int m1, int m2, double* u, double* v, double* w);
void rshRotationRepresentation(int index, int l, rsh_representations_t *lrs);
double rshRotationP(int index, int l, int i, int m1, int m2, rsh_representations_t *lrs);
double rshRotationU(int index, int l, int m1, int m2, rsh_representations_t *lrs);
double rshRotationV(int index, int l, int m1, int m2, rsh_representations_t *lrs);
double rshRotationW(int index, int l, int m1, int m2, rsh_representations_t *lrs);


msym_error_t generateRSHRepresentations(int sopsl, msym_symmetry_operation_t sops[sopsl], int lmax, rsh_representations_t *lrs){
    msym_error_t ret = MSYM_SUCCESS;
    for(int l = 0;l <= lmax;l++){
        if(lrs[l].d != 2*l+1){
            ret = MSYM_INVALID_BASIS_FUNCTIONS;
            msymSetErrorDetails("Invalid dimension of real spherical harmonic (expected %d, got %d)",2*l+1, lrs[l].d);
            goto err;
        }
        for(int i = 0; i < sopsl;i++){
            rshSymmetryOperationRepresentation(sops,i,l,lrs);
        }
    }
    
err:
    return ret;
}



void rshSymmetryOperationRepresentation(msym_symmetry_operation_t *sops, int index, int l, rsh_representations_t *lrs){
    if(0 == l){
        int d = lrs[0].d;
        double (*st)[d][d] = lrs[0].t;
        st[index][0][0] = 1;
    } else if(1 == l){
        int d = lrs[1].d;
        double (*st)[d][d] = lrs[1].t;
        double (*r)[d] = st[index];
        double t[3][3];
        symmetryOperationMatrix(&sops[index], t);
        // x,y,z -> py, pz, px
        r[0][0] = t[1][1];
        r[0][1] = t[1][2];
        r[0][2] = t[1][0];
        r[1][0] = t[2][1];
        r[1][1] = t[2][2];
        r[1][2] = t[2][0];
        r[2][0] = t[0][1];
        r[2][1] = t[0][2];
        r[2][2] = t[0][0];
    } else {
        int d = lrs[l].d;
        double (*st)[d][d] = lrs[l].t;
        double (*r)[d] = st[index];        
        switch (sops[index].type) {
            case INVERSION:
                if(l & 1){
                    memset(r,0,sizeof(double[d][d]));
                    for(int i = 0;i < d;i++) r[i][i] = -1;
                    break;
                } //fallthrough
            case IDENTITY:
                mleye(d,r);
                break;
            case REFLECTION:
            case IMPROPER_ROTATION:
            case PROPER_ROTATION:
            default:
                rshRotationRepresentation(index, l, lrs);
                break;
                
        }
    }
}

double rshRotationP(int index, int l, int i, int m1, int m2, rsh_representations_t *lrs){
    int d1 = lrs[1].d, dl = lrs[l-1].d, ol = (dl - 1)/2;
    double (*st1)[d1][d1] = lrs[1].t, (*stl)[dl][dl] = lrs[l-1].t;
    double (*r1)[d1] = st1[index], (*rl)[dl] = stl[index];
    double ret = 0;
    
    if (m2 == l) {
        ret = r1[i + 1][2]*rl[m1 + ol][l - 1 + ol] - r1[i + 1][0]*rl[m1 + ol][1 - l + ol];
    } else if (m2 == -l) {
        ret = r1[i + 1][2]*rl[m1 + ol][1 - l + ol] + r1[i + 1][0]*rl[m1 + ol][l - 1 + ol];
    } else {
        ret = r1[i + 1][1]*rl[m1 + ol][m2 + ol];
    }
    
    return ret;
}

double rshRotationU(int index, int l, int m1, int m2, rsh_representations_t *lrs){
    return rshRotationP(index, l, 0, m1, m2, lrs);
}

double rshRotationV(int index, int l, int m1, int m2, rsh_representations_t *lrs){
    double ret = 0;
    if (m1 == 0) {
        ret = rshRotationP(index, l, 1, 1, m2, lrs) + rshRotationP(index, l, -1, -1, m2, lrs);
    } else if (m1 == 1){
        ret = M_SQRT2*rshRotationP(index, l, 1, 0, m2, lrs);
    } else if (m1 == -1){
        ret = M_SQRT2*rshRotationP(index, l, -1, 0, m2, lrs);
    } else if (m1 > 0){
        ret = rshRotationP(index, l, 1, m1 - 1, m2, lrs) - rshRotationP(index,l, -1, -m1 + 1, m2, lrs);
    } else {
        ret = rshRotationP(index, l, 1, m1 + 1, m2, lrs) + rshRotationP(index,l, -1, -m1 - 1, m2, lrs);
    }
    
    return ret;
}

double rshRotationW(int index, int l, int m1, int m2, rsh_representations_t *lrs){
    double ret = 0;
    if (m1 > 0) {
        ret = rshRotationP(index, l, 1, m1 + 1, m2, lrs) + rshRotationP(index, l, -1, -m1 - 1, m2, lrs);
    } else {
        ret = rshRotationP(index, l, 1, m1 - 1, m2, lrs) - rshRotationP(index, l, -1, -m1 + 1, m2, lrs);
    }
    
    return ret;
}

void rshRotationRepresentation(int index, int l, rsh_representations_t *lrs){
    int d = lrs[l].d;
    double (*st)[d][d] = lrs[l].t;
    double (*r)[d] = st[index];
    for (int m1 = -l; m1 <= l; m1++) {
        for (int m2 = -l; m2 <= l; m2++) {
            double u, v, w;
            rshCalculateUVWCoefficients(l, m1, m2, &u, &v, &w);
            if(u != 0){
                u *= rshRotationU(index, l, m1, m2, lrs);
            }
            if(v != 0){
                v *= rshRotationV(index, l, m1, m2, lrs);
            }
            if(w != 0){
                w *= rshRotationW(index, l, m1, m2, lrs);
            }
            r[m1 + l][m2 + l] = u + v + w;
        }
    }
}

void rshCalculateUVWCoefficients(int l, int m1, int m2, double *u, double *v, double *w) {
    if(m1 == 0){ // d = 1
        if(abs(m2) == l){
            *u =  sqrt(l/(4*l - 2.0));
            *v = -0.5*sqrt((l - 1.0)/(2*l - 1.0));
        } else {
            double m22 = SQR(m2), l2 = SQR(l), l2m22 = l2 - m22;
            *u =  sqrt(l2/l2m22);
            *v = -0.5*sqrt(2*(l2 - l)/l2m22);
        }
        *w = 0;
    } else { // d = 0
        double am1 = abs(m1);
        double div = (abs(m2) == l ? 2.0*l*(2.0*l - 1) : (l + m2)*(l - m2));
        *u =  sqrt((l + m1)*(l - m1)/div);
        *v =  0.5*sqrt((l + am1 - 1)*(l + am1)/div);
        *w = -0.5*sqrt((l - am1 - 1)*(l - am1)/div);
    }
}
