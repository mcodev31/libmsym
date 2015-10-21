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
#include "msym.h"
#include "symop.h"

#define SQR(x) ((x)*(x))

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880
#endif

typedef struct _rsh_representations {
    int d;
    void *t;
} rsh_representations_t;

void rshSymmetryOperationRepresentation(msym_point_group_t *pg, int index, int l, rsh_representations_t *lts);
void rshCalculateUVWCoefficients(int m1, int m2, int l, double* u, double* v, double* w);
void rshProperRotationRepresentation(int index, int l, rsh_representations_t *lrs);

void generateTransforms(){
    int lmax = 2;
    msym_point_group_t mpg;
    msym_point_group_t *pg = &mpg;
    pg->order = 1;
    msym_symmetry_operation_t sop = {.order = 3, .power = 1, .v = {0,0,1}, .type = PROPER_ROTATION};
    pg->sops = &sop;
    rsh_representations_t *lrs = calloc(lmax+1,sizeof(*lrs));
    for(int l = 0; l <= lmax;l++){
        int d = 2*l+1;
        lrs[l].t = calloc(pg->order+1,sizeof(double[d][d]));
        lrs[l].d = d;
    }
    
    
    for(int l = 0;l <= lmax;l++){
        for(int i = 0; i < pg->order;i++){
            rshSymmetryOperationRepresentation(pg,i,l,lrs);
        }
    }
    
    
    for(int i = 0; i < pg->order;i++){
        printSymmetryOperation(&pg->sops[i]);
        for(int l = 0;l <= lmax;l++){
            printTransform(lrs[l].d,lrs[l].d,lrs[l].t);
        }
    }
    
}



void rshSymmetryOperationRepresentation(msym_point_group_t *pg, int index, int l, rsh_representations_t *lrs){
    if(0 == l){
        int d = lrs[0].d;
        double (*st)[d][d] = lrs[0].t;
        st[index][0][0] = 1;
    } else if(1 == l){
        int d = lrs[1].d;
        double (*st)[d][d] = lrs[1].t;
        symmetryOperationMatrix(&pg->sops[index], st[index]);
    } else {
        int d = lrs[l].d;
        double (*t)[d][d] = calloc(1, sizeof(*t));
        generateBasisFunctionTransforms(1,pg->sops,l,t);
        printTransform(d,d,t[0]);
        free(t);
        switch (pg->sops[index].type) {
            case PROPER_ROTATION:
                rshProperRotationRepresentation(index, l, lrs);
                break;
                
            default:
                printf("can only handle proper rotation");
                exit(1);
                break;
        }
    }
}

double rshProperRotationP(int index, int l, int i, int m1, int m2, rsh_representations_t *lrs){
    int d1 = lrs[1].d, dl = lrs[l-1].d, o1 = 1, ol = (dl - 1)/2;
    double (*st1)[d1][d1] = lrs[1].t, (*stl)[dl][dl] = lrs[l-1].t;
    double (*r1)[d1] = st1[index], (*rl)[dl] = stl[index];
    double ret = 0;
    
    if (m2 == l) {
        ret = r1[i + o1][1 + o1]*rl[m1 + ol][l - 1 + ol] - r1[i + o1][-1 + o1]*rl[m1 + ol][1 - l + ol];
    } else if (m2 == -l) {
        ret = r1[i + o1][1 + o1]*rl[m1 + ol][1 - l + ol] - r1[i + o1][-1 + o1]*rl[m1 + ol][l - 1 + ol];
    } else {
        ret = r1[i + o1][o1]*rl[m1 + ol][m2 + ol];
    }
    
    return ret;
}

double rshProperRotationU(int index, int l, int m1, int m2, rsh_representations_t *lrs){
    return rshProperRotationP(index, l, 0, m1, m2, lrs);
}

double rshProperRotationV(int index, int l, int m1, int m2, rsh_representations_t *lrs){
    double ret = 0;
    if (m1 == 0) {
        ret = rshProperRotationP(index, l, 1, 1, m2, lrs) + rshProperRotationP(index, l, -1, -1, m2, lrs);
    } else if (m1 == 1){
        ret = M_SQRT2*rshProperRotationP(index, l, 1, 0, m2, lrs);
    } else if (m1 == -1){
        ret = M_SQRT2*rshProperRotationP(index, l, -1, 0, m2, lrs);
    } else if (m1 > 0){
        ret = rshProperRotationP(index, l, 1, m1 - 1, m2, lrs) + rshProperRotationP(index,l, -1,-m1 + 1, m2, lrs);
    } else {
        ret = rshProperRotationP(index, l, 1, m1 + 1, m2, lrs) + rshProperRotationP(index,l, -1,-m1 - 1, m2, lrs);
    }
    
    return ret;
}

double rshProperRotationW(int index, int l, int m1, int m2, rsh_representations_t *lrs){
    double ret = 0;
    if (m1 > 0) {
        ret = rshProperRotationP(index, l, 1, m1 + 1, m2, lrs) + rshProperRotationP(index, l, -1, -m1 - 1, m2, lrs);
    } else {
        ret = rshProperRotationP(index, l, 1, m1 - 1, m2, lrs) + rshProperRotationP(index, l, -1, -m1 + 1, m2, lrs);
    }
    
    return ret;
}

void rshProperRotationRepresentation(int index, int l, rsh_representations_t *lrs){
    int d = lrs[l].d;
    double (*st)[d][d] = lrs[l].t;
    double (*r)[d] = st[index];
    for (int m1 = -l; m1 <= l; m1++) {
        for (int m2 = -l; m2 <= l; m2++) {
            double u, v, w;
            rshCalculateUVWCoefficients(m1, m2, l, &u, &v, &w);
            if(u != 0){
                u *= rshProperRotationU(index, l, m1, m2, lrs);
            }
            if(v != 0){
                v *= rshProperRotationV(index, l, m1, m2, lrs);
            }
            if(w != 0){
                w *= rshProperRotationW(index, l, m1, m2, lrs);
            }
            
            r[m1 + l][m2 + l] = u + v + w;
        }
    }
}

#define RSH_EPSILON 16

void rshCalculateUVWCoefficients(int m1, int m2, int l, double* ou, double* ov, double* ow) {
    double u = 0, v = 0, w = 0;
    if(m1 == 0){ // d = 1
        if(abs(m2) == l){
            u =  sqrt(l/(4*l - 2.0));
            v = -0.5*sqrt((l - 1.0)/(2*l - 1.0));
        } else {
            double m22 = SQR(m2), l2 = SQR(l), l2m22 = l2 - m22;
            u =  sqrt(l2/l2m22);
            v = -0.5*sqrt(2*(l2 - l)/l2m22);
        }
        w = 0;
    } else { // d = 0
        double am1 = abs(m1);
        double div = (abs(m2) == l ? 2.0*l*(2.0*l - 1) : (l + m2)*(l - m2));
        u =  sqrt((l + m1)*(l - m1)/div);
        v =  0.5*sqrt((l + am1 - 1)*(l + am1)/div);
        w = -0.5*sqrt((l - am1 - 1)*(l - am1)/div);
    }
    
    *ou = fabs(u) > DBL_EPSILON*RSH_EPSILON ? u : 0.0;
    *ov = fabs(v) > DBL_EPSILON*RSH_EPSILON ? v : 0.0;
    *ow = fabs(w) > DBL_EPSILON*RSH_EPSILON ? w : 0.0;
}
