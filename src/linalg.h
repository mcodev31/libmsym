//
//  linalg.h
//  libmsym
//
//  Created by Marcus Johansson on 13/04/14.
//  Copyright (c) 2014 Marcus Johansson. 
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#ifndef __MSYM_LINALG_h
#define __MSYM_LINALG_h

void meye(double E[3][3]);
int vzero(const double v[3], double t);
int vparallel(const double v1[3], const double v2[3], double t);
int vperpendicular(const double v1[3], const double v2[3], double t);
double vnorm(double v[3]);
double vnorm2(const double v1[3],double v2[3]);
double vlnorm(int l, double *v);
double vlnorm2(int l, const double *v1, double *v2);
double vabs(const double v[3]);
double vlabs(int l, const double *v);
void vinv(double v[3]);
void vcopy(const double vi[3], double vo[3]);
void vlcopy(int l, const double *vi, double *vo);
void vcross(const double v1i[3], const double v2i[3], double vr[3]);
double vcrossnorm(const double[3],const double[3], double[3]);
double vdot(const double[3], const double[3]);
double vldot(int l, const double *v1, const double *v2);
int vequal(const double v1[3], const double v2[3], double t);
void vadd(const double[3], const double[3], double[3]);
void vladd(int l, const double *v1, const double *v2, double *vr);
void madd(const double A[3][3], const double B[3][3], double C[3][3]);
double vlsum(int l, double *v);
double vlsumsqr(int l, double *v);
void vsub(const double[3], const double[3], double[3]);
void vlsub(int l, const double *v1, const double *v2, double *vr);
void vscale(double, const double[3], double[3]);
void vlscale(double s,int l, const double *v, double *vr);
void mscale(double s,const double m[3][3], double mr[3][3]);
void vproj_plane(double v[3], double plane[3], double proj[3]);
void vproj(const double v[3], const double u[3], double vo[3]);
void vlproj(int l, const double *v, const double *u, double *vo);
void vcomplement(const double v1[3], double v2[3]);
double vangle(const double[3], const double[3]);
void vrotate(double theta, const double v[3], const double axis[3], double vr[3]);
void mrotate(double theta, const double axis[3], double m[3][3]);
void vreflect(const double v[3], const double axis[3], double vr[3]);
void mreflect(const double axis[3], double m[3][3]);
void mproject(double const v[3], double m[3][3]);
void mvmul(const double v[3], const double m[3][3], double r[3]);
void mmmul(const double A[3][3], const double B[3][3], double C[3][3]);
void minv(const double M[3][3], double I[3][3]);
double mdet(const double M[3][3]);
void mcopy(const double A[3][3], double B[3][3]);
int mequal(const double A[3][3], const double B[3][3], double t);
void malign(const double v[3], const double axis[3], double m[3][3]);
int ipow(int b, int e);

void jacobi(double m[6], double e[3], double ev[3][3], double threshold);

#ifndef __LIBMSYM_NO_VLA__
void mleye(int l, double (*E)[l]);
void mladd(int l, const double A[l][l], const double B[l][l], double C[l][l]);
void mlscale(double s,int l, const double m[l][l], double mr[l][l]);
void mvlmul(int r, int c, const double M[r][c], const double v[c], double vo[r]);
void mmlmul(int rla, int cla, const double A[rla][cla], int clb, const double B[cla][clb], double C[rla][clb]);
void mmtlmul(int rla, int cla, const double A[rla][cla], int rlb, const double B[rlb][cla], double C[rla][rlb]);
void mmlsymmul(int dim, const double m1[dim][dim], const double m2[dim][dim], double mr[dim][dim]);
double mltrace(int l, const double M[l][l]);
void mltranspose(int rl, int cl, const double A[rl][cl], double B[cl][rl]);
void mlcopy(int l, const double A[l][l], double B[l][l]);
int mgs(int l, const double M[l][l], double O[l][l], int n, double t);
int mgs2(int l, int lm, const double m[l][l], double o[l][l], int n, double t);
void kron(int al, const double A[al][al], int bl, const double B[bl][bl], int cl, double C[cl][cl]);
void kron2(int ar, int ac, const double A[ar][ac], int br, int bc, const double B[br][bc], double C[ar*br][ac*bc]);
void mlFilterSmall(int l, double A[l][l]);
#endif /* ifndef __LIBMSYM_NO_VLA__ */

#endif /* defined(__MSYM_LINALG_h) */
