//
//  orbital.c
//  libmsym
//
//  Created by Marcus Johansson on 07/11/14.
//  Copyright (c) 2014 Marcus Johansson. 
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdarg.h>

#include "msym.h"
#include "orbital.h"
#include "linalg.h"
#include "symop.h"
#include "permutation.h"
#include "point_group.h"
#include "subspace.h"

void printTransform(int r, int c, double M[r][c]);
void tabPrintTransform(int r, int c, double M[r][c],int indent);
void printSubspaceTree(CharacterTable *ct, msym_subspace_t *ss,int indent);
void tabprintf(char *format, int indent, ...);

msym_error_t getOrbitalSubspaceCoefficients(msym_subspace_t *ss, int basisl, msym_orbital_t basis[basisl], int *offset, double c[basisl][basisl]);
int filterSubspace(msym_subspace_t *ss);



//These are not the real orbitals but linear combinations:
//lm+ = p(z)((x+iy)^m+(x-iy)^m)
msym_error_t orbitalFromQuantumNumbers(int n, int l, int m, msym_orbital_t *o){
    
    if(l > n || abs(m) > l) goto err;
    
    o->n = n;
    o->l = l;
    o->m = m;
    
    memset(o->name,0,sizeof(o->name));
        
    switch(l) {
        case 0 :
            snprintf(o->name, sizeof(o->name), "%ds",n);
            //o->v = spolynomial[m+l];
            break;
        case 1 : {
            char *d = "?";
            //o->v = ppolynomial[m+l];
            switch(m) {
                case -1 : //Not really the case but makes things easier to generate
                    d = "y";
                    break;
                case 1 :
                    d = "x";
                    break;
                case 0 :
                    d = "z";
                    break;
            }
            snprintf(o->name, sizeof(o->name), "%dp%s",n,d);
            break;
        }
        case 2 : {
            //o->v = dpolynomial[m+l];
            char *d = (signbit(m) == 1 ? "-" : "+");
            snprintf(o->name, sizeof(o->name), "%dd%d%s",n,abs(m),d);
            break;
        }
        default : {
            char t = 'f' - 3 + l;
            char *d = (signbit(m) == 1 ? "-" : "+");
            snprintf(o->name, sizeof(o->name), "%d%c%d%s",n,t,abs(m),d);
        }
    }
    return MSYM_SUCCESS;
err:
    msymSetErrorDetails("Invalid orbital quantum numbers n:%d l:%d m:%d",n,l,m);
    return MSYM_INVALID_ORBITALS;
}

msym_error_t orbitalFromName(char *name, msym_orbital_t *o){
    int n, l, m;
    char cl, cm1, cm2;
    
    sscanf(name,"%d%c%c%c",&n,&cl,&cm1,&cm2);
    
    switch(cl) {
        case 's' : l = m = 0; break;
        case 'p' : {
            l = 1;
            switch(cm1) {
                case 'x' : m = 1; break;
                case 'y' : m = -1; break;
                case 'z' : m = 0; break;
                default : goto err;
            }
            break;
        }
        default :
            if(cl < 'd' || cl == 'e' || cl > 'z') goto err;
            l = cl == 'd' ? 2 : cl + 3 - 'f';
            m = (cm1 - (int)'0')*(cm2 == '-' ? -1 : 1);
    }

    return orbitalFromQuantumNumbers(n,l,m,o);
    
err:
    msymSetErrorDetails("Invalid orbital name %s",name);
    return MSYM_INVALID_ORBITALS;

}


msym_error_t orbitalPolynomial(int l, int m, double *poly){
    msym_error_t ret = MSYM_SUCCESS;
    int pdim = ipow(3,l);
    if(abs(m) > l) {ret = MSYM_INVALID_ORBITALS; goto err;}
    switch (l) {
        case 0 : vlcopy(pdim, spolynomial[m+l],poly); break;
        case 1 : vlcopy(pdim, ppolynomial[m+l],poly); break;
        case 2 : vlcopy(pdim, dpolynomial[m+l],poly); break;
        default: {
            msymSetErrorDetails("Cannot handle azimithal %d",l);
            ret = MSYM_INVALID_ORBITALS;
        }
    }
err:
    return ret;
}

//We can split this into a part for each l and just build the subspaces, and we already have the permutation matrix so
msym_error_t findProjection(CharacterTable *ct, int sopsl, msym_symmetry_operation_t sops[sopsl], msym_permutation_t perm[sopsl], int l, msym_orbital_t *basis[2*l+1]){
    msym_error_t ret = MSYM_SUCCESS;
    int kdim = ipow(3,l), setl = perm[0].p_length;
    double (*mkron)[kdim] = malloc(sizeof(double[kdim][kdim]));
    double (*mperm)[setl] = malloc(sizeof(double[setl][setl]));
    
    for(int m = 0; m < 2*l+1;m++){
        permutationMatrix(&perm[m], mperm);
    }
    
    free(mperm);
    free(mkron);
    return ret;
}

msym_error_t generateOrbitalTransforms(int sopsl, msym_symmetry_operation_t sops[sopsl], int l, double transform[sopsl][2*l+1][2*l+1]){
    msym_error_t ret = MSYM_SUCCESS;
    int kdim = ipow(3,l), norbs = 2*l+1;
    double (*mkron)[kdim][kdim] = malloc(sizeof(double[2][kdim][kdim]));
    double (*poly)[kdim] = malloc(sizeof(double[norbs][kdim]));
    
    for(int m = -l; m <= l;m++){
        if(MSYM_SUCCESS != (ret = orbitalPolynomial(l,m,poly[m+l]))) goto err;
        
        //Normalization
        vlnorm(kdim, poly[m+l]);
    }
    
    
    for(int i = 0;i < sopsl;i++){
        double M[3][3];
        mkron[0][0][0] = 1.0;
        symmetryOperationMatrix(&sops[i], M);
        for(int j = 0, d = 1;j < l;j++, d *= 3){
            kron(3,M,d,mkron[j%2],3*d,mkron[(j+1)%2]);
        }
        mmlmul(norbs, kdim, poly, kdim, mkron[l%2], mkron[(l+1)%2]);
        mmtlmul(norbs, kdim, mkron[(l+1)%2], norbs, poly, transform[i]);
        
        /* Scaling
        for(int j = 0; j < norbs;j++){
            double scale = vldot(kdim, poly[j], poly[j]);
            vlscale(1.0/scale, norbs, transform[i][j], transform[i][j]);
        }*/
    }
    
    
err:
    free(mkron);
    free(poly);
    return ret;
}

msym_error_t testSpan(msym_point_group_t *pg, int esl, msym_equivalence_set_t *es, msym_permutation_t **perm, int basisl, msym_orbital_t basis[basisl], msym_thresholds_t *thresholds, int *subspacel, msym_subspace_t **subspace, int **ospan){
    msym_error_t ret = MSYM_SUCCESS;
    
    int lmax = -1;
    for(int i = 0;i < basisl;i++){
        lmax = basis[i].l > lmax ? basis[i].l : lmax;
    }
    double (*bspan)[pg->ct->l] = calloc(lmax+1, sizeof(double[pg->ct->l])); // span of individual basis functions
    double (*pspan)[pg->ct->l] = calloc(esl, sizeof(double[pg->ct->l]));    // span of permutation operators
    double (*lspan)[pg->ct->l] = calloc(esl, sizeof(double[pg->ct->l]));    // total span of basis function on each ES
    double *rspan = calloc(pg->ct->l, sizeof(double));                      // decoposed total span of symmetrized basis
    double *mspan = calloc(pg->ct->l, sizeof(double));                      // memory for calculations
    
    int (*les)[lmax+1] = calloc(esl, sizeof(int[lmax+1]));                  // number of l-type basis functions in each ES
    
    if(lmax < 0){ret = MSYM_INVALID_ORBITALS; return ret;} //if we goto err here, code will get ugly due to scope
    
    /* determine number of l-type basis functions in each ES */
    for(int i = 0;i < esl;i++){
        msym_element_t *e = es[i].elements[0];
        for(int j = 0;j < e->aol;j++){
            les[i][e->ao[j]->l] += e->ao[j]->m == 0;
        }
    }
    
    /* calculate span of irreducible representations for basis functions and permutations */
    for(int s = 0; s < pg->sopsl;s++){
        for(int l = 0; l <= lmax;l++){
            for(int k = 0;k < pg->ct->l;k++){
                bspan[l][k] += pg->ct->irrep[k].v[pg->sops[s].cla]*symmetryOperationYCharacter(&pg->sops[s],l); //TODO: dynamic orbital/vibration
            }
        }
        for(int i = 0;i < esl;i++){
            int uma = 0;
            for(int j = 0; j < perm[i][s].c_length;j++) uma += perm[i][s].c[j].l == 1; //this is why we loop over irreps twice
            for(int k = 0;k < pg->ct->l;k++){
                pspan[i][k] += uma*pg->ct->irrep[k].v[pg->sops[s].cla];
            }
        }
    }
    
    /* scale and calculate total basis function span on each ES */
    for(int l = 0;l <= lmax;l++) vlscale(1.0/pg->order, pg->ct->l, bspan[l], bspan[l]);
    for(int i = 0;i < esl;i++){
        vlscale(1.0/pg->order, pg->ct->l, pspan[i], pspan[i]);
        for(int l = 0; l <= lmax;l++){
            //printf("lesl[%d][%d] = %d\n",i,l,lesl[i][l]);
            if(les[i][l] == 0) continue;
            vlscale((double) les[i][l], pg->ct->l, bspan[l], mspan);
            vladd(pg->ct->l, mspan, lspan[i], lspan[i]);
        }
    }
    /*
    for(int l = 0; l <= lmax;l++){
        printf("spherical harmonics span %d\n",l);
        for(int k = 0;k < pg->ct->l;k++) {
            int ssvl = (int)round(bspan[l][k]);
            if(ssvl) printf(" + %d%s",ssvl,pg->ct->irrep[k].name);
        }
        printf("\n");
    }
    
    for(int i = 0; i < esl;i++){
        printf("spherical harmonics total span %d\n",i);
        for(int k = 0;k < pg->ct->l;k++) {
            int ssvl = (int)round(lspan[i][k]);
            if(ssvl) printf(" + %d%s",ssvl,pg->ct->irrep[k].name);
        }
        printf("\n");
    }
    
    for(int i = 0; i < esl;i++){
        printf("permutation  span %d\n",i);
        for(int k = 0;k < pg->ct->l;k++) {
            int ssvl = (int)round(pspan[i][k]);
            if(ssvl) printf(" + %d%s",ssvl,pg->ct->irrep[k].name);
        }
        printf("\n");
    }*/
    
    /* calculate direct product of irreducible representations spanned by basis functions and permutations on each ES */
    for(int i = 0;i < esl;i++){
        for(int k = 0;k < pg->ct->l;k++){
            for(int j = 0;j < pg->ct->l && round(pspan[i][k]) > 0;j++){
                if(round(lspan[i][j]) == 0) continue;
                directProduct(pg->ct->l, &pg->ct->irrep[k], &pg->ct->irrep[j], mspan);
                
                //printf("%lf%s x %lf%s\n",pspan[i][k],pg->ct->irrep[k].name,lspan[i][j],pg->ct->irrep[j].name);
                //for(int i = 0;i < pg->ct->l;i++) printf(" | %lf", mspan[i]);
                //printf("\n");
                vlscale(pspan[i][k]*lspan[i][j], pg->ct->l, mspan, mspan);
                //for(int i = 0;i < pg->ct->l;i++) printf(" | %lf", mspan[i]);
                //printf("\n\n");
                
                vladd(pg->ct->l, mspan, rspan, rspan);
            }
        }
    }
    /* decompose direct product into irreducible representations */
    decomposeRepresentation(pg->ct, rspan, mspan);
    
    printf("decomposition of ");
    
    for(int i = 0;i < pg->ct->l;i++) printf(" + %lfx%d", rspan[i], i);
    printf(" = ");
    for(int i = 0;i < pg->ct->l;i++) printf(" + %d%s", (int) round(mspan[i]), pg->ct->irrep[i].name);
    printf("\n");
    
    
    /*
    for(int k = 0;k < pg->ct->l;k++){
        iss[ssl].type = MASS_WEIGHTED_COORDINATES;
        iss[ssl].irrep = k;
        for(int i = 0;i < esl;i++){
            for(int s = 0; s < pg->sopsl;s++){
                int uma = 0;
                for(int j = 0; j < perm[i][s].c_length;j++) uma += perm[i][s].c[j].l == 1;
                pspan[i][k] += uma*pg->ct->irrep[k].v[pg->sops[s].cla];
            }
            
            int ssvl = round(pg->ct->irrep[k].d*pspan[i][k]/pg->order);
            iss[ssl].subspacel += (ssvl > 0);
            ispan[i][k] += ssvl;
            aspan[k] += ssvl;
        }
        
        if(iss[ssl].subspacel){
            iss[ssl].subspace = calloc(iss[ssl].subspacel, sizeof(msym_subspace_t));
            ssl++;
        }
        
        printf("permutation span %s = %d\n",pg->ct->irrep[k].name,aspan[k]);
        for(int i = 0; i < esl;i++) printf("\tspan[%d] %s = %d\n",i,pg->ct->irrep[k].name,ispan[i][k]);
        
    }*/
    
    return ret;
}

msym_error_t generateOrbitalSubspaces(msym_point_group_t *pg, int esl, msym_equivalence_set_t *es, msym_permutation_t **perm, int basisl, msym_orbital_t basis[basisl], msym_thresholds_t *thresholds, int *subspacel, msym_subspace_t **subspace, int **pspan){
    msym_error_t ret = MSYM_SUCCESS;
    
    testSpan(pg,esl,es,perm,basisl,basis,thresholds,subspacel,subspace,pspan);
    exit(1);
    
    int lmax = -1, nmax = -1, eslmax = -1;
    for(int i = 0;i < basisl;i++){
        lmax = basis[i].l > lmax ? basis[i].l : lmax;
        nmax = basis[i].n > nmax ? basis[i].n : nmax;
    }
    
    if(lmax < 0){ret = MSYM_INVALID_ORBITALS; return ret;} //if we goto err here, code will get ugly due to scope
    
    for(int i = 0;i < esl;i++) eslmax = es[i].length > eslmax ? es[i].length : eslmax;
    
    struct _ltransforms {int d; void *t;} *lts = calloc(lmax+1,sizeof(struct _ltransforms));
    double (*mkron)[(2*lmax+1)*pg->order] = malloc(sizeof(double[(2*lmax+1)*pg->order][(2*lmax+1)*pg->order]));
    double (*mperm)[pg->order] = malloc(sizeof(double[pg->order][pg->order]));
    
    double (*mproj)[pg->ct->l+1][(2*lmax+1)*pg->order][(2*lmax+1)*pg->order] = malloc(sizeof(double[lmax+1][pg->ct->l+1][(2*lmax+1)*pg->order][(2*lmax+1)*pg->order]));
    double (*lspan)[pg->ct->l] = malloc(sizeof(double[lmax+1][pg->ct->l]));
    int (*ispan) = calloc(pg->ct->l,sizeof(int));
    int *aspan = calloc(pg->ct->l,sizeof(int));
    int *nl = malloc(sizeof(int[lmax+1]));
    
    msym_orbital_t *(*omap)[nmax][lmax][2*lmax+1] = malloc(sizeof(msym_orbital_t *[eslmax][nmax+1][lmax+1][2*lmax+1]));
    
    *subspace = NULL;
    
    msym_subspace_t *iss = calloc(pg->ct->l, sizeof(msym_subspace_t));
    for(int k = 0;k < pg->ct->l;k++){
        iss[k].type = ATOMIC_ORBITAL;
        iss[k].irrep = k;
        iss[k].subspacel = esl;
        iss[k].subspace = calloc(esl, sizeof(msym_subspace_t));
    }
    
    for(int l = 0; l <= lmax;l++){
        lts[l].d = 2*l+1;
        lts[l].t = malloc(sizeof(double[pg->sopsl][lts[l].d][lts[l].d]));
        if(MSYM_SUCCESS != (ret = generateOrbitalTransforms(pg->sopsl, pg->sops, l, lts[l].t))) goto err;
    }
    
    for(int i = 0; i < esl;i++){
        int esilmax = -1, esinmax = -1;
        
        memset(nl,0,sizeof(int[lmax+1]));
        for(int j = 0;j < es[i].elements[0]->aol;j++){
            esilmax = esilmax < es[i].elements[0]->ao[j]->l ? es[i].elements[0]->ao[j]->l : esilmax;
            esinmax = esinmax < es[i].elements[0]->ao[j]->n ? es[i].elements[0]->ao[j]->n : esinmax;
            nl[es[i].elements[0]->ao[j]->l] += es[i].elements[0]->ao[j]->m == 0;
        }
        
        msym_orbital_t *(*esomap)[esinmax+1][esilmax+1][2*esilmax+1] = omap;
        
        memset(esomap,0,sizeof(msym_orbital_t *[es->length][esinmax+1][esilmax+1][2*esilmax+1]));
        for(int a = 0;a < es[i].length;a++){
            for(int ao = 0;ao < es[i].elements[a]->aol;ao++){
                msym_orbital_t *o = es[i].elements[a]->ao[ao];
                esomap[a][o->n][o->l][o->m+o->l] = o;
            }
        }
        
        memset(lspan,0,sizeof(double[lmax+1][pg->ct->l]));
        
        for(int l = 0;l <= esilmax;l++){
            int d = es[i].length*lts[l].d;
            double (*mlproj)[d][d] = mproj[l];
            memset(mlproj,0,sizeof(double[pg->ct->l][d][d]));
            memset(ispan,0,sizeof(int[pg->ct->l]));
            
            for(int s = 0;s < pg->sopsl;s++){
                double (*lt)[lts[l].d][lts[l].d] = lts[l].t;
                permutationMatrix(&perm[i][s], mperm);
                kron(perm[i][s].p_length,mperm,lts[l].d,lt[s],d,mkron);
                
                char buf[16];
                symmetryOperationShortName(&pg->sops[s], 16, buf);
                printf("%s%d=\n",buf,s);
                //printSymmetryOperation(&pg->sops[s]);
                
                printTransform(d, d, mkron);
                
                printf("Sops = union(Sops,(%s%d,))\n",buf,s);
                
                for(int k = 0;k < pg->ct->l;k++){
                    lspan[l][k] += pg->ct->irrep[k].v[pg->sops[s].cla]*mltrace(d, mkron);
                    mlscale(pg->ct->irrep[k].v[pg->sops[s].cla], d, mkron, mlproj[pg->ct->l]); //mproj[pg->ct.l] is a workspace
                    mladd(d, mlproj[pg->ct->l], mlproj[k], mlproj[k]); //Could do this based on the span later, but it's not a huge amount of work
                }
            }
            memset(mlproj[pg->ct->l],0,sizeof(double[d][d]));
            int nirrepl = 0;
            for(int k = 0;k < pg->ct->l;k++){
                int lirrepl = nirrepl;
                msym_subspace_t *ss = &iss[k].subspace[i];
                ispan[k] = (((int)round(lspan[l][k]))/pg->order);
                mlscale(((double) pg->ct->irrep[k].d)/pg->order, d, mlproj[k], mlproj[k]);
                nirrepl = mgs(d, mlproj[k], mlproj[pg->ct->l], nirrepl, thresholds->orthogonalization/basisl);
                //printTransform(d, d, mlproj[pg->ct->l]);
                if(nirrepl - lirrepl != ispan[k]*pg->ct->irrep[k].d){
                    //printTransform(d, d, mlproj[k]);
                    ret = MSYM_SUBSPACE_ERROR;
                    msymSetErrorDetails("Ortogonal subspace of dimension (%d) inconsistent with span (%d) in %s",nirrepl - lirrepl,ispan[k]*pg->ct->irrep[k].d,pg->ct->irrep[k].name);
                    goto err;
                    
                }
                
                ss->type = ATOMIC_ORBITAL;
                ss->irrep = k;
                
                if(ispan[k] > 0){
                    ss->subspace = realloc(ss->subspace, sizeof(msym_subspace_t[ss->subspacel+nl[l]]));
                    for(int n = l+1; n <= esinmax;n++){
                        if(esomap[0][n][l][0] == NULL) continue;
                        
                        aspan[k] += ispan[k]*pg->ct->irrep[k].d;
                        msym_subspace_t *nss = &ss->subspace[ss->subspacel];
                        memset(nss,0,sizeof(msym_subspace_t));

                        nss->type = ATOMIC_ORBITAL;
                        nss->irrep = ss->irrep;
                        nss->d = ispan[k]*pg->ct->irrep[k].d;
                        double (*space)[d] = malloc(sizeof(double[nss->d][d]));
                        for(int dim = 0; dim < ss->subspace[ss->subspacel].d;dim++){
                            vlnorm2(d, mlproj[pg->ct->l][lirrepl+dim], space[dim]);
                        }
                        nss->space = (double*) space;
                        nss->basisl = 0;
                        nss->basis.o = malloc(sizeof(msym_orbital_t *[d]));
                        for(int e = 0;e < es[i].length;e++){
                            for(int m = -l;m <= l;m++){
                                nss->basis.o[nss->basisl++] = esomap[e][n][l][m+l];
                            }
                        }
                        ss->subspacel++;
                        if(nss->basisl != d) {
                            ret = MSYM_SUBSPACE_ERROR;
                            msymSetErrorDetails("Basis length (%d) inconsistent with projection operator dimensions (%d)",nss->basisl,d);
                            goto err;
                        }
                    }
                }
            }
        }
    }

    printf("Subspace span (vectors) = ");
    for(int k = 0;k < pg->ct->l;k++){
        printf(" + %d%s",aspan[k],pg->ct->irrep[k].name);
    }
    printf("\n");
    
    msym_subspace_t tss = {.subspacel = pg->ct->l, .subspace = iss, .d = 0, .basisl = 0, .space = NULL};
    filterSubspace(&tss);
    *subspace = tss.subspace;
    *subspacel = tss.subspacel;
    *pspan = aspan;
    
    printSubspace(pg->ct, &tss);
    
    free(omap);
    free(nl);
    free(ispan);
    free(lspan);
    free(mproj);
    free(mperm);
    free(mkron);
    for(int l = 0;l <= lmax;l++){
        free(lts[l].t);
    }
    free(lts);
    
    return ret;
    
err:
    free(aspan);
    free(omap);
    free(nl);
    free(ispan);
    free(lspan);
    free(mproj);
    free(mperm);
    free(mkron);
    for(int l = 0;l < lmax;l++){
        free(lts[l].t);
    }
    free(lts);
    for(int k = 0;k < pg->ct->l;k++){
        freeSubspace(&iss[k]);
    }
    free(iss);
    return ret;
}

msym_error_t getOrbitalSubspaces(int ssl, msym_subspace_t ss[ssl], int basisl, msym_orbital_t basis[basisl], double c[basisl][basisl]){
    msym_error_t ret = MSYM_SUCCESS;
    int index = 0;
    memset(c,0,sizeof(double[basisl][basisl]));
    for(int i = 0;i < ssl;i++){
        if(MSYM_SUCCESS != (ret = getOrbitalSubspaceCoefficients(&ss[i],basisl,basis,&index,c))) goto err;
    }
    
    if(index != basisl){
        msymSetErrorDetails("Subspace index (%d) does not match basis length (%d)",index,basisl);
        ret = MSYM_INVALID_SUBSPACE;
        goto err;
    }
    
    return ret;
err:
    return ret;
}

msym_error_t getOrbitalSubspaceCoefficients(msym_subspace_t *ss, int basisl, msym_orbital_t basis[basisl], int *offset, double c[basisl][basisl]){
    msym_error_t ret = MSYM_SUCCESS;
    
    int index = *offset;
    if(index >= basisl) {
        msymSetErrorDetails("Subspace index (%d) is larger than basis length (%d)",index,basisl);
        ret = MSYM_INVALID_SUBSPACE;
        goto err;
    }
    
    if(ss->subspacel == 0){
        double (*space)[ss->basisl] = (double (*)[ss->basisl]) ss->space;
        if(index+ss->d > basisl) {
            msymSetErrorDetails("Generated subspaces (%d) is larger than basis length (%d)",index+ss->d,basisl);
            ret = MSYM_INVALID_SUBSPACE;
            goto err;
        }
        for(int d = 0;d < ss->d;d++){
            for(int b = 0;b < ss->basisl;b++){
                c[index][ss->basis.o[b]-basis] = space[d][b];
            }
            //printf("orbital %d is in irrep %d\n",index,ss->irrep);
            index++;
        }
    } else {
        for(int i = 0;i < ss->subspacel;i++){
            if(MSYM_SUCCESS != (ret = getOrbitalSubspaceCoefficients(&ss->subspace[i],basisl,basis,&index,c))) goto err;
        }
    }
    
    *offset = index;
    
    return ret;
err:
    return ret;
}

//Density matrix without occupation numbers
void densityMatrix(int l, double M[l][l], double D[l][l]){
    memset(D,0,sizeof(double[l][l]));
    for(int i = 0; i < l;i++){
        for(int j = 0;j < l;j++){
            for(int k = 0;k < l;k++){
                D[i][j] += M[k][i]*M[k][j];
            }
        }
    }
}

void printOrbital(msym_orbital_t *orb){
    printf("Orbital(%d,%d,%d) : %s\n",orb->n, orb->l, orb->m, orb->name);
}



