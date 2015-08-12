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
        
        
        ////
        char buf[16];
        symmetryOperationShortName(&sops[i], 16, buf);
        printf("%s%d=\n",buf,i);
        
        printTransform(norbs, norbs, transform[i]);
        
        printf("Sops = union(Sops,(%s%d,))\n",buf,i);
        /////
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




msym_error_t testSpan2(msym_point_group_t *pg, int esl, msym_equivalence_set_t *es, msym_permutation_t **perm, int basisl, msym_basis_function_t basis[basisl], msym_element_t *elements, msym_equivalence_set_t **esmap, msym_thresholds_t *thresholds, int *subspacel, msym_subspace_2_t **subspace){
    msym_error_t ret = MSYM_SUCCESS;
    
    int lmax = -1, nmax = 0;
    msym_basis_type_2_t ftype = basis[0].type;
    for(int i = 0;i < basisl;i++){
        if(basis[i].type != ftype) {nmax = -1; break;}
        lmax = basis[i].f.sh.l > lmax ? basis[i].f.sh.l : lmax;
        nmax = basis[i].f.sh.n > nmax ? basis[i].f.sh.n : nmax;
    }
    
    if(lmax < 0 || nmax < 1){
        if(nmax == -1) msymSetErrorDetails("Basis functions are not of the same type");
        else msymSetErrorDetails("Invalid sperical harmonics quantum numbers");
        ret = MSYM_INVALID_ORBITALS;
        return ret;
    } else if (ftype != SPHERICAL_HARMONIC){
        msymSetErrorDetails("Basis function type not supported");
        ret = MSYM_INVALID_ORBITALS;
        return ret;
    }
    
    int projm = (2*lmax+1)*pg->order;
    
    double (*bspan)[pg->ct2->d] = calloc(lmax+1, sizeof(double[pg->ct2->d]));     // span of individual basis functions
    double (*pspan)[pg->ct2->d] = calloc(esl, sizeof(double[pg->ct2->d]));        // span of permutation operators
    double (*lspan)[pg->ct->l] = calloc(esl, sizeof(double[pg->ct->l]));        // total span of basis function on each ES
    double *rspan = calloc(pg->ct->l, sizeof(double));                          // total direct product symmetrized basis
    double *dspan = calloc(pg->ct->l, sizeof(double));                          // decoposed total span of symmetrized basis
    double *mspan = calloc(pg->ct->l, sizeof(double));                          // span decomposition memory
    double (*mproj)[projm] = calloc(projm, sizeof(double[projm]));              // projection operator memory
    double (*mscal)[projm] = calloc(projm, sizeof(double[projm]));              // projection scaling memory
    
    
    struct _ltransforms {int d; void *t;} *lts = calloc(lmax+1,sizeof(struct _ltransforms)); // transformation matrices for basis functions
    
    int (*les)[lmax+1] = calloc(esl, sizeof(int[lmax+1]));                      // number of l-type basis functions in each ES
    
    msym_basis_function_t dbf = {.type = ftype};
    double (*ctable)[pg->ct2->d] = pg->ct2->table;
    
    /* determine number of l-type basis functions in each ES */
    for(int o = 0;o < basisl;o++){
        les[esmap[basis[o].element - elements] - es][basis[o].f.sh.l] += basis[o].f.sh.m == 0;
    }
    
    
    int (*lesold)[lmax+1] = calloc(esl, sizeof(int[lmax+1]));
    for(int i = 0;i < esl;i++){
        msym_element_t *e = es[i].elements[0];
        for(int j = 0;j < e->aol;j++){
            lesold[i][e->ao[j]->l] += e->ao[j]->m == 0;
        }
    }

    
    /* calculate span of irreducible representations for basis functions and permutations */
    for(int s = 0; s < pg->sopsl;s++){
        for(int l = 0; l <= lmax;l++){
            dbf.f.sh.l = l; //TODO: function type based
            double c = symmetryOperationCharacter(&pg->sops[s], &dbf);
            printSymmetryOperation(&pg->sops[s]);
            for(int k = 0;k < pg->ct2->d;k++){
                bspan[l][k] += c*ctable[k][pg->sops[s].cla];
            }
        }
        for(int i = 0;i < esl;i++){
            int uma = 0;
            for(int j = 0; j < perm[i][s].c_length;j++) uma += perm[i][s].c[j].l == 1; //this is why we loop over irreps twice
            for(int k = 0;k < pg->ct2->d;k++){
                pspan[i][k] += uma*ctable[k][pg->sops[s].cla];
            }
        }
    }
    
    /* scale basis span, calculate basis function transforms and symmetry species vectors */
    for(int l = 0;l <= lmax;l++){
        int d = 2*l+1;
        double (*lproj)[d] = mproj;
        double (*lscal)[d] = mscal;
        
        vlscale(1.0/pg->order, pg->ct2->d, bspan[l], bspan[l]);
        
        lts[l].d = d;
        lts[l].t = malloc(sizeof(double[pg->order+1][d][d]));
        
        double (*st)[d][d] = lts[l].t;
        
        if(MSYM_SUCCESS != (ret = generateOrbitalTransforms(pg->order, pg->sops, l, lts[l].t))) goto err; //TODO: generalize basis function concept
        
        memset(st[pg->order], 0, sizeof(double[d][d]));
        
        for(int k = 0, oirl = 0, nirl = 0;k < pg->ct2->d;k++, oirl = nirl){
            
            int vspan = pg->ct2->s[k].d*((int) round(bspan[l][k]));
            if(vspan == 0) continue;
            
            memset(lproj, 0, sizeof(double[d][d]));
            for(int s = 0;s < pg->order;s++){
                mlscale(ctable[k][pg->sops[s].cla], d, st[s], lscal);
                mladd(d, lscal, lproj, lproj);
            }
            
            mlscale(((double) pg->ct2->s[k].d)/pg->order, d, lproj, lproj);
            nirl = mgs(d, lproj, st[pg->order], oirl, thresholds->orthogonalization/basisl);
            
            if(nirl - oirl != vspan){
                printf("bspan[%d][%d] = %lf\n",l,k,bspan[l][k]);
                ret = MSYM_SUBSPACE_ERROR;
                msymSetErrorDetails("Ortogonal subspace of dimension (%d) inconsistent with span (%d) in %s",nirl - oirl,vspan,pg->ct2->s[k].name);
                goto err;
                
            }
        }
        
        for(int i = 0; i < d;i++) vlnorm(d, st[pg->order][i]);
        
        ///////
        for(int i = 0;i < pg->ct2->d;i++) printf(" + %d%s", (int) round(bspan[l][i]), pg->ct2->s[i].name);
        printf("\n");
        printTransform(d, d, st[pg->order]);
        ///////
        
    }
    
    /* scale permutation span and calculate total basis function span on each ES */
    for(int i = 0;i < esl;i++){
        vlscale(1.0/pg->order, pg->ct->l, pspan[i], pspan[i]);
        for(int l = 0; l <= lmax;l++){
            //printf("lesl[%d][%d] = %d\n",i,l,lesl[i][l]);
            if(les[i][l] == 0) continue;
            les[i][l] /= es[i].length;
            vlscale((double) les[i][l], pg->ct->l, bspan[l], mspan);
            vladd(pg->ct->l, mspan, lspan[i], lspan[i]);
        }
    }
    
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
    }
    
    for(int i = 0;i < esl;i++){
        for(int l = 0; l <= lmax;l++) {
            if(les[i][l] != lesold[i][l]){
                printf("les calc error %d != %d\n",les[i][l], lesold[i][l]);
                exit(1);
            }
        }
    }
    
    /* calculate direct product of irreducible representations spanned by basis functions and permutations on each ES (don't really need to do this) */
    for(int i = 0;i < esl;i++){
        for(int k = 0;k < pg->ct2->d;k++){
            for(int j = 0;j < pg->ct2->d && round(pspan[i][k]) > 0;j++){
                if(round(lspan[i][j]) == 0) continue;
                printf("need direct product calc");
                exit(1);
                directProduct(pg->ct2->d, &pg->ct->irrep[k], &pg->ct->irrep[j], mspan);
                
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
    decomposeRepresentation(pg->ct, rspan, dspan);
    
    printf("decomposed into\n");
    for(int prk = 0;prk < pg->ct->l;prk++) printf(" + %d%s", (int) round(dspan[prk]), pg->ct->irrep[prk].name);
    printf("\n");
    int ssl = 0;


err:
    return ret;
}

msym_error_t testSpanConvert(msym_point_group_t *pg, int esl, msym_equivalence_set_t *es, msym_permutation_t **perm, int basisl, msym_orbital_t basis[basisl], msym_thresholds_t *thresholds, int *subspacel, msym_subspace_t **subspace, int **ospan){
    
    msym_character_table_2_t ct;
    ct.d = pg->ct->l;
    double t[ct.d][ct.d];
    msym_symmetry_species_2_t ssp[ct.d];
    ct.s = ssp;
    ct.table = t;
    for(int i = 0;i < pg->ct->l;i++){
        memcpy(t[i], pg->ct->irrep[i].v, sizeof(double[ct.d]));
        ssp[i].d = pg->ct->irrep[i].d;
        sprintf(ssp[i].name, "%s",pg->ct->irrep[i].name);
    }
    
    pg->ct2 = &ct;
    double (*ctable)[pg->ct2->d] = pg->ct2->table;
    
    for(int i = 0;i < pg->ct->l;i++){
        for(int j = 0;j < pg->ct->l;j++){
            printf("%lf ",ctable[i][j]);
        }
        printf("\n");
    }
    
    msym_basis_function_t *bs = calloc(basisl, sizeof(msym_basis_function_t));
    

    
    for(int i = 0;i < basisl;i++){
        bs[i].type = SPHERICAL_HARMONIC;
        sprintf(bs[i].name, "%s",basis[i].name);
        bs[i].f.sh.n = basis[i].n;
        bs[i].f.sh.l = basis[i].l;
        bs[i].f.sh.m = basis[i].m;
        for(int j = 0;j < esl;j++){
            for(int k = 0; k < es[j].length;k++){
                msym_element_t *e = es[j].elements[k];
                for(int l = 0;l < e->aol;l++){
                    if(e->ao[l] == &basis[i]){
                        bs[i].element = e;
                        printf("bs[%d] = es[%d].elements[%d]->ao[%d]\n",i,j,k,l);
                        break;
                    }
                }
            }
        }
    }
    
    msym_element_t *start = NULL;
    
    int nelement = 0;
    
    
    for(int i = 0;i < esl;i++){
        for(int j = 0;j < es[i].length;j++){
            if(start == NULL || es[i].elements[j] < start) start = es[i].elements[j];
        }
        nelement += es[i].length;
    }
    
    msym_equivalence_set_t **esmap = calloc(nelement, sizeof(msym_equivalence_set_t *));
    
    for(int i = 0;i < esl;i++){
        for(int j = 0;j < es[i].length;j++){
            if(start == NULL || es[i].elements[j] < start) start = es[i].elements[j];
        }
    }
    
    for(int i = 0;i < esl;i++){
        for(int j = 0;j < es[i].length;j++){
            esmap[es[i].elements[j] - start] = &es[i];
        }
    }

    
    msym_subspace_2_t *ss = NULL;
    
    testSpan2(pg, esl, es, perm, basisl, bs, start, esmap, thresholds, subspacel, &ss);
    
    return MSYM_INVALID_CHARACTER_TABLE;
}

msym_error_t testSpan(msym_point_group_t *pg, int esl, msym_equivalence_set_t *es, msym_permutation_t **perm, int basisl, msym_orbital_t basis[basisl], msym_thresholds_t *thresholds, int *subspacel, msym_subspace_t **subspace, int **ospan){
    msym_error_t ret = MSYM_SUCCESS;
    
    int lmax = -1, nmax = -1;
    for(int i = 0;i < basisl;i++){
        lmax = basis[i].l > lmax ? basis[i].l : lmax;
        nmax = basis[i].n > nmax ? basis[i].n : nmax;
    }
    
    if(lmax < 0){ret = MSYM_INVALID_ORBITALS; return ret;} //if we goto err here, code will get ugly due to scope
    
    int projm = (2*lmax+1)*pg->order;
    double (*bspan)[pg->ct->l] = calloc(lmax+1, sizeof(double[pg->ct->l]));     // span of individual basis functions
    double (*pspan)[pg->ct->l] = calloc(esl, sizeof(double[pg->ct->l]));        // span of permutation operators
    double (*lspan)[pg->ct->l] = calloc(esl, sizeof(double[pg->ct->l]));        // total span of basis function on each ES
    double *rspan = calloc(pg->ct->l, sizeof(double));                          // total direct product symmetrized basis
    double *dspan = calloc(pg->ct->l, sizeof(double));                          // decoposed total span of symmetrized basis
    double *mspan = calloc(pg->ct->l, sizeof(double));                          // span decomposition memory
    double (*mproj)[projm] = calloc(projm, sizeof(double[projm]));              // projection operator memory
    double (*mscal)[projm] = calloc(projm, sizeof(double[projm]));              // projection scaling memory
    double (*mperm)[pg->order] = calloc(pg->sopsl, sizeof(double[pg->order]));  // permutation memory
    double (*morth)[pg->order] = calloc(pg->sopsl, sizeof(double[pg->order]));  // permutation orthoginalization memory
    double (*mbasis)[projm] = calloc(basisl, sizeof(double[projm]));            // basis function coefficients
    double (*mdec)[projm] = calloc(basisl, sizeof(double[projm]));              // directo product decomposition memory
    int *ispan = calloc(pg->ct->l, sizeof(int));                                // integer decoposed total span of symmetrized basis

    
    
    int (*les)[lmax+1] = calloc(esl, sizeof(int[lmax+1]));                      // number of l-type basis functions in each ES
    
    msym_orbital_t *(*omap)[nmax][lmax][2*lmax+1] = malloc(sizeof(msym_orbital_t *[pg->order][nmax+1][lmax+1][2*lmax+1]));
    
    msym_subspace_t *iss = NULL;
    int ssi = 0;
    
    
    struct _ltransforms {int d; void *t;} *lts = calloc(lmax+1,sizeof(struct _ltransforms)); // transformation matrices for basis functions
    
    
    
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
            printf("l=%d -> %lf ",l, symmetryOperationYCharacter(&pg->sops[s],l));
            printSymmetryOperation(&pg->sops[s]);
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
    
    /* scale basis span, calculate basis function transforms and symmetry species vectors */
    for(int l = 0;l <= lmax;l++){
        int d = 2*l+1;
        double (*lproj)[d] = mproj;
        double (*lscal)[d] = mscal;
        
        vlscale(1.0/pg->order, pg->ct->l, bspan[l], bspan[l]);
        
        lts[l].d = d;
        lts[l].t = malloc(sizeof(double[pg->sopsl+1][d][d]));
        
        double (*st)[d][d] = lts[l].t;
        
        if(MSYM_SUCCESS != (ret = generateOrbitalTransforms(pg->sopsl, pg->sops, l, lts[l].t))) goto err; //TODO: generalize basis function concept
        
        memset(st[pg->sopsl], 0, sizeof(double[d][d]));
        
        for(int k = 0, oirl = 0, nirl = 0;k < pg->ct->l;k++, oirl = nirl){
            int vspan = pg->ct->irrep[k].d*((int) round(bspan[l][k]));
            if(vspan == 0) continue;
            
            memset(lproj, 0, sizeof(double[d][d]));
            for(int s = 0;s < pg->sopsl;s++){
                mlscale(pg->ct->irrep[k].v[pg->sops[s].cla], d, st[s], lscal);
                mladd(d, lscal, lproj, lproj);
            }
            
            mlscale(((double) pg->ct->irrep[k].d)/pg->order, d, lproj, lproj);
            nirl = mgs(d, lproj, st[pg->sopsl], oirl, thresholds->orthogonalization/basisl);
            
            if(nirl - oirl != vspan){
                printf("bspan[%d][%d] = %lf\n",l,k,bspan[l][k]);
                ret = MSYM_SUBSPACE_ERROR;
                msymSetErrorDetails("Ortogonal subspace of dimension (%d) inconsistent with span (%d) in %s",nirl - oirl,vspan,pg->ct->irrep[k].name);
                goto err;
                
            }
        }
        
        for(int i = 0; i < d;i++) vlnorm(d, st[pg->sopsl][i]);
        
        ///////
        for(int i = 0;i < pg->ct->l;i++) printf(" + %d%s", (int) round(bspan[l][i]), pg->ct->irrep[i].name);
        printf("\n");
        printTransform(d, d, st[pg->sopsl]);
        ///////
        
    }
    
    /* scale permutation span and calculate total basis function span on each ES */
    for(int i = 0;i < esl;i++){
        vlscale(1.0/pg->order, pg->ct->l, pspan[i], pspan[i]);
        for(int l = 0; l <= lmax;l++){
            //printf("lesl[%d][%d] = %d\n",i,l,lesl[i][l]);
            if(les[i][l] == 0) continue;
            vlscale((double) les[i][l], pg->ct->l, bspan[l], mspan);
            vladd(pg->ct->l, mspan, lspan[i], lspan[i]);
        }
    }
    
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
    }
    
    /* calculate direct product of irreducible representations spanned by basis functions and permutations on each ES (don't really need to do this) */
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
    decomposeRepresentation(pg->ct, rspan, dspan);
    
    printf("decomposed into\n");
    for(int prk = 0;prk < pg->ct->l;prk++) printf(" + %d%s", (int) round(dspan[prk]), pg->ct->irrep[prk].name);
    printf("\n");
    int ssl = 0;
    
    // Allocate basis length and realloc?
    for(int k = 0;k < pg->ct->l;k++) ssl += dspan[k]; // count number of subspaces
    printf("resulting in %d total subspaces\n",ssl);
    iss = calloc(ssl, sizeof(msym_subspace_t));
    
    for(int i = 0; i < esl; i++){
        int d = es[i].length;
        double (*pproj)[d] = mproj;
        double (*pscal)[d] = mscal;
        double (*porth)[d] = morth;
        
        memset(porth, 0, sizeof(double[d][d]));
        
        int esilmax = -1, esinmax = -1;
        
        for(int j = 0;j < es[i].elements[0]->aol;j++){
            esilmax = esilmax < es[i].elements[0]->ao[j]->l ? es[i].elements[0]->ao[j]->l : esilmax;
            esinmax = esinmax < es[i].elements[0]->ao[j]->n ? es[i].elements[0]->ao[j]->n : esinmax;
        }
        
        /* orbital mapping for ES */
        msym_orbital_t *(*esomap)[esinmax+1][esilmax+1][2*esilmax+1] = omap;
        memset(esomap,0,sizeof(msym_orbital_t *[es->length][esinmax+1][esilmax+1][2*esilmax+1]));
        for(int a = 0;a < es[i].length;a++){
            for(int ao = 0;ao < es[i].elements[a]->aol;ao++){
                msym_orbital_t *o = es[i].elements[a]->ao[ao];
                esomap[a][o->n][o->l][o->m+o->l] = o;
            }
        }
        
        for(int k = 0, oirl = 0, nirl = 0;k < pg->ct->l;k++, oirl = nirl){
            int vspan = pg->ct->irrep[k].d*((int) round(pspan[i][k]));
            if(vspan == 0) continue;
            
            memset(pproj, 0, sizeof(double[d][d]));
            for(int s = 0;s < pg->sopsl;s++){
                if(pg->ct->irrep[k].v[pg->sops[s].cla] == 0) continue;
                permutationMatrix(&perm[i][s], mperm);
                mlscale(pg->ct->irrep[k].v[pg->sops[s].cla], d, mperm, pscal);
                mladd(d, pscal, pproj, pproj);
            }
            
            mlscale(((double) pg->ct->irrep[k].d)/pg->order, d, pproj, pproj);
            nirl = mgs(d, pproj, porth, oirl, thresholds->orthogonalization/basisl);
            
            if(nirl - oirl != vspan){
                ret = MSYM_SUBSPACE_ERROR;
                msymSetErrorDetails("Ortogonal ES subspace of dimension (%d) inconsistent with span (%d) in %s",nirl - oirl,vspan,pg->ct->irrep[k].name);
                goto err;
                
            }
            
            for(int oi = oirl; oi < nirl;oi++) vlnorm(d, porth[oi]);
            
            for(int l = 0;l <= lmax;l++){
                if(les[i][l] <= 0) continue;
                int li = 0, ld = lts[l].d;
                double (*lst)[ld][ld] = lts[l].t;
                for(int lk = 0;lk < pg->ct->l;lk++){
                    int lvspan = pg->ct->irrep[lk].d*((int) round(bspan[l][lk])), dd = d*ld;
                    if(lvspan == 0) continue;
                    kron2(vspan, d, &porth[oirl], lvspan, ld, &lst[pg->sopsl][li], mbasis);
                    li += (int) round(bspan[l][lk]);
                    printf("%d-dimensional %s part of permutaion and %d-dimensional %s part of l=%d will produce lcao:\n",vspan,pg->ct->irrep[k].name,lvspan,pg->ct->irrep[lk].name,l);
                    
                    printf("%d %lf %p %p %p\n",li, (&lst[pg->sopsl][li])[0][0],&lst[pg->sopsl][li],&(lst[pg->sopsl][li]),lst[pg->sopsl]);
                    printTransform(ld,ld,lst[pg->sopsl]);
                    printTransform(vspan, d, &porth[oirl]);
                    printTransform(lvspan, ld, &(lst[pg->sopsl][li]));
                    printTransform(vspan*lvspan, dd, mbasis);
                    
                    directProduct(pg->ct->l, &pg->ct->irrep[k], &pg->ct->irrep[lk], rspan);
                    vlscale(pspan[i][k]*bspan[l][lk], pg->ct->l, rspan, rspan);
                    decomposeRepresentation(pg->ct, rspan, mspan);
                    if(pg->ct->irrep[k].d > 1 && pg->ct->irrep[lk].d > 1){
                        memcpy(mdec, mbasis, sizeof(double[vspan*lvspan][dd]));
                        //printf("mbasis\n");
                        //printTransform(vspan*lvspan, dd, mbasis);
                        memset(mbasis, 0, sizeof(double[vspan*lvspan][dd]));
                        printf("This direct product needs to be reduced\n");
                        for(int dk = 0, doirl = 0, dnirl = 0;dk < pg->ct->l;dk++, doirl = dnirl){
                            int dvspan = pg->ct->irrep[dk].d*((int) round(mspan[dk]));
                            if(dvspan == 0) continue;
                            
                            double (*dproj)[dd] = mproj; // projection operator for direct product
                            double (*dscal)[dd] = mscal;
                            
                            memset(dproj, 0, sizeof(double[dd][dd]));
                            
                            for(int s = 0;s < pg->sopsl;s++){
                                //if(pg->ct->irrep[dk].v[pg->sops[s].cla] == 0) continue;
                                permutationMatrix(&perm[i][s], mperm);
                                kron(d, mperm, ld, lst[s], dd, dscal);
                                
                                mlscale(pg->ct->irrep[dk].v[pg->sops[s].cla], dd, dscal, dscal);
                                mladd(dd, dscal, dproj, dproj);
                            }
                            
                            mlscale(((double) pg->ct->irrep[dk].d)/pg->order, dd, dproj, dproj);
                            
                            //printf("projection operator for symmetry species %s",pg->ct->irrep[dk].name);
                            //printTransform(dd, dd, dproj);
                            
                            mmtlmul(dd,dd,dproj,vspan*lvspan,mdec,dscal);
                            
                            //printf("mscal\n");
                            //printTransform(dd, vspan*lvspan, dscal);
                            
                            memset(dproj, 0, sizeof(double[dd][dd]));
                            mltranspose(dd, vspan*lvspan, dscal, dproj);
                
                            /////
                            //printf("mdec\n");
                            //printTransform(vspan*lvspan, dd, mdec);
                            //printf("dproj\n");
                            //printTransform(vspan*lvspan, dd, dproj);
                            //////
                            
                            dnirl = mgs(dd, dproj, mbasis, doirl, thresholds->orthogonalization/basisl);
                            if(dnirl - doirl != dvspan){
                                ret = MSYM_SUBSPACE_ERROR;
                                //printf("Ortogonal subspace decomposition of dimension (%d) inconsistent with span (%d) in %s",dnirl - doirl,dvspan,pg->ct->irrep[dk].name);
                                //printTransform(vspan*lvspan, dd, mbasis);
                                msymSetErrorDetails("Ortogonal subspace decomposition of dimension (%d) inconsistent with span (%d) in %s",dnirl - doirl,vspan,pg->ct->irrep[dk].name);
                                goto err;
                                
                            }
                            
                            
                            printf(" + %d%s", (int) round(mspan[dk]), pg->ct->irrep[dk].name);
                            
                        }
                        printf("\n");
                    }
                    /////
                    printf("decomposed %d%s x %d%s = ",(int) round(pspan[i][k]),pg->ct->irrep[k].name,(int) round(bspan[l][lk]),pg->ct->irrep[lk].name);
                    for(int prk = 0;prk < pg->ct->l;prk++) printf(" + %d%s", (int) round(mspan[prk]), pg->ct->irrep[prk].name);
                    printf("\n");
                    printTransform(vspan*lvspan, dd, mbasis);
                    //////
                    
                    for(int sk = 0,si=0;sk < pg->ct->l;sk++){
                        int svspan = pg->ct->irrep[sk].d*((int) round(mspan[sk])), dd = d*ld;
                        if(svspan == 0) continue;
                        double (*sbasis)[dd] = mbasis;
                        double (*sdec)[dd] = mdec;
                        if(pg->ct->irrep[sk].d > 1){
                        //if(round(mspan[sk]) > 1 && pg->ct->irrep[sk].d > 1){
                            printf("multi-dimensional partitioning required, running C3v test\n");
                            double (*dproj)[dd] = mproj; // projection operator for direct product
                            double (*dscal)[dd] = mscal;
                            
                            memset(dproj, 0, sizeof(double[dd][dd]));
                            memset(mdec, 0, sizeof(double[vspan*lvspan][dd]));
                            //memcpy(mdec, &sbasis[si], sizeof(double[vspan*lvspan][dd]));
                            //memset(&sbasis[si], 0, sizeof(double[vspan*lvspan][dd]));
                            //TODO: placeholder
                            
                            
                            
                            double c3v_char[2][6] =
                                {
                                    [0] = {1,0,0, 1,0,0},
                                    [1] = {1,0,0,-1,0,0}
                                };
                            
                            double d4h_char[2][16] =
                            {
                                
                                [0] = {1,0,0,0,0,0,0,0, 1,0,0,0,0,0,0,0},
                                [1] = {1,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0}
                            };
                            
                            char *cs_name[2] = {"A'","A''"};
                            
                            double (*cs_char)[pg->sopsl] = NULL;
                            
                            if(pg->type == POINT_GROUP_Cnv && pg->n == 3){
                                cs_char = c3v_char;
                                for(int s = 0;s < pg->sopsl;s++){
                                    printf("%lf, %lf", cs_char[0][s], cs_char[1][s]);
                                    printSymmetryOperation(&pg->sops[s]);
                                }
                            } else if(pg->type == POINT_GROUP_Dnh && pg->n == 4){
                                cs_char = d4h_char;
                                for(int s = 0;s < pg->sopsl;s++){
                                    printf("%lf, %lf", cs_char[0][s], cs_char[1][s]);
                                    printSymmetryOperation(&pg->sops[s]);
                                }
                            } else {
                                printf("%s need subgroup characters for operations\n", pg->name);
                                for(int s = 0;s < pg->sopsl;s++){
                                    printf("%d ", pg->sops[s].cla);
                                    printSymmetryOperation(&pg->sops[s]);
                                }
                                exit(1);
                                
                            }
                            
                            for(int dim = 0, doirl = 0, dnirl = 0;dim < pg->ct->irrep[sk].d;dim++, doirl = dnirl){
                                memset(dproj, 0, sizeof(double[dd][dd]));
                                for(int s = 0;s < pg->sopsl;s++){
                                    if(cs_char[dim][s] == 0) continue;
                                    permutationMatrix(&perm[i][s], mperm);
                                    kron(d, mperm, ld, lst[s], dd, dscal);
                                    mlscale(cs_char[dim][s], dd, dscal, dscal);
                                    mladd(dd, dscal, dproj, dproj);
                                }
                                mlscale(1.0/2.0, dd, dproj, dproj);
                                
                                printf("constructed Cs projection operator %s\n",cs_name[dim]);
                                printTransform(dd, dd, dproj);
                                
                                //mmtlmul(dd,dd,dproj,svspan,mdec,dscal);
                                mmtlmul(dd,dd,dproj,svspan,&sbasis[si],dscal);
                                
                                memset(dproj, 0, sizeof(double[dd][dd]));
                                mltranspose(dd, svspan, dscal, dproj);
                                dnirl = mgs(dd, dproj, mdec, doirl, thresholds->orthogonalization/basisl);
                                if(dnirl - doirl != round(mspan[sk])){
                                    ret = MSYM_SUBSPACE_ERROR;
                                    msymSetErrorDetails("Multi-dimensional subspace decomposition of dimension (%d) inconsistent with representaion (%d) in %s",dnirl - doirl,round(mspan[sk]),cs_name[dim]);
                                    goto err;
                                }
                                
                                for(int oi = doirl; oi < dnirl;oi++) vlnorm(dd, sdec[oi]);
                            }
                            //do apply symmetry here
                            
                            
                            memcpy(&sbasis[si], mdec, sizeof(double[svspan][dd]));
                            int mdim = round(mspan[sk]);
                            int found[3] = {0,0,0}, ifound = 0; //found[mdim]
                            //TODO: placeholder
                            double tmp1[2][3][60], mtmp[60][60], mtmp2[60][60];
                            //tmp1[mdim][mdim]
                            //mtmp2[dd][mdim]
                            
                            double (*g)[mdim][mdim] = tmp1;
                            double (*mt)[mdim] = mtmp;
                            
                            for(int s = 0;s < pg->sopsl && ifound < mdim;s++){
                                
                                permutationMatrix(&perm[i][s], mperm);
                                kron(d, mperm, ld, lst[s], dd, dscal);
                                for(int dim = 0;dim < mdim;dim++){
                                    //for(int bla -> irrep.d) more dimensions;
                                    if(found[dim]) continue;
                                    mmtlmul(dd, dd, dscal, mdim, &sdec[1*mdim], mt);
                                    mmlmul(1,dd,&sdec[dim], mdim, mt, &g[1][dim]);
                                    printf("v = sdec[%d]'*Sops[%d]*sdec[%d-%d] = ",dim,s,mdim,2*mdim-1);
                                    printTransform(1, mdim, &g[1][dim]);
                                    if(vlabs(mdim, g[1][dim]) > thresholds->zero) found[dim] = 1;
                                }
                            }
                            
                            printf("projected vectors before\n");
                            printTransform(svspan, dd, &sbasis[si]);
                            
                            //just 2 dimensions at the moment
                            for(int cdim = 0;cdim < mdim;cdim++){
                                vlnorm2(dd, sdec[cdim], sbasis[si+pg->ct->irrep[sk].d*cdim]);
                            }
                            for(int dim = 1; dim < pg->ct->irrep[sk].d && dim < 2;dim++){
                                mmlmul(mdim, mdim, g[dim], dd, &sdec[dim*mdim], dproj);
                                
                                for(int cdim = 0;cdim < mdim;cdim++){
                                    vlnorm2(dd, dproj[cdim], sbasis[si+dim+pg->ct->irrep[sk].d*cdim]);
                                }
                            }
                            
                            //for i=1:3 nv[:,2i-1] = S1[:,i];nv[:,2i] = Tmp[:,i]/norm(Tmp[:,i])
                            
                            printf("projected vectors after\n");
                            printTransform(svspan, dd, &sbasis[si]);
                            //exit(1);
                        }
                        for(int ir = 0;ir < svspan;ir += pg->ct->irrep[sk].d){
                            for(int n = l+1; n <= esinmax;n++){
                                if(esomap[0][n][l][0] == NULL) continue;
                                printf("adding %s subspace %d of %d for n = %d, l = %d with %d vectors\n",pg->ct->irrep[sk].name,ssi,ssl,n,l,svspan);
                                
                                msym_subspace_t *ss = &iss[ssi++];
                                
                                ss->type = ATOMIC_ORBITAL;
                                ss->irrep = sk;
                                ss->d = pg->ct->irrep[sk].d;
                                double (*space)[dd] = malloc(sizeof(double[ss->d][dd]));
                                //this will be incorrect since we have calculated that there should be partitioning of multidimensional spaces
                                for(int dim = 0; dim < ss->d;dim++){
                                    vlnorm2(dd, sbasis[si+dim+ir], space[dim]);
                                }
                                printTransform(ss->d, dd, space);
                                ss->space = (double*) space;
                                ss->basisl = 0;
                                ss->basis.o = malloc(sizeof(msym_orbital_t *[dd]));
                                for(int e = 0;e < es[i].length;e++){
                                    for(int m = -l;m <= l;m++){
                                        if(NULL == (ss->basis.o[ss->basisl++] = esomap[e][n][l][m+l])){
                                            ret = MSYM_SUBSPACE_ERROR;
                                            msymSetErrorDetails("Missing expected basis function for n = %d, l = %d, m = %d on atom %d when generating subspaces",n,l,m,e);
                                            goto err;
                                        }
                                    }
                                }
                                printSubspace(pg->ct, ss);
                            }
                        }
                        si += svspan;
                    }
                }
            }
        }
        
        //for(int i = 0; i < d;i++) vlnorm(d, porth[i]);
        
        printf("permutation es[%d] combine with spherical here (or perhaps above)\n",i);
        
        ///////
        for(int k = 0;k < pg->ct->l;k++) printf(" + %d%s", (int) round(pspan[i][k]), pg->ct->irrep[k].name);
        printf("\n");
        printTransform(d, d, porth);
        ///////
    }
    
    for(int k = 0;k < pg->ct->l;k++) ispan[k] = pg->ct->irrep[k].d*((int) round(dspan[k]));
    for(int i = 0;i < ssi;i++){printf("%d = %s\n",i,pg->ct->irrep[iss[i].irrep].name);};
    
    msym_subspace_t tss = {.subspacel = ssi, .subspace = iss, .d = 0, .basisl = 0, .space = NULL};
    filterSubspace(&tss);
    *subspace = tss.subspace;
    *subspacel = tss.subspacel;
    *ospan = ispan;
    
    printSubspace(pg->ct, &tss);
    
    
    
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
    
err:
    return ret;
}

msym_error_t generateOrbitalSubspaces(msym_point_group_t *pg, int esl, msym_equivalence_set_t *es, msym_permutation_t **perm, int basisl, msym_orbital_t basis[basisl], msym_thresholds_t *thresholds, int *subspacel, msym_subspace_t **subspace, int **pspan){
    msym_error_t ret = MSYM_SUCCESS;
    
    return testSpanConvert(pg,esl,es,perm,basisl,basis,thresholds,subspacel,subspace,pspan);
    
    //exit(1);
    
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
                    msymSetErrorDetails("Ortogonal ES subspace of dimension (%d) inconsistent with span (%d) in %s",nirrepl - lirrepl,ispan[k]*pg->ct->irrep[k].d,pg->ct->irrep[k].name);
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



