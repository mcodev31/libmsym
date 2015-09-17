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
void printNewSubspace(msym_character_table_t *ct, int l, msym_subspace_2_t ss[l]);
void tabPrintTransform(int r, int c, double M[r][c],int indent);
void tabprintf(char *format, int indent, ...);


//These are not the real orbitals but linear combinations:
//lm+ = p(z)((x+iy)^m+(x-iy)^m)
msym_error_t basisFunctionFromQuantumNumbers(int n, int l, int m, msym_basis_function_t *bf){
    
    if(l > n || abs(m) > l) goto err;
    
    bf->f.sh.n = n;
    bf->f.sh.l = l;
    bf->f.sh.m = m;
    
    memset(bf->name,0,sizeof(bf->name));
    
    switch(l) {
        case 0 :
            snprintf(bf->name, sizeof(bf->name), "%ds",n);
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
            snprintf(bf->name, sizeof(bf->name), "%dp%s",n,d);
            break;
        }
        case 2 : {
            //o->v = dpolynomial[m+l];
            char *d = (signbit(m) == 1 ? "-" : "+");
            snprintf(bf->name, sizeof(bf->name), "%dd%d%s",n,abs(m),d);
            break;
        }
        default : {
            char t = 'f' - 3 + l;
            char *d = (signbit(m) == 1 ? "-" : "+");
            snprintf(bf->name, sizeof(bf->name), "%d%c%d%s",n,t,abs(m),d);
        }
    }
    return MSYM_SUCCESS;
err:
    msymSetErrorDetails("Invalid orbital quantum numbers n:%d l:%d m:%d",n,l,m);
    return MSYM_INVALID_ORBITALS;
}

msym_error_t basisFunctionFromName(char *name, msym_basis_function_t *bf){
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
    
    return basisFunctionFromQuantumNumbers(n,l,m,bf);
    
err:
    msymSetErrorDetails("Invalid orbital name %s",name);
    return MSYM_INVALID_ORBITALS;
    
}

msym_error_t sphericalHarmoicPolynomial(int l, int m, double *poly){
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

msym_error_t generateBasisFunctionTransforms(int sopsl, msym_symmetry_operation_t sops[sopsl], int l, double transform[sopsl][2*l+1][2*l+1]){
    msym_error_t ret = MSYM_SUCCESS;
    int kdim = ipow(3,l), norbs = 2*l+1;
    double (*mkron)[kdim][kdim] = malloc(sizeof(double[2][kdim][kdim]));
    double (*poly)[kdim] = malloc(sizeof(double[norbs][kdim]));
    
    for(int m = -l; m <= l;m++){
        if(MSYM_SUCCESS != (ret = sphericalHarmoicPolynomial(l,m,poly[m+l]))) goto err;
        
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
        
        
        /* julia output of symmetry operations
        char buf[16];
        symmetryOperationShortName(&sops[i], 16, buf);
        printf("%s%d=\n",buf,i);
        
        printTransform(norbs, norbs, transform[i]);
        
        printf("Sops = union(Sops,(%s%d,))\n",buf,i);
        */
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

/* move this */
msym_error_t findDecentSubgroup(msym_point_group_t *pg, int irrep, int sgl, msym_subgroup_t sg[sgl], msym_thresholds_t *thresholds, msym_subgroup_t **osg){
    msym_error_t ret = MSYM_SUCCESS;
    *osg = NULL;
    msym_character_table_t *ct = pg->ct;
    switch(pg->type){
        case MSYM_POINT_GROUP_TYPE_Cn :
        {
            ret = MSYM_POINT_GROUP_ERROR;
            printf("Separate by imaginary component, what we have is the real components\n");
            msymSetErrorDetails("Point %s group has complex characters in symmetry species %s",pg->name, ct->s[irrep].name);
            goto err;
        }
        case MSYM_POINT_GROUP_TYPE_Cnh :
        {
            printf("Separate by imaginary component, what we have is the real components\n");
            ret = MSYM_POINT_GROUP_ERROR;
            msymSetErrorDetails("Lowering of symmetry would require complex characters for point group %s symmetry species %s",pg->name, ct->s[irrep].name);
            goto err;
        }
        case MSYM_POINT_GROUP_TYPE_Ci :
        case MSYM_POINT_GROUP_TYPE_Cs :
        {
            ret = MSYM_POINT_GROUP_ERROR;
            msymSetErrorDetails("Cannot lower symmetry of point group %s symmetry species %s",pg->name, ct->s[irrep].name);
            goto err;
        }
        case MSYM_POINT_GROUP_TYPE_I :
        case MSYM_POINT_GROUP_TYPE_Ih :
        {
            ret = MSYM_POINT_GROUP_ERROR;
            msymSetErrorDetails("Lowering of symmetry not implemented for point group %s symmetry species %s",pg->name, ct->s[irrep].name);
            goto err;
        }
        case MSYM_POINT_GROUP_TYPE_Dnd :
        case MSYM_POINT_GROUP_TYPE_Cnv : {
            if(ct->s[irrep].d != 2){
                ret = MSYM_POINT_GROUP_ERROR;
                msymSetErrorDetails("Cannot lower symmetry of point group %s symmetry species %s",pg->name, ct->s[irrep].name);
                goto err;
            }
            for(int i = 0;i < sgl;i++){
                if(sg[i].type == MSYM_POINT_GROUP_TYPE_Cs){
                    *osg = &sg[i];
                    break;
                }
            }
            break;
        }
        case MSYM_POINT_GROUP_TYPE_Dn : {
            if(ct->s[irrep].d != 2){
                ret = MSYM_POINT_GROUP_ERROR;
                msymSetErrorDetails("Cannot lower symmetry of point group %s symmetry species %s",pg->name, ct->s[irrep].name);
                goto err;
            }
            for(int i = 0;i < sgl;i++){
                if(sg[i].type == MSYM_POINT_GROUP_TYPE_Cn && sg[i].n == 2){
                    int h = 0;
                    for(int j = 0;j < sg[i].order;j++){
                        msym_symmetry_operation_t *sop = sg[i].sops[j];
                        if(sop->type == PROPER_ROTATION && sop->order == 2 && sop->orientation == HORIZONTAL){
                            h = 1;
                            printf("skipping h C2 subgroup\n");
                            break;
                        }
                    }
                    if(!h){
                        *osg = &sg[i];
                        break;
                    }
                }
            }
            break;
        }
        case MSYM_POINT_GROUP_TYPE_Dnh : {
            if(ct->s[irrep].d != 2){
                ret = MSYM_POINT_GROUP_ERROR;
                msymSetErrorDetails("Cannot lower symmetry of point group %s symmetry species %s",pg->name, ct->s[irrep].name);
                goto err;
            }
            for(int i = 0;i < sgl;i++){
            
                if(sg[i].type == MSYM_POINT_GROUP_TYPE_Cs){
                    int h = 0;
                    for(int j = 0;j < sg[i].order;j++){
                        msym_symmetry_operation_t *sop = sg[i].sops[j];
                        if(sop->type == REFLECTION && sop->orientation == HORIZONTAL){
                            h = 1;
                            printf("skipping h Cs subgroup\n");
                            break;
                        }
                    }
                    if(!h){
                        *osg = &sg[i];
                        break;
                    }
                }
            }
            break;
        }
        
        case MSYM_POINT_GROUP_TYPE_Oh :
        case MSYM_POINT_GROUP_TYPE_O  :
            ret = MSYM_POINT_GROUP_ERROR;
            msymSetErrorDetails("Not implemented for octahedral symmetry of symmetry species %s, requires Dn orientation determination",pg->name, ct->s[irrep].name);
            goto err;
        case MSYM_POINT_GROUP_TYPE_Th :
        case MSYM_POINT_GROUP_TYPE_T : {
            if(ct->s[irrep].d == 2){
                ret = MSYM_POINT_GROUP_ERROR;
                msymSetErrorDetails("Lowering of symmetry would require complex characters for point group %s symmetry species %s",pg->name, ct->s[irrep].name);
                goto err;
            }
            //fall through
        }
        case MSYM_POINT_GROUP_TYPE_Td : {
            if(ct->s[irrep].d == 2){
                for(int i = 0;i < sgl;i++){
                    if(sg[i].type == MSYM_POINT_GROUP_TYPE_Cs){
                        *osg = &sg[i];
                        break;
                    }
                }
            } else if(ct->s[irrep].d == 3){
                for(int i = 0;i < sgl;i++){
                    if(sg[i].type == MSYM_POINT_GROUP_TYPE_Dn && sg[i].n == 2){
                        *osg = &sg[i];
                        break;
                    }
                }
            } else {
                ret = MSYM_POINT_GROUP_ERROR;
                msymSetErrorDetails("Cannot lower symmetry of point group %s symmetry species %s",pg->name, ct->s[irrep].name);
                goto err;
            }
            break;
        }
        default:
        {
            ret = MSYM_POINT_GROUP_ERROR;
            msymSetErrorDetails("Unknown pointgroup %s symmetry species %s",pg->name, ct->s[irrep].name);
            goto err;
        }
    }
    if(*osg == NULL){
        ret = MSYM_POINT_GROUP_ERROR;
        msymSetErrorDetails("Could not find subgroup for point group %s symmetry species %s",pg->name, ct->s[irrep].name);
    }
err:
    return ret;
}

msym_error_t getDecentSubgroupCharacters(msym_point_group_t *pg, msym_subgroup_t *sg, double (*c)[pg->order]){
    msym_error_t ret = MSYM_SUCCESS;
    int e = 0;
    if(sg->type == MSYM_POINT_GROUP_TYPE_Dn && sg->n == 2){
        int index = 0;
        double d2c[3][3] = {
            [0] = { 1, -1, -1},
            [1] = {-1, -1,  1},
            [2] = {-1,  1, -1}
        };
        memset(c, 0, sizeof(double[3][pg->order]));
        for(int s = 0;s < pg->order && !((index == 3) && e);s++){
            for(int i = 0;i < sg->order;i++){
                if(&pg->sops[s] != sg->sops[i]) continue;
                if(index == 3 && e) break;
                if(pg->sops[s].type == IDENTITY){
                    e = 1;
                    c[0][s] = c[1][s] = c[2][s] = 1;
                    if(index == 3) break;
                } else {
                    c[0][s] = d2c[0][index];
                    c[1][s] = d2c[1][index];
                    c[2][s] = d2c[2][index];
                    index++;
                    if(index == 3 && e) break;
                }
            }
        }
    
    } else if((sg->type == MSYM_POINT_GROUP_TYPE_Cs) || (sg->type == MSYM_POINT_GROUP_TYPE_Cn && sg->n == 2)){
        int faxis = 0;
        memset(c, 0, sizeof(double[2][pg->order]));
        for(int s = 0;s < pg->order && !(e && faxis);s++){
            for(int i = 0;i < sg->order;i++){
                if(&pg->sops[s] != sg->sops[i]) continue;
                if(pg->sops[s].type == IDENTITY){
                    e = 1;
                    c[0][s] = c[1][s] = 1;
                    if(faxis) break;
                } else {
                    faxis = 1;
                    c[0][s] = 1;
                    c[1][s] = -1;
                    if(e) break;
                }
            }
        }

    } else {
        ret = MSYM_POINT_GROUP_ERROR;
        msymSetErrorDetails("Cannot determine symmetry decent character of subgroup %s",sg->name);
    }

err:
    return ret;
}



msym_error_t generateSALCSubspaces(msym_point_group_t *pg, int sgl, msym_subgroup_t sg[sgl], int esl, msym_equivalence_set_t *es, msym_permutation_t **perm, int basisl, msym_basis_function_t basis[basisl], msym_element_t *elements, msym_equivalence_set_t **esmap, msym_thresholds_t *thresholds, int *ossl, msym_subspace_2_t **oss, int **ospan){
    msym_error_t ret = MSYM_SUCCESS;
    msym_character_table_t *ct = pg->ct;
    int lmax = -1, nmax = 0;
    enum _msym_basis_type_t ftype = basis[0].type;
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
    } else if (ftype != MSYM_BASIS_TYPE_REAL_SPHERICAL_HARMONIC){
        msymSetErrorDetails("Basis function type not supported");
        ret = MSYM_INVALID_ORBITALS;
        return ret;
    }
    
    int projm = (2*lmax+1)*pg->order;
    
    double (*bspan)[ct->d] = calloc(lmax+1, sizeof(double[ct->d]));   // span of individual basis functions
    double (*pspan)[ct->d] = calloc(esl, sizeof(double[ct->d]));      // span of permutation operators
    
    double (*lspan)[ct->d] = calloc(esl, sizeof(double[ct->d]));        // total span of basis function on each ES
    double *rspan = calloc(ct->d, sizeof(double));                          // total direct product symmetrized basis
    double *dspan = calloc(ct->d, sizeof(double));                          // decoposed total span of symmetrized basis (double)
    
    double *mspan = calloc(ct->d, sizeof(double));                          // span decomposition memory
    double (*mproj)[projm] = calloc(projm, sizeof(double[projm]));              // projection operator memory
    double (*mscal)[projm] = calloc(projm, sizeof(double[projm]));              // projection scaling memory
    double (*mperm)[pg->order] = calloc(pg->order, sizeof(double[pg->order]));  // permutation memory
    double (*morth)[pg->order] = calloc(pg->order, sizeof(double[pg->order]));  // permutation orthoginalization memory
    double (*mbasis)[projm] = calloc(basisl, sizeof(double[projm]));            // basis function coefficients
    double (*mdec)[projm] = calloc(basisl, sizeof(double[projm]));              // directo product decomposition memory
    double (*sgc)[pg->order] = calloc(5,sizeof(double[pg->order]));
    
    double *mdcomp = NULL;
    double *mdproj = NULL;
    int *mdfound = NULL;
    
    int *ispan = calloc(ct->d, sizeof(int));                               // decoposed total span of symmetrized basis (int)
    int *isalc = calloc(ct->d, sizeof(int));                               // number of added salcs to irrep
    int *esnmax = calloc(esl, sizeof(int));                                     // max n in eqset
    
    msym_basis_function_t *(*esbfmap)[pg->order][nmax+1][lmax+1][2*lmax+1] = calloc(esl,sizeof(msym_basis_function_t *[pg->order][nmax+1][lmax+1][2*lmax+1]));
    
    struct _ltransforms {int d; void *t;} *lts = calloc(lmax+1,sizeof(struct _ltransforms)); // transformation matrices for basis functions
    
    int (*les)[lmax+1] = calloc(esl, sizeof(int[lmax+1]));                      // number of l-type basis functions in each ES
    
    msym_basis_function_t dbf = {.type = ftype};
    double (*ctable)[ct->d] = ct->table;
    
    msym_subspace_2_t *ss = calloc(ct->d, sizeof(msym_subspace_2_t));
    
    /* determine number of l-type basis functions in each ES */
    for(int o = 0;o < basisl;o++){
        les[esmap[basis[o].element - elements] - es][basis[o].f.sh.l] += basis[o].f.sh.m == 0;
    }
    
    for(int o = 0;o < basisl;o++){
        msym_basis_function_t *bf = &basis[o];
        int ei = (int)(bf->element - elements), esi;
        msym_equivalence_set_t *e = esmap[ei];
        for(esi = 0;esi < e->length && e->elements[esi] != bf->element;esi++){}; //could improve perf with a map here
        if(esi >= e->length){
            ret = MSYM_INVALID_ORBITALS;
            msymSetErrorDetails("Basis function does not map to any equivalence set");
            goto err;
        }
        if(bf->f.sh.n > esnmax[e - es]) esnmax[e - es] = bf->f.sh.n;
        esbfmap[e - es][esi][bf->f.sh.n][bf->f.sh.l][bf->f.sh.m+bf->f.sh.l] = bf;
    }

    /* calculate span of irreducible representations for basis functions and permutations */
    for(int s = 0; s < pg->order;s++){
        for(int l = 0; l <= lmax;l++){
            dbf.f.sh.l = l; //TODO: function type based
            double c = symmetryOperationCharacter(&pg->sops[s], &dbf);
            for(int k = 0;k < ct->d;k++){
                bspan[l][k] += c*ctable[k][pg->sops[s].cla];
            }
        }
        for(int i = 0;i < esl;i++){
            int uma = 0;
            for(int j = 0; j < perm[i][s].c_length;j++) uma += perm[i][s].c[j].l == 1; //this is why we loop over irreps twice
            for(int k = 0;k < ct->d;k++){
                pspan[i][k] += uma*ctable[k][pg->sops[s].cla];
            }
        }
    }
    
    
    
    /* scale basis span, calculate basis function transforms and symmetry species vectors */
    for(int l = 0;l <= lmax;l++){
        int d = 2*l+1;
        double (*lproj)[d] = mproj;
        double (*lscal)[d] = mscal;
        
        vlscale(1.0/pg->order, ct->d, bspan[l], bspan[l]);
        
        lts[l].d = d;
        lts[l].t = malloc(sizeof(double[pg->order+1][d][d]));
        
        double (*st)[d][d] = lts[l].t;
        
        if(MSYM_SUCCESS != (ret = generateBasisFunctionTransforms(pg->order, pg->sops, l, lts[l].t))) goto err; //TODO: generalize basis function concept
        
        memset(st[pg->order], 0, sizeof(double[d][d]));
        
        for(int k = 0, oirl = 0, nirl = 0;k < ct->d;k++, oirl = nirl){
            printf("bspan[%d][%d] = %lf this is only correct if the representation is irreducible\n",l,k,bspan[l][k]);
            int vspan = ct->s[k].d*((int) round(bspan[l][k]));
            if(vspan == 0) continue;
            
            memset(lproj, 0, sizeof(double[d][d]));
            for(int s = 0;s < pg->order;s++){
                mlscale(ctable[k][pg->sops[s].cla], d, st[s], lscal);
                mladd(d, lscal, lproj, lproj);
            }
            
            mlscale(((double) ct->s[k].d)/pg->order, d, lproj, lproj);
            nirl = mgs(d, lproj, st[pg->order], oirl, thresholds->orthogonalization/basisl);
            
            if(nirl - oirl != vspan){
                
                printTransform(d, d, st[pg->order]);
                ret = MSYM_SUBSPACE_ERROR;
                msymSetErrorDetails("Ortogonal subspace of dimension (%d) inconsistent with span (%d) in %s",nirl - oirl,vspan,ct->s[k].name);
                goto err;
                
            }
        }
        
        for(int i = 0; i < d;i++) vlnorm(d, st[pg->order][i]);
        
    }

    /* scale permutation span and calculate total basis function span on each ES */
    for(int i = 0;i < esl;i++){
        vlscale(1.0/pg->order, ct->d, pspan[i], pspan[i]);
        for(int l = 0; l <= lmax;l++){
            if(les[i][l] == 0) continue;
            les[i][l] /= es[i].length;
            vlscale((double) les[i][l], ct->d, bspan[l], mspan);
            vladd(ct->d, mspan, lspan[i], lspan[i]);
        }
    }

    
    /* calculate direct product of irreducible representations spanned by basis functions and permutations on each ES (don't really need to do this) */
    for(int i = 0;i < esl;i++){
        for(int k = 0;k < ct->d;k++){
            for(int j = 0;j < ct->d && round(pspan[i][k]) > 0;j++){
                if(round(lspan[i][j]) == 0) continue;
                directProduct(ct->d, ctable[k], ctable[j], mspan);
                vlscale(pspan[i][k]*lspan[i][j], ct->d, mspan, mspan);
                vladd(ct->d, mspan, rspan, rspan);
            }
        }
    }
    
    /* decompose direct product into irreducible representations */
    decomposeRepresentation(ct, rspan, dspan);
    int ddim_max = 0;
    for(int k = 0;k < ct->d;k++){
        ispan[k] = (int)round(dspan[k]);
        ddim_max = ddim_max > ispan[k] ? ddim_max : ispan[k];
        ss[k].s = k;
        ss[k].salcl = ispan[k];
        ss[k].salc = calloc(ss[k].salcl, sizeof(msym_salc_t));
    }
    
    //printf("ddim_max = %d\n",ddim_max);
    mdcomp = malloc(sizeof(double[5][ddim_max][ddim_max]));
    mdproj = malloc(sizeof(double[projm][projm]));
    mdfound = malloc(sizeof(int[5][ddim_max]));

    
    printf("decomposed\n");
    for(int prk = 0;prk < ct->d;prk++) printf(" + %d%s", ispan[prk], ct->s[prk].name);
    printf("\n");
    
    
    
    for(int i = 0; i < esl; i++){
        int d = es[i].length;
        double (*pproj)[d] = mproj, (*pscal)[d] = mscal, (*porth)[d] = morth;
        
        memset(porth, 0, sizeof(double[d][d]));
        for(int k = 0, oirl = 0, nirl = 0;k < ct->d;k++, oirl = nirl){
            int vspan = ct->s[k].d*((int) round(pspan[i][k]));
            if(vspan == 0) continue;
            
            memset(pproj, 0, sizeof(double[d][d]));
            for(int s = 0;s < pg->order;s++){
                if(ctable[k][pg->sops[s].cla] == 0) continue;
                permutationMatrix(&perm[i][s], mperm);
                mlscale(ctable[k][pg->sops[s].cla], d, mperm, pscal);
                mladd(d, pscal, pproj, pproj);
            }
            
            mlscale(((double) ct->s[k].d)/pg->order, d, pproj, pproj);
            nirl = mgs(d, pproj, porth, oirl, thresholds->orthogonalization/basisl);
            
            if(nirl - oirl != vspan){
                ret = MSYM_SUBSPACE_ERROR;
                msymSetErrorDetails("Ortogonal ES subspace of dimension (%d) inconsistent with span (%d) in %s",nirl - oirl,vspan,ct->s[k].name);
                goto err;
                
            }
            
            for(int oi = oirl; oi < nirl;oi++) vlnorm(d, porth[oi]);
            for(int l = 0;l <= lmax;l++){
                if(les[i][l] <= 0) continue;
                int li = 0, ld = lts[l].d;
                double (*lst)[ld][ld] = lts[l].t;
                for(int lk = 0;lk < ct->d;lk++){
                    int lvspan = ct->s[lk].d*((int) round(bspan[l][lk])), dd = d*ld;
                    if(lvspan == 0) continue;
                    kron2(vspan, d, &porth[oirl], lvspan, ld, &lst[pg->order][li], mbasis);
                    li += lvspan;
                    
                    directProduct(ct->d, ctable[k], ctable[lk], rspan);
                    vlscale(pspan[i][k]*bspan[l][lk], ct->d, rspan, rspan);
                    decomposeRepresentation(ct, rspan, mspan);
                    if(ct->s[k].d > 1 && ct->s[lk].d > 1){
                        memcpy(mdec, mbasis, sizeof(double[vspan*lvspan][dd]));
                        memset(mbasis, 0, sizeof(double[vspan*lvspan][dd]));
                        for(int dk = 0, doirl = 0, dnirl = 0;dk < ct->d;dk++, doirl = dnirl){
                            int dvspan = ct->s[dk].d*((int) round(mspan[dk]));
                            if(dvspan == 0) continue;
                            
                            double (*dproj)[dd] = mproj; // projection operator for direct product
                            double (*dscal)[dd] = mscal;
                            
                            memset(dproj, 0, sizeof(double[dd][dd]));
                            
                            for(int s = 0;s < pg->order;s++){
                                if(ctable[dk][pg->sops[s].cla] == 0) continue;
                                permutationMatrix(&perm[i][s], mperm);
                                kron(d, mperm, ld, lst[s], dd, dscal);
                                
                                mlscale(ctable[dk][pg->sops[s].cla], dd, dscal, dscal);
                                mladd(dd, dscal, dproj, dproj);
                            }
                            
                            mlscale(((double) ct->s[dk].d)/pg->order, dd, dproj, dproj);
                            mmtlmul(dd,dd,dproj,vspan*lvspan,mdec,dscal);
                            memset(dproj, 0, sizeof(double[dd][dd]));
                            mltranspose(dd, vspan*lvspan, dscal, dproj);
                            
                            dnirl = mgs(dd, dproj, mbasis, doirl, thresholds->orthogonalization/basisl);
                            if(dnirl - doirl != dvspan){
                                ret = MSYM_SUBSPACE_ERROR;
                                msymSetErrorDetails("Ortogonal subspace decomposition of dimension (%d) inconsistent with span (%d) in %s",dnirl - doirl,vspan,ct->s[dk].name);
                                goto err;
                                
                            }
                        }
                    }
                    for(int sk = 0,si=0;sk < ct->d;sk++){
                        int svspan = ct->s[sk].d*((int) round(mspan[sk])), dd = d*ld;
                        if(svspan == 0) continue;
                        double (*sbasis)[dd] = mbasis;
                        double (*sdec)[dd] = mdec;
                        if(ct->s[sk].d > 1){
                            printf("multi-dimensional partitioning required\n");
                            double (*dproj)[dd] = mproj; // projection operator for direct product
                            double (*dscal)[dd] = mscal;
                        
                            memset(dproj, 0, sizeof(double[dd][dd]));
                            memset(mdec, 0, sizeof(double[vspan*lvspan][dd]));

                            msym_subgroup_t *rsg = NULL;
                            if(MSYM_SUCCESS != (ret = findDecentSubgroup(pg, sk, sgl, sg, thresholds, &rsg))) goto err;
                            printf("using subgroup %s\n",rsg->name);
                            
                            if(MSYM_SUCCESS != (ret = getDecentSubgroupCharacters(pg, rsg, sgc))) goto err;
                            
                            //printTransform(5, pg->order, sgc);

                            
                            
                            for(int dim = 0, doirl = 0, dnirl = 0;dim < ct->s[sk].d;dim++, doirl = dnirl){
                                memset(dproj, 0, sizeof(double[dd][dd]));
                                for(int s = 0;s < pg->order;s++){
                                    if(sgc[dim][s] == 0) continue;
                                    permutationMatrix(&perm[i][s], mperm);
                                    kron(d, mperm, ld, lst[s], dd, dscal);
                                    mlscale(sgc[dim][s], dd, dscal, dscal);
                                    mladd(dd, dscal, dproj, dproj);
                                }

                                mlscale(1.0/2.0, dd, dproj, dproj);
                                
                                //printf("constructed subgroup projection operator %s\n",cs_name[dim]);
                                //printTransform(dd, dd, dproj);
                                
                                mmtlmul(dd,dd,dproj,svspan,&sbasis[si],dscal);
                                
                                memset(dproj, 0, sizeof(double[dd][dd]));
                                mltranspose(dd, svspan, dscal, dproj);
                                dnirl = mgs(dd, dproj, sdec, doirl, thresholds->orthogonalization/basisl);
                                if(dnirl - doirl != round(mspan[sk])){
                                    ret = MSYM_SUBSPACE_ERROR;
                                    msymSetErrorDetails("Multi-dimensional subspace decomposition of dimension (%d) inconsistent with representaion (%d) in subgroup irrep %d",dnirl - doirl,round(mspan[sk]),dim);
                                    goto err;
                                }
                                
                                for(int oi = doirl; oi < dnirl;oi++) vlnorm(dd, sdec[oi]);

                            }
                            //do apply symmetry here
                            memcpy(&sbasis[si], sdec, sizeof(double[svspan][dd]));
                            
                            int mdim = round(mspan[sk]);
                            int (*found)[mdim] = (int (*)[mdim]) mdfound;
                            double (*g)[mdim][mdim] = (double (*)[mdim][mdim]) mdcomp;
                            double (*mt)[mdim] = (double (*)[mdim]) mdproj;
                            
                            memset(found,0,sizeof(int[ct->s[sk].d][mdim]));
                            //printf("mdim = %d\n",mdim);
                            for(int s = 0;s < pg->order ;s++){
                                permutationMatrix(&perm[i][s], mperm);
                                kron(d, mperm, ld, lst[s], dd, dscal);
                                for(int dim = 0;dim < mdim;dim++){
                                    for(int sg = 1;sg < ct->s[sk].d;sg++){
                                        if(found[sg][dim]) continue;
                                        mmtlmul(dd, dd, dscal, mdim, &sdec[sg*mdim], mt);
                                        mmlmul(1,dd,&sdec[dim], mdim, mt, &g[sg][dim]);
                                        //printf("v = sdec[%d]'*Sops[%d]*sdec[%d-%d] = ",dim,s,mdim,2*mdim-1);
                                        //printTransform(1, mdim, &g[sg][dim]);
                                        if(vlabs(mdim, g[sg][dim]) > thresholds->zero) found[sg][dim] = 1;
                                    }
                                }
                            }
                            //printf("projected vectors before\n");
                            //printTransform(svspan, dd, &sbasis[si]);
                            
                            //just 2 dimensions at the moment
                            for(int cdim = 0;cdim < mdim;cdim++){
                                vlnorm2(dd, sdec[cdim], sbasis[si+ct->s[sk].d*cdim]);
                            }
                            for(int dim = 1; dim < ct->s[sk].d;dim++){
                                mmlmul(mdim, mdim, g[dim], dd, &sdec[dim*mdim], dproj);
                                
                                for(int cdim = 0;cdim < mdim;cdim++){
                                    vlnorm2(dd, dproj[cdim], sbasis[si+dim+ct->s[sk].d*cdim]);
                                }
                            }
                            
                            //for i=1:3 nv[:,2i-1] = S1[:,i];nv[:,2i] = Tmp[:,i]/norm(Tmp[:,i])
                            
                            //printf("projected vectors after\n");
                            //printTransform(svspan, dd, &sbasis[si]);
                            //exit(1);
                        }
                        for(int ir = 0;ir < svspan;ir += ct->s[sk].d){
                            for(int n = l+1; n <= esnmax[i];n++){
                                if(esbfmap[i][0][n][l][0] == NULL) continue;
                                if(isalc[sk] >= ispan[sk]){
                                    ret = MSYM_SUBSPACE_ERROR;
                                    msymSetErrorDetails("Exceeded calculated number of SALCs %d >= %d",isalc[sk],ispan[sk]);
                                    goto err;
                                }
                                //printf("adding salc %d of %d %d\n",isalc[sk],ispan[sk],sk);
                                msym_salc_t *salc = &ss[sk].salc[isalc[sk]];
                                salc->d = ct->s[sk].d;
                                double (*pf)[dd] = calloc(salc->d,sizeof(double[dd]));
                                //this will be incorrect since we have calculated that there should be partitioning of multidimensional spaces
                                for(int dim = 0; dim < salc->d;dim++){
                                    vlnorm2(dd, sbasis[si+dim+ir], pf[dim]);
                                }
                                //printTransform(ct->s[sk].d, dd, pf);
                                salc->pf = (double*) pf;
                                salc->fl = 0;
                                salc->f = calloc(dd,sizeof(msym_basis_function_t *));
                                for(int e = 0;e < es[i].length;e++){
                                    for(int m = -l;m <= l;m++){
                                        /*if(salc->fl >= dd){
                                            printf("added too many salcs %d\n",dd);
                                        }*/
                                        if(NULL == (salc->f[salc->fl++] = esbfmap[i][e][n][l][m+l])){
                                            ret = MSYM_SUBSPACE_ERROR;
                                            msymSetErrorDetails("Missing expected basis function for n = %d, l = %d, m = %d on atom %d when generating subspaces",n,l,m,e);
                                            goto err;
                                        }
                                    }
                                }
                                isalc[sk]++;
                                //printf("added %d functions\n",salc->fl);
                            }
                        }
                        si += svspan;
                    }
                }
            }
        }
        
        //for(int i = 0; i < d;i++) vlnorm(d, porth[i]);
        
        
        /*
        for(int k = 0;k < ct->d;k++) printf(" + %d%s", (int) round(pspan[i][k]), ct->s[k].name);
        printf("\n");
        printTransform(d, d, porth);
        */
    }
    
    //printNewSubspace(ct,ct->d,ss);
    *ospan = ispan;
    *ossl = ct->d;
    *oss = ss;

    
    free(bspan);
    free(pspan);
    free(lspan);
    free(rspan);
    free(dspan);
    free(mspan);
    free(mproj);
    free(mscal);
    free(mperm);
    free(morth);
    free(mbasis);
    free(mdec);
    free(sgc);
    free(mdcomp);
    free(mdproj);
    free(mdfound);
    //free(ispan); used
    free(isalc);
    free(esnmax);
    free(esbfmap);
    for(int l = 0;l <= lmax;l++){
        free(lts[l].t);
    }
    free(lts);
    free(les);
    //free(ss); used
    
    return ret;

err:
    free(bspan);
    free(pspan);
    free(lspan);
    free(rspan);
    free(dspan);
    free(mspan);
    free(mproj);
    free(mscal);
    free(mperm);
    free(morth);
    free(mbasis);
    free(mdec);
    free(sgc);
    free(mdcomp);
    free(mdproj);
    free(mdfound);
    free(ispan);
    free(isalc);
    free(esnmax);
    free(esbfmap);
    for(int l = 0;l <= lmax;l++){
        free(lts[l].t);
    }
    free(lts);
    free(les);
    for(int k = 0;k < ct->d;k++){
        for(int i = 0;i < ss[k].salcl;i++){
            free(ss[k].salc[i].f);
            free(ss[k].salc[i].pf);
        }
        free(ss[k].salc);
    }
    free(ss);
    return ret;
}

/*
msym_error_t generateOrbitalSubspaces(msym_point_group_t *pg, int esl, msym_equivalence_set_t *es, msym_permutation_t **perm, int basisl, msym_orbital_t basis[basisl], msym_thresholds_t *thresholds, int *subspacel, msym_subspace_t **subspace, int **pspan){
    msym_error_t ret = MSYM_SUCCESS;
    
    //return testSpan(pg,esl,es,perm,basisl,basis,thresholds,subspacel,subspace,pspan);
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
        lts[l].t = malloc(sizeof(double[pg->order][lts[l].d][lts[l].d]));
        if(MSYM_SUCCESS != (ret = generateOrbitalTransforms(pg->order, pg->sops, l, lts[l].t))) goto err;
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
            
            for(int s = 0;s < pg->order;s++){
                double (*lt)[lts[l].d][lts[l].d] = lts[l].t;
                permutationMatrix(&perm[i][s], mperm);
                kron(perm[i][s].p_length,mperm,lts[l].d,lt[s],d,mkron);
                
                
                char buf[16];
                symmetryOperationShortName(&pg->sops[s], 16, buf);
                printf("%s%d=\n",buf,s);
                
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
*/

/*
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
}*/

/*
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
}*/

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

void printNewSubspace(msym_character_table_t *ct, int l, msym_subspace_2_t ss[l]){
    for(int k = 0;k < l;k++){
        printf("Subspace %d %s\n",k,ct->s[ss[k].s].name);
        for(int i = 0;i < ss[k].salcl;i++){
            for(int j = 0;j < ss[k].salc[i].fl;j++){
                msym_basis_function_t *bf = ss[k].salc[i].f[j];
                if(bf == NULL){
                    printf("error bf\n");
                    exit(1);
                }
                printf("\t  %s%s\t\t",bf->element->name,bf->name);
            }
            printf("\n");

            double (*space)[ss[k].salc[i].fl] = (double (*)[ss[k].salc[i].fl]) ss[k].salc[i].pf;
            if(space == NULL){
                printf("error space\n");
                exit(1);
            }
            tabPrintTransform(ss[k].salc[i].d,ss[k].salc[i].fl,space,1);
        }
    }
}

