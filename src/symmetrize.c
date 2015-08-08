//
//  symmetrize.c
//  Symmetry
//
//  Created by Marcus Johansson on 04/02/15.
//  Copyright (c) 2015 Marcus Johansson. 
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "symmetrize.h"
#include "orbital.h"
#include "linalg.h"

#define SQR(x) ((x)*(x))

msym_error_t addProjectionOntoSubspace(int d, double orb[d], msym_subspace_t *ss, msym_orbital_t basis[d], double mem[d], double proj[d]);

msym_error_t symmetrizeMoleculeProject(msym_point_group_t *pg, int esl, msym_equivalence_set_t *es, msym_permutation_t **perm, msym_thresholds_t *thresholds, double *err);
msym_error_t symmetrizeMoleculeLinear(msym_point_group_t *pg, int esl, msym_equivalence_set_t *es, msym_permutation_t **perm, msym_thresholds_t *thresholds, double *err);

msym_error_t symmetrizeMolecule(msym_point_group_t *pg, int esl, msym_equivalence_set_t *es, msym_permutation_t **perm, msym_thresholds_t *thresholds, double *err){
    msym_error_t ret = MSYM_SUCCESS;
    if((pg->type == POINT_GROUP_Cnv || pg->type == POINT_GROUP_Dnh) && pg->n == 0){
        ret = symmetrizeMoleculeLinear(pg,esl,es,perm,thresholds,err);
    } else {
        ret = symmetrizeMoleculeProject(pg,esl,es,perm,thresholds,err);
    }
    
    return ret;
}
/* This is a projection into the fully symmetric space.
 * A little more computation than if we just recreate it from one atom,
 * but it is independant of the chosen atom and we can get the size
 * of the fully symmetric component.
 * The sizes of the individual equivalence sets are rather small anyways.
 */

msym_error_t symmetrizeMoleculeProject(msym_point_group_t *pg, int esl, msym_equivalence_set_t *es, msym_permutation_t **perm, msym_thresholds_t *thresholds, double *err){
    msym_error_t ret = MSYM_SUCCESS;
    double e = 0.0;
    double (*v)[3] = malloc(sizeof(double[pg->order][3]));
    for(int i = 0; i < esl;i++){
        if(es[i].length > pg->order){
            ret = MSYM_SYMMETRIZATION_ERROR;
            msymSetErrorDetails("Equivalence set (%d elements) larger than order of point group (%d)",es[i].length,pg->order);
            goto err;
        }
        memset(v, 0, sizeof(double[pg->order][3]));
        for(int j = 0; j < pg->sopsl;j++){
            for(int k = 0; k < es[i].length;k++){
                int p = perm[i][j].p[k];
                double sv[3];
                applySymmetryOperation(&pg->sops[j], es[i].elements[k]->v, sv);
                vadd(sv, v[p], v[p]);
            }
        }
        double sl = 0.0, ol = 0.0;
        for(int j = 0; j < es[i].length;j++){
            ol += vdot(es[i].elements[j]->v,es[i].elements[j]->v);
            sl += vdot(v[j],v[j]);
            vscale(1.0/((double)pg->order), v[j], es[i].elements[j]->v);
        }
        sl /= SQR((double)pg->order);
        if(!(es[i].length == 1 && ol <= thresholds->zero)) e += (ol-sl)/ol; //e = fmax(e,(ol-sl)/ol);
    }
    
    *err = sqrt(fmax(e,0.0)); //should never be < 0, but it's a dumb way to die
err:
    free(v);
    return ret;
}

msym_error_t symmetrizeMoleculeLinear(msym_point_group_t *pg, int esl, msym_equivalence_set_t *es, msym_permutation_t **perm, msym_thresholds_t *thresholds, double *err){
    msym_error_t ret = MSYM_SUCCESS;
    double e = 0.0;
    double (*v)[3] = malloc(sizeof(double[pg->order][3]));
    double (*vinf)[3] = malloc(sizeof(double[pg->order][3]));
    msym_symmetry_operation_t *cinf = NULL;
    
    for(int i = 0; i < pg->sopsl;i++){
        if(pg->sops[i].type == PROPER_ROTATION && pg->sops[i].order == 0) {
            cinf = &pg->sops[i];
            break;
        }
    }
    
    if(cinf == NULL){
        ret = MSYM_SYMMETRIZATION_ERROR;
        msymSetErrorDetails("Cannot find Cinf operation in linear point group");
        goto err;
    }
    
    for(int i = 0; i < esl;i++){
        if(es[i].length > pg->order){
            ret = MSYM_SYMMETRIZATION_ERROR;
            msymSetErrorDetails("Equivalence set (%d elements) larger than order of point group (%d)",es[i].length,pg->order);
            goto err;
        }
        
        memset(v, 0, sizeof(double[pg->order][3]));
        
        for(int k = 0; k < es[i].length;k++){
            vproj(es[i].elements[k]->v, cinf->v, vinf[k]);
        }
        
        for(int j = 0; j < pg->sopsl;j++){
            for(int k = 0; k < es[i].length;k++){
                int p = perm[i][j].p[k];
                double sv[3];
                applySymmetryOperation(&pg->sops[j], vinf[k], sv);
                vadd(sv, v[p], v[p]);
            }
        }
        double sl = 0.0, ol = 0.0;
        for(int j = 0; j < es[i].length;j++){
            ol += vdot(es[i].elements[j]->v,es[i].elements[j]->v);
            sl += vdot(v[j],v[j]);
            vscale(1.0/((double)pg->order), v[j], es[i].elements[j]->v);
        }
        sl /= SQR((double)pg->order);
        if(!(es[i].length == 1 && ol <= thresholds->zero)) e = fmax(e,(ol-sl)/ol);

        
    }
    
    *err = sqrt(e);
err:
    free(v);
    free(vinf);
    return ret;
}

msym_error_t symmetrizeOrbitals2(msym_point_group_t *pg, int ssl, msym_subspace_t *ss, int *span, int basisl, msym_orbital_t basis[basisl], msym_thresholds_t *thresholds, double orb[basisl][basisl],double symorb[basisl][basisl]){
    msym_error_t ret = MSYM_SUCCESS;
    double (*proj)[pg->ct->l][basisl] = malloc(sizeof(double[basisl][pg->ct->l][basisl]));
    double *mem = malloc(sizeof(double[basisl]));
    double (*comp)[pg->ct->l] = malloc(sizeof(double[basisl][pg->ct->l]));
    int *icomp = calloc(basisl,sizeof(int));
    int (*ispan) = calloc(pg->ct->l,sizeof(int));
    memset(proj,0,sizeof(double[basisl][pg->ct->l][basisl]));
    
    printf("SUBSPACES\n");
    msym_subspace_t tss = {.subspacel = ssl, .subspace = ss, .d = 0, .basisl = 0, .space = NULL};
    printSubspace(pg->ct, &tss);
    
    for(int o = 0;o < basisl;o++){
        double mcomp = -1.0;
        for(int k = 0;k < pg->ct->l;k++){
            for(int s = 0;s < ssl;s++){
                if(ss[s].irrep == k){
                    if(MSYM_SUCCESS != (ret = addProjectionOntoSubspace(basisl, orb[o], &ss[s], basis, mem, proj[o][k]))) goto err;
                }
            }
            comp[o][k] = vlabs(basisl, proj[o][k]);
            //printf("orbital %d compinent in %s = %lf\n",o,pg->ct->irrep[k].name,comp[o][k]);
            if(comp[o][k] > mcomp){
                icomp[o] = k;
                mcomp = comp[o][k];
            }
        }
        ispan[icomp[o]]++;
        printf("o = %d: ", o);
        printTransform(1,pg->ct->l,comp[o]);
    }
    
    
    
    for(int o = 0;o < basisl;o++){
        //ispan[icomp[o]]++;
        //printf("orbital %d (%lf) has largest component (%lf) in %s\n",o,vlabs(basisl,orb[o]),vlabs(basisl,proj[o][icomp[o]]),pg->ct->irrep[icomp[o]].name);
        //scale back to full length, this is a more reasonable option, but will look at that later
        //vlnorm2(basisl, proj[o][icomp[o]], symorb[o]);
        //vlscale(vlabs(basisl, orb[o]), basisl, symorb[o], symorb[o]);
        
        //just throw away
        vlcopy(basisl, proj[o][icomp[o]], symorb[o]);
    }
    
    //printf("Orbital span (vectors) = ");
    for(int k = 0;k < pg->ct->l;k++){
        if(ispan[k] != span[k]){
            msymSetErrorDetails("Projected orbitals do not span the expected irredicible representations. Expected %d%s, got %d",span[k],pg->ct->irrep[k].name,ispan[k]);
            ret = MSYM_SYMMETRIZATION_ERROR;
            goto err;
        }
        //printf(" + %d%s",ispan[k],pg->ct->irrep[k].name);
    }
    //printf("\n");
    
    
    free(ispan);
    free(icomp);
    free(comp);
    free(mem);
    free(proj);
    return ret;
err:
    free(ispan);
    free(icomp);
    free(comp);
    free(mem);
    free(proj);
    return ret;
}

/* TODO: lots of room for optimization in this code, pressed for time 
 * This code can no longer handle tree structured subpaces, they should be removed.
 * Way too complicated (and pointless) to do recursive multidimensional averaging
 */
msym_error_t symmetrizeOrbitals(msym_point_group_t *pg, int ssl, msym_subspace_t *ss, int *span, int basisl, msym_orbital_t basis[basisl], msym_thresholds_t *thresholds, double orb[basisl][basisl],double symorb[basisl][basisl]){
    msym_error_t ret = MSYM_SUCCESS;
    double (*proj)[pg->ct->l][basisl] = malloc(sizeof(double[basisl][pg->ct->l][basisl]));
    double (*dproj)[ssl] = calloc(basisl,sizeof(double[ssl]));
    
    double (*comp)[pg->ct->l] = calloc(basisl,sizeof(double[pg->ct->l]));
    int *icomp = calloc(basisl,sizeof(int));
    int (*ispan) = calloc(pg->ct->l,sizeof(int));
    memset(proj,0,sizeof(double[basisl][pg->ct->l][basisl]));
    memset(symorb,0,sizeof(double[basisl][basisl]));
    
    int md = 0;
    //could deduce from pg type but can't be bothered
    for(int k = 0;k < pg->ct->l;k++) md = (md > pg->ct->irrep[k].d ? md : pg->ct->irrep[k].d);
    
    double (*mem)[basisl] = malloc(sizeof(double[md*md+1][basisl]));
    int (*pf)[md] = calloc(basisl+1,sizeof(int[md]));
    double (*dmem)[md+1] = calloc(md,sizeof(double[md+1]));
    double *dmpf = (double *) dmem;

    
    /*printf("SUBSPACES\n");
    msym_subspace_t tss = {.subspacel = ssl, .subspace = ss, .d = 0, .basisl = 0, .space = NULL};
    printSubspace(pg->ct, &tss);*/
    
    /* not really needed anymore, we have can do this in the next loop */
    for(int o = 0;o < basisl;o++){
        double mcomp = -1.0;
        for(int k = 0;k < pg->ct->l;k++){
            for(int s = 0;s < ssl;s++){
                if(ss[s].irrep == k){
                    if(MSYM_SUCCESS != (ret = addProjectionOntoSubspace(basisl, orb[o], &ss[s], basis, mem[0], proj[o][k]))) goto err;
                }
            }
            comp[o][k] = vlabs(basisl, proj[o][k]);
            //printf("orbital %d compinent in %s = %lf\n",o,pg->ct->irrep[k].name,comp[o][k]);
            if(comp[o][k] > mcomp){
                icomp[o] = k;
                mcomp = comp[o][k];
            }
        }
        
        ispan[icomp[o]]++;
        printf("o = %d: ", o);
        printTransform(1,pg->ct->l,comp[o]);
    }
    
    for(int k = 0;k < pg->ct->l;k++){
        if(ispan[k] != span[k]){
            msymSetErrorDetails("Projected orbitals do not span the expected irredicible representations. Expected %d%s, got %d",span[k],pg->ct->irrep[k].name,ispan[k]);
            ret = MSYM_SYMMETRIZATION_ERROR;
            goto err;
        }
        //printf(" + %d%s",ispan[k],pg->ct->irrep[k].name);
    }
    
    
    for(int o = 0;o < basisl;o++){
        int ko = icomp[o];
        printf("basis %d components = ",o);
        for(int k = 0;k < pg->ct->l;k++){
            printf("%lf%s + ",comp[o][k],pg->ct->irrep[k].name);
        }
        printf("\n");
        
        if(pg->ct->irrep[ko].d > 1){
            for(int s = 0;s < ssl;s++){
                memset(mem[1], 0, sizeof(double[basisl]));
                if(MSYM_SUCCESS != (ret = addProjectionOntoSubspace(basisl, orb[o], &ss[s], basis, mem[0], mem[1]))) goto err;
                dproj[o][s] = vlabs(basisl, mem[1]);
            }
            
            printf("basis %s components = ",pg->ct->irrep[ko].name);
            for(int s = 0;s < ssl;s++){
                printf("%lf + ",dproj[o][s]);
            }
            printf("\n");
        }
    }
    
    for(int o = 0;o < basisl;o++){
        int ko = icomp[o], dim = pg->ct->irrep[ko].d, found = 0;
        memset(pf[o], -1, sizeof(int[md])); //2s complement
        
        
        for(int i = 0;i < o && !found;i++){
            for(int j = 0;j < md && !found;j++){
                found = pf[i][j] == o;
            }
        }
        
        if(found) continue;
        
        for(int i = 0;i < md;i++){dmpf[i] = DBL_MAX;}
        
        for(int po = 0; po < basisl && dim > 1;po++){
            if(icomp[po] != ko || o == po) continue;
            vlsub(ssl, dproj[o], dproj[po], mem[0]);
            double c = vlabs(ssl, mem[0]), mc = 0.0;
            int mic = 0;
            for(int i = 1;i < dim;i++){
                double diff = fabs(dmpf[i] - c);
                if(c < dmpf[i] && (diff > mc)){
                    mic = i;
                    mc = diff;
                }
            }
            dmpf[mic] = c;
            pf[o][mic] = po;
        }
        
        pf[o][0] = o;
        
        printf("basis %d (%s) partner functions = ",o,pg->ct->irrep[ko].name);
        for(int i = 0;i < dim;i++){
            printf("%d,",pf[o][i]);
        }
        printf("\n");
    }
    
    
    //should validate pf, only 1 orb of each
    
    for(int o = 0;o < basisl;o++){
        int dim = pg->ct->irrep[icomp[o]].d, md2 = md*md;
        if(pf[o][0] == -1 || dim <= 1) continue;
        printf("partner functions: ");
        for(int i = 0;i < dim;i++){
            printf("%d,",pf[o][i]);
        }
        printf(" require averaging and orthogonalization WARNING need to choose same component in all\n");
        for(int s = 0;s < ssl;s++){
            
            double avg = 0;
            msym_subspace_t *oss = &ss[s];
            if(oss->irrep != icomp[o]) {
                printf("skipping subspace %d (%s) only projecting into %s\n",s,pg->ct->irrep[oss->irrep].name,pg->ct->irrep[icomp[o]].name);
                continue;
            }
            memset(pf[basisl], -1, sizeof(int[md]));
            printf("components in subspace: %d: ",s);
            for(int d = 0;d < dim;d++){
                int pfo = pf[o][d];
                printf("%d = %lf, ",pfo,dproj[pfo][s]);
                avg += dproj[pfo][s];
                if(oss->d != dim){
                    printf("ERRORROO subspace dimension %d != %d\n",oss->d,dim);
                    printSubspace(pg->ct, oss);
                    exit(1);
                }
                
                for(int i = 0; i < dim;i++){
                    double (*space)[oss->basisl] = (double (*)[oss->basisl]) oss->space, c = -1.0;
                    memset(mem[md2], 0, sizeof(double[basisl]));
                    //printf("projecting onto basis %d of %d with basis length: %d\n",i,dim,oss->basisl);
                    for(int j = 0; j < oss->basisl;j++){
                        //printf("index = %ld -> %lf\n",oss->basis.o[j] - basis,space[i][j]);
                        mem[md2][oss->basis.o[j] - basis] = space[i][j];
                    }
                    vlproj(basisl, orb[pfo], mem[md2], mem[d*md+i]);
                    double dabs = vlabs(basisl, mem[d*md+i]);
                    dmem[d][i] = dabs;
                    if(dabs > c){
                        //choose max component so that no two have the same
                        int fi = 0;
                        for(fi = 0;fi < d;fi++){
                            if(pf[basisl][fi] == i) break;
                        }
                        if(fi == d){
                            c = dabs;
                            pf[basisl][d] = i;
                            //printf("orbital %d max component %lf = %d\n",d,c,i);
                        }
                    }
                    //printTransform(1,basisl,orb[pfo]);
                    //printTransform(1,basisl,mem[md]);
                    //printTransform(1,basisl,mem[d]);
                }
            }
            avg /= dim;
            
            
            printf(" average = %lf\n",avg);
            /* stupid error check;
            int pfm[2] = {1,1};
            for(int d = 0;d < dim;d++){
                pfm[0] *= d+1;
                pfm[1] *= pf[basisl][d]+1;
                printf("*%d",pf[basisl][d]+1);
            }
            
            if(pfm[0] != pfm[1]){
                printf("error %d != %d \n",pfm[0], pfm[1]);
                exit(1);
            }
             end stupid error check */
            int zspace = 0;
            for(int d = 0;d < dim;d++){
                if(dmem[d][pf[basisl][d]] <= thresholds->zero){
                    printf("found zero max component for partner function %d attempting switch\n",d);
                    int fim = -1, fid = pf[basisl][d];
                    double fimd = 0.0;
                    for(int fi = 0;fi < dim;fi++){
                        int fii = pf[basisl][fi];
                        if(dmem[fi][fid] > thresholds->zero && dmem[d][fii] > thresholds->zero && dmem[fi][fid] > fimd){
                            fim = fi;
                        }
                    }
                    if(fim >= 0){
                        
                        printf("found replacement orbital %d(%d) <-> %d(%d)\n",fim,pf[basisl][fim],d,fid);
                        for(int pr = 0;pr < dim;pr++){
                            printf("was max %d(%lf)\n",pf[basisl][pr],dmem[pr][pf[basisl][pr]]);
                        }
                        pf[basisl][d] = pf[basisl][fim];
                        pf[basisl][fim] = fid;
                    } else {
                        printf("could not find replacement orbital, removing from subspace\n");
                        for(int pr = 0;pr < dim;pr++){
                            printf("was max %d(%lf)\n",pf[basisl][pr],dmem[pr][pf[basisl][pr]]);
                        }
                        zspace = 1;
                        break;
                    }
                }
                printf("orbital %d max component in %d = %lf\n",d,pf[basisl][d],dmem[d][pf[basisl][d]]);
            }
            
            if(!zspace){
                printf("adding symmetrized orbitals %d\n",o);
                for(int d = 0;d < dim;d++){
                    
                    vlnorm2(basisl, mem[d*md+pf[basisl][d]], mem[md2]);
                    vlscale(avg, basisl, mem[md2], mem[md2]);
                    printf("adding degenerate component %d\n",d);
                    printTransform(1,basisl,mem[md2]);
                    vladd(basisl, symorb[pf[o][d]], mem[md2], symorb[pf[o][d]]);
                }
            }
            
        }
    }
    
    for(int o = 0;o < basisl;o++){
        int dim = pg->ct->irrep[icomp[o]].d;
        if(dim != 1) continue;
        vlcopy(basisl, proj[o][icomp[o]], symorb[o]);
    }
    
    //printf("Orbital span (vectors) = ");
    
    //printf("\n");
    
    
    free(ispan);
    free(icomp);
    free(comp);
    free(mem);
    free(proj);
    return ret;
err:
    free(ispan);
    free(icomp);
    free(comp);
    free(mem);
    free(proj);
    return ret;
}

msym_error_t symmetrizeTranslation(msym_point_group_t *pg, msym_equivalence_set_t *es, msym_permutation_t *perm, int pi, double translation[3]){
    msym_error_t ret = MSYM_SUCCESS;
    double (*v)[3] = calloc(es->length,sizeof(double[3]));
    
    for(int j = 0; j < pg->sopsl;j++){
        int p = perm[j].p[pi];
        double stranslation[3];
        applySymmetryOperation(&pg->sops[j], translation, stranslation);
        vadd(stranslation, v[p], v[p]);
    }
    
    double scale = ((double)es->length)/pg->order;
    
    for(int i = 0;i < es->length;i++){
        vscale(scale, v[i], v[i]);
        vadd(es->elements[i]->v,v[i],es->elements[i]->v);
    }
    
err:
    free(v);
    return ret;
}

msym_error_t addProjectionOntoSubspace(int d, double orb[d], msym_subspace_t *ss, msym_orbital_t basis[d], double mem[d], double proj[d]){
    msym_error_t ret = MSYM_SUCCESS;
    if(ss->subspacel){
        for(int i = 0;i < ss->subspacel;i++){
            if(MSYM_SUCCESS != (ret = addProjectionOntoSubspace(d, orb, &ss->subspace[i], basis, mem, proj))) goto err;
        }
    } else {
        for(int i = 0; i < ss->d;i++){
            double (*space)[ss->basisl] = (double (*)[ss->basisl]) ss->space;
            memset(mem, 0, sizeof(double[d]));
            for(int j = 0; j < ss->basisl;j++){
                mem[ss->basis.o[j] - basis] = space[i][j];
            }
            vlproj(d, orb, mem, mem);
            vladd(d, mem, proj, proj);
        }
    }

err:
    return ret;
}


