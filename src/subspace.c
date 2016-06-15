//
//  subspace.c
//  libmsym
//
//  Created by Marcus Johansson on 28/05/15.
//  Copyright (c) 2014 Marcus Johansson.
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

#include "msym.h"
#include "linalg.h"
#include "subspace.h"
#include "permutation.h"
#include "rsh.h"

#include "debug.h"

#define SQR(x) ((x)*(x))

#define PARTNER_THRESHOLD 1.0e-6

msym_error_t projectLinearlyIndependent(int dim, int vdim, double v[vdim][dim], int udim, double u[udim][dim], msym_thresholds_t *thresholds, double cmem[dim], double mem[dim][dim], double o[dim][dim], int *oirl);



void decomposeSubRepresentation(msym_point_group_t *pg, const msym_subgroup_t **rsg, double (*sgc)[5][pg->order], int span[pg->ct->d], int (*sgd)[5]){
    msym_character_table_t *ct = pg->ct;
    msym_symmetry_operation_t *sops = pg->sops;
    int sopsl = pg->order;
    int icosahedral = MSYM_POINT_GROUP_TYPE_I == pg->type || MSYM_POINT_GROUP_TYPE_Ih == pg->type;
    double (*ctable)[ct->d] = ct->table;
    memset(sgd, 0, ct->d*sizeof(*sgd));
    for(int ck = 0;ck < ct->d;ck++){
        for(int sk = 0;sk < ct->d;sk++){
            if(NULL == rsg[sk]) continue;
            int irrepd = ct->s[sk].d;
            if(!(icosahedral && irrepd == 5)){
                for(int d = 0; d < irrepd;d++){
                    double prod = 0;
                    for(int s = 0; s < sopsl;s++){
                        prod += sgc[sk][d][s]*ctable[ck][sops[s].cla];
                    }
                    //printf("decomposition of irrep %s into subdimension %d irrep %d = %lf\n",ct->s[ck].name, d, sk, prod/rsg[sk]->order);
                    sgd[sk][d] += (int) round(span[ck]*prod/rsg[sk]->order);
                }
            } else {
                int o = rsg[sk]->order;
                int order[] = {o,o,o,2,2};
                int dim[] = {1,2,2,1,1};
                for(int d = 0; d < irrepd;d++){
                    double prod = 0;
                    for(int s = 0; s < sopsl;s++){
                        prod += dim[d]*sgc[sk][d][s]*ctable[ck][sops[s].cla];
                    }
                    //printf("decomposition of irrep %s into subdimension %d irrep %d = %lf\n",ct->s[ck].name, d, sk, prod/order[d]);
                    sgd[sk][d] += (int) round(span[ck]*prod/order[d]);
                }
            }
        }
    }
}

msym_error_t generateBasisRepresentations(int n, int sopsl, msym_symmetry_operation_t sops[sopsl], int lmax, rsh_representations_t *lrsh){
    msym_error_t ret = MSYM_SUCCESS;
    for(int l = 0;l <= lmax;l++){
        int d = 2*l+1;
        lrsh[l].d = d;
        lrsh[l].t = malloc(sizeof(double[n][d][d]));
    }
    
    if(MSYM_SUCCESS != (ret = generateRSHRepresentations(sopsl, sops, lmax, lrsh))) goto err;
    
    return ret;
err:
    for(int l = 0;l <= lmax;l++){
        free(lrsh[l].t);
        lrsh[l].t = NULL;
        lrsh[l].d = 0;
    }
    return ret;
    
}


msym_error_t generateProjectionOperator(int d, int sopsl, double c[sopsl], msym_permutation_t perm[sopsl], int ld, double (*lsops)[ld][ld], double proj[perm->p_length*ld][perm->p_length*ld]){
    msym_error_t ret = MSYM_SUCCESS;
    
    int pd = perm->p_length;
    
    memset(proj,0,sizeof(double[pd*ld][pd*ld]));
    
    for(int s = 0;s < sopsl;s++){
        if(c[s] == 0) continue;
        for(int pi = 0, po = 0;pi < pd;pi++, po += ld){
            int pr = perm[s].p[pi]*ld;
            for(int li = 0;li < ld;li++){
                int r = pr + li;
                for(int lj = 0;lj < ld;lj++){
                    proj[r][po+lj] += lsops[s][li][lj]*c[s];
                }
            }
        }
    }
    
    mlscale(((double)d)/sopsl, pd*ld, proj, proj);
    
    return ret;
}


msym_error_t generatePermutationSubspaces(msym_point_group_t *pg, msym_permutation_t perm[pg->order], int span[pg->ct->d], msym_thresholds_t *thresholds, double pmem[4][perm->p_length][perm->p_length], double (**pss)[pg->ct->d], double ss[perm->p_length][perm->p_length]){
    msym_error_t ret = MSYM_SUCCESS;
    
    int dim = perm->p_length, sopsl = pg->order;
    msym_character_table_t *ct = pg->ct;
    msym_symmetry_operation_t *sops = pg->sops;
    double (*ctable)[ct->d] = ct->table;
    double (*proj)[dim] = pmem[0];
    
    memset(ss, 0, dim*sizeof(*ss));
    memset(pss, 0, ct->d*sizeof(*pss));
    
    for(int k = 0, oirl = 0, nirl = 0;k < ct->d;k++, oirl = nirl){
        int irrepd = ct->s[k].d, vspan = irrepd*span[k];
        
        if(vspan == 0) continue;
        
        memset(proj,0,dim*sizeof(*proj));
        for(int s = 0;s < sopsl;s++){
            double c = ctable[k][sops[s].cla];
            if(c == 0) continue;
            for(int i = 0;i < dim;i++){
                proj[perm[s].p[i]][i] += c;
            }
        }

        nirl = mgs2(dim, vspan,proj, ss, oirl, thresholds->orthogonalization);
        if(nirl - oirl != vspan){
            debug_printTransform(dim, dim, ss);
            ret = MSYM_SUBSPACE_ERROR;
            msymSetErrorDetails("Ortogonal permutation subspace of dimension (%d) inconsistent with span (%d) in %s",nirl - oirl,vspan,ct->s[k].name);
            goto err;
        }
        
        pss[k] = &ss[oirl];
        
    }
    
    //for(int i = 0; i < dim;i++) vlnorm(dim, ss[i]);

err:
    return ret;
}

msym_error_t generateSubspaces(msym_point_group_t *pg, msym_permutation_t perm[pg->order], int ld, double (*lrsops)[ld][ld], int span[pg->ct->d], double (*sgc)[5][pg->order], int (*sgd)[5], msym_thresholds_t *thresholds, double cmem[pg->order], double pmem[4][perm->p_length*ld][perm->p_length*ld], double (*(*pss)[5])[pg->ct->d], double ss[perm->p_length*ld][perm->p_length*ld]){
    msym_error_t ret = MSYM_SUCCESS;
    
    int pd = perm->p_length, dim = pd*ld;
    msym_character_table_t *ct = pg->ct;
    double (*ctable)[ct->d] = ct->table;
    double (*proj)[dim] = pmem[0], (*sspg)[dim] = pmem[1], (*sssg)[dim] = pmem[2], (*mem)[dim] = pmem[3];
    
    msym_symmetry_operation_t *sops = pg->sops;
    int sopsl = pg->order;
    int icosahedral = MSYM_POINT_GROUP_TYPE_I == pg->type || MSYM_POINT_GROUP_TYPE_Ih == pg->type;
    
    memset(ss, 0, dim*sizeof(*ss));
    memset(pss, 0, ct->d*sizeof(*pss));
    
    for(int k = 0, oirl = 0, nirl = 0;k < ct->d;k++, oirl = nirl){
        int irrepd = ct->s[k].d, pgvspan = irrepd*span[k];
        
        if(pgvspan == 0) continue;
        
        for(int s = 0;s < pg->order;s++) cmem[s] = ctable[k][sops[s].cla];
        
        if(MSYM_SUCCESS != (ret = generateProjectionOperator(irrepd,sopsl,cmem,perm,ld,lrsops,proj))) goto err;
        
        if(irrepd == 1){
            nirl = mgs2(dim, pgvspan, proj, ss, oirl, thresholds->orthogonalization);
            if(nirl - oirl != pgvspan){
                debug_printTransform(dim, dim, ss);
                ret = MSYM_SUBSPACE_ERROR;
                msymSetErrorDetails("Ortogonal subspace of dimension (%d) inconsistent with span (%d) in %s",nirl - oirl,pgvspan,ct->s[k].name);
                goto err;
            }
            pss[k][0] = &ss[oirl];
        } else if(!(icosahedral && irrepd == 5)){
            
            int pgnirl = mgs2(dim, pgvspan, proj, sspg, 0, thresholds->orthogonalization);
            for(int d = 0; d < irrepd;d++,oirl = nirl){
                
                if(MSYM_SUCCESS != (ret = generateProjectionOperator(1,sopsl,sgc[k][d],perm,ld,lrsops,proj))) goto err;
                
                int sgnirl = mgs2(dim, sgd[k][d], proj, sssg, 0, thresholds->orthogonalization);
                
                if(MSYM_SUCCESS != (ret = projectLinearlyIndependent(dim, pgnirl, sspg, sgnirl, sssg, thresholds, cmem, mem, ss, &nirl))) goto err;
                
                if(nirl - oirl != span[k]){
                    ret = MSYM_SUBSPACE_ERROR;
                    msymSetErrorDetails("Ortogonal subsubspace of dimension (%d) inconsistent with span (%d) in %s",nirl - oirl,span[k],ct->s[k].name);
                    goto err;
                }
                pss[k][d] = &ss[oirl];
                
            }
        } else {
            int idim[] = {1,2,2}, sdim[] = {3,4}, ssd = 0;
            int pgnirl = mgs2(dim, pgvspan, proj, sspg, 0, thresholds->orthogonalization);
            for(int d = 0; d < 3;d++,oirl = nirl){
                if(MSYM_SUCCESS != (ret = generateProjectionOperator(idim[d],sopsl,sgc[k][d],perm,ld,lrsops,proj))) goto err;
                
                int sgnirl = mgs2(dim, sgd[k][d], proj, sssg, 0, thresholds->orthogonalization);
                
                int id = idim[d];
                if(id > 1){
                    int n = 0;
                    if(MSYM_SUCCESS != (ret = projectLinearlyIndependent(dim, pgnirl, sspg, sgnirl, sssg, thresholds, cmem, mem, sssg, &n))) goto err;
                    sgnirl = n;
                    for(int sd = 0; sd < id;sd++,oirl = nirl){
                        int sid = sdim[sd];
                        if(MSYM_SUCCESS != (ret = generateProjectionOperator(1,sopsl,sgc[k][sid],perm,ld,lrsops,proj))) goto err;
                        
                        // sspg, sssg and proj are taken, use mem, and take proj as mem after
                        int ignirl = mgs2(dim, sgd[k][sid], proj, mem, 0, thresholds->orthogonalization);
                        
                        if(MSYM_SUCCESS != (ret = projectLinearlyIndependent(dim, sgnirl, sssg, ignirl, mem, thresholds, cmem, proj, ss, &nirl))) goto err;
                        
                        if(nirl - oirl != span[k]){
                            debug_printTransform(sgnirl, dim, sssg);
                            ret = MSYM_SUBSPACE_ERROR;
                            msymSetErrorDetails("Ortogonal icosahedral subsubspace of dimension (%d) inconsistent with span (%d) in %s",nirl - oirl,span[k],ct->s[k].name);
                            goto err;
                        }
                        pss[k][ssd] = &ss[oirl];
                        ssd++;
                    }
                } else {
                    if(MSYM_SUCCESS != (ret = projectLinearlyIndependent(dim, pgnirl, sspg, sgnirl, sssg, thresholds, cmem, mem, ss, &nirl))) goto err;
                    if(nirl - oirl != span[k]){
                        debug_printTransform(dim, dim, ss);
                        ret = MSYM_SUBSPACE_ERROR;
                        msymSetErrorDetails("Ortogonal icosahedral subspace of dimension (%d) inconsistent with span (%d) in %s",nirl - oirl,span[k],ct->s[k].name);
                        goto err;
                    }
                    pss[k][d] = &ss[oirl];
                    
                    ssd++;
                }
            }
        }
    }
    
    //for(int i = 0; i < dim;i++) vlnorm(dim, ss[i]);
    
err:
    
    return ret;
}

msym_error_t generateSubspacesMatrix(msym_point_group_t *pg, msym_permutation_t perm[pg->order], int ld, double (*lrsops)[ld][ld], int span[pg->ct->d], double (*sgc)[5][pg->order], msym_thresholds_t *thresholds, double cmem[pg->order], double pmem[4][perm->p_length*ld][perm->p_length*ld], double (*(*pss)[5])[pg->ct->d], double ss[perm->p_length*ld][perm->p_length*ld]){
    msym_error_t ret = MSYM_SUCCESS;
    
    int pd = perm->p_length, dim = pd*ld;
    msym_character_table_t *ct = pg->ct;
    double (*ctable)[ct->d] = ct->table;
    double (*projpg)[dim] = pmem[0], (*projsg)[dim] = pmem[1], (*projig)[dim] = pmem[2], (*mem)[dim] = pmem[3];
    double trace = 0.0;
    msym_symmetry_operation_t *sops = pg->sops;
    int sopsl = pg->order;
    int icosahedral = MSYM_POINT_GROUP_TYPE_I == pg->type || MSYM_POINT_GROUP_TYPE_Ih == pg->type;
    
    memset(ss, 0, dim*sizeof(*ss));
    memset(pss, 0, ct->d*sizeof(*pss));
    
    for(int k = 0, oirl = 0, nirl = 0;k < ct->d;k++, oirl = nirl){
        int irrepd = ct->s[k].d, pgvspan = irrepd*span[k];
        
        if(pgvspan == 0) continue;

        for(int s = 0;s < pg->order;s++) cmem[s] = ctable[k][sops[s].cla];
        
        if(MSYM_SUCCESS != (ret = generateProjectionOperator(irrepd,sopsl,cmem,perm,ld,lrsops,projpg))) goto err;
        
        if(irrepd == 1){
            nirl = mgs2(dim, pgvspan, projpg, ss, oirl, thresholds->orthogonalization);
            if(nirl - oirl != pgvspan){
                debug_printTransform(dim, dim, ss);
                ret = MSYM_SUBSPACE_ERROR;
                msymSetErrorDetails("Ortogonal subspace of dimension (%d) inconsistent with span (%d) in %s",nirl - oirl,pgvspan,ct->s[k].name);
                goto err;
            }
            pss[k][0] = &ss[oirl];
        } else if(!(icosahedral && irrepd == 5)){
            for(int d = 0; d < irrepd;d++,oirl = nirl){
                
                if(MSYM_SUCCESS != (ret = generateProjectionOperator(1,sopsl,sgc[k][d],perm,ld,lrsops,projsg))) goto err;
                clean_debug_printf("mmlmul %dx%d %d\n",dim,dim,__LINE__);
                mmlsymmul(dim, projsg, projpg, mem);
                //mmlmul(dim, dim, projsg, dim, projpg, mem);
                clean_debug_printf("done mmlmul %d\n",__LINE__);
                trace = mltrace(dim, mem);
                mlscale(span[k]/trace, dim, mem, mem);
                
                nirl = mgs2(dim, span[k], mem, ss, oirl, thresholds->orthogonalization);
                if(nirl - oirl != span[k]){
                    debug_printTransform(dim, dim, ss);
                    ret = MSYM_SUBSPACE_ERROR;
                    msymSetErrorDetails("Ortogonal subsubspace of dimension (%d) inconsistent with span (%d) in %s",nirl - oirl,span[k],ct->s[k].name);
                    goto err;
                }
                pss[k][d] = &ss[oirl];
                
            }
        } else {
            int idim[] = {1,2,2}, sdim[] = {3,4}, ssd = 0;
            for(int d = 0; d < 3;d++,oirl = nirl){
                if(MSYM_SUCCESS != (ret = generateProjectionOperator(idim[d],sopsl,sgc[k][d],perm,ld,lrsops,projsg))) goto err;
                clean_debug_printf("mmlmul %d\n",__LINE__);
                mmlsymmul(dim, projsg, projpg, projig);
                //mmlmul(dim, dim, projsg, dim, projpg, projig);
                clean_debug_printf("done mmlmul %d\n",__LINE__);
                int id = idim[d];
                if(id > 1){
                    for(int sd = 0; sd < id;sd++,oirl = nirl){
                        int sid = sdim[sd];
                        if(MSYM_SUCCESS != (ret = generateProjectionOperator(1,sopsl,sgc[k][sid],perm,ld,lrsops,projsg))) goto err;
                        clean_debug_printf("mmlmul %d\n",__LINE__);
                        mmlsymmul(dim, projsg, projig, mem);
                        //mmlmul(dim, dim, projsg, dim, projig, mem);
                        clean_debug_printf("done mmlmul %d\n",__LINE__);
                        trace = mltrace(dim, mem);
                        mlscale(span[k]/trace, dim, mem, mem); // We might have small components in these subspaces
                        
                        nirl = mgs2(dim, span[k], mem, ss, oirl, thresholds->orthogonalization);
                        if(nirl - oirl != span[k]){
                            debug_printTransform(dim, dim, mem);
                            ret = MSYM_SUBSPACE_ERROR;
                            msymSetErrorDetails("Ortogonal icosahedral subsubspace of dimension (%d) inconsistent with span (%d) in %s",nirl - oirl,span[k],ct->s[k].name);
                            goto err;
                        }
                        pss[k][ssd] = &ss[oirl];
                        ssd++;
                    }
                } else {
                    
                    nirl = mgs2(dim, span[k], projig, ss, oirl, thresholds->orthogonalization);
                    if(nirl - oirl != span[k]){
                        debug_printTransform(dim, dim, ss);
                        ret = MSYM_SUBSPACE_ERROR;
                        msymSetErrorDetails("Ortogonal icosahedral subspace of dimension (%d) inconsistent with span (%d) in %s",nirl - oirl,span[k],ct->s[k].name);
                        goto err;
                    }
                    pss[k][d] = &ss[oirl];
                    
                    ssd++;
                }
            }
        }
    }
    
    //for(int i = 0; i < dim;i++) vlnorm(dim, ss[i]);
    
err:
    
    return ret;
}

msym_error_t findSplittingFieldSubgroup(msym_point_group_t *pg, int irrep, int sgl, const msym_subgroup_t sg[sgl], msym_thresholds_t *thresholds, const msym_subgroup_t **osg){
    msym_error_t ret = MSYM_SUCCESS;
    *osg = NULL;
    msym_character_table_t *ct = pg->ct;
    
    switch(ct->s[irrep].d){
        case 2 : { //2-dimensional
            switch(pg->type){
                case MSYM_POINT_GROUP_TYPE_Cnv :
                case MSYM_POINT_GROUP_TYPE_Td  : {
                    for(int i = 0;i < sgl;i++){
                        if(sg[i].type == MSYM_POINT_GROUP_TYPE_Cs){
                            *osg = &sg[i];
                            break;
                        }
                    }
                    break;
                }
                case MSYM_POINT_GROUP_TYPE_Dn  :
                case MSYM_POINT_GROUP_TYPE_Dnd :
                case MSYM_POINT_GROUP_TYPE_Dnh :
                case MSYM_POINT_GROUP_TYPE_O   :
                case MSYM_POINT_GROUP_TYPE_Oh  : {
                    for(int i = 0;i < sgl;i++){
                        if(sg[i].type == MSYM_POINT_GROUP_TYPE_Cn && sg[i].n == 2){
                            int h = 0;
                            for(int j = 0;j < sg[i].order;j++){
                                msym_symmetry_operation_t *sop = sg[i].sops[j];
                                if(sop->type == PROPER_ROTATION && sop->order == 2 && sop->orientation == HORIZONTAL){
                                    h = 1;
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
                case MSYM_POINT_GROUP_TYPE_Cn  :
                case MSYM_POINT_GROUP_TYPE_Cnh :
                case MSYM_POINT_GROUP_TYPE_T   :
                case MSYM_POINT_GROUP_TYPE_Th  : {
                    ret = MSYM_SUBSPACE_ERROR;
                    msymSetErrorDetails("Cannot construct splitting field, point group %s has complex characters in symmetry species %s",pg->name, ct->s[irrep].name);
                    goto err;
                }
                default: break;
                    
            }
            break;
        }
        case 3 :{ //3-dimensional
            switch(pg->type){
                case MSYM_POINT_GROUP_TYPE_T  :
                case MSYM_POINT_GROUP_TYPE_Td :
                case MSYM_POINT_GROUP_TYPE_Th :
                case MSYM_POINT_GROUP_TYPE_O  :
                case MSYM_POINT_GROUP_TYPE_Oh :
                case MSYM_POINT_GROUP_TYPE_I  :
                case MSYM_POINT_GROUP_TYPE_Ih : {
                    for(int i = 0;i < sgl;i++){
                        if(sg[i].type == MSYM_POINT_GROUP_TYPE_Dn && sg[i].n == 2){
                            int h = 1;
                            for(int j = 0; j < sg[i].order; j++){
                                msym_symmetry_operation_t *sop = sg[i].sops[j];
                                if(!(sop->orientation == HORIZONTAL || sop->orientation == NONE)){
                                    h = 0;
                                    break;
                                }
                            }
                            if(h){
                                *osg = &sg[i];
                                break;
                            }
                        }
                    }
                }
                default: break;
            }
            break;
        }
        case 4 : { //4-dimensional
            if(!(MSYM_POINT_GROUP_TYPE_I == pg->type || MSYM_POINT_GROUP_TYPE_Ih == pg->type)){
                ret = MSYM_SUBSPACE_ERROR;
                msymSetErrorDetails("Invalid irrep dimension (4) when getting splitting field subgroup for non icosahedral point group %s", pg->name);
                goto err;
            }
            
            for(int i = 0;i < sgl;i++){
                if(sg[i].type == MSYM_POINT_GROUP_TYPE_Dn && sg[i].n == 2){
                    *osg = &sg[i];
                    break;
                }
            }
            
            break;
        }
        case 5 : { //5-dimensional
            if(!(MSYM_POINT_GROUP_TYPE_I == pg->type || MSYM_POINT_GROUP_TYPE_Ih == pg->type)){
                ret = MSYM_SUBSPACE_ERROR;
                msymSetErrorDetails("Invalid irrep dimension (5) when getting splitting field subgroup for non icosahedral point group %s", pg->name);
                goto err;
            }
            
            for(int i = 0;i < sgl;i++){
                if(sg[i].type == MSYM_POINT_GROUP_TYPE_Dn && sg[i].n == 5){
                    *osg = &sg[i];
                }
            }
            
            break;
        }
        default : break;
    }
    
    if(*osg == NULL){
        ret = MSYM_SUBSPACE_ERROR;
        msymSetErrorDetails("Could not find splitting filed subgroup in dimension %d for point group %s symmetry species %s",ct->s[irrep].d, pg->name, ct->s[irrep].name);
    }
err:
    return ret;
    
}

#define PHI 1.618033988749894848204586834

msym_error_t getSplittingFieldCharacters(msym_point_group_t *pg, const msym_subgroup_t *sg, double (*c)[pg->order], int *cd){
    msym_error_t ret = MSYM_SUCCESS;
    int e = 0;
    
    for(int i = 0;i < 5;i++){cd[i] = 1;}
    if((sg->type == MSYM_POINT_GROUP_TYPE_Cs) || (sg->type == MSYM_POINT_GROUP_TYPE_Cn && sg->n == 2)){
        int faxis = 0;
        memset(c, 0, sizeof(double[pg->order]));
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
        
    } else if(sg->type == MSYM_POINT_GROUP_TYPE_Dn && sg->n == 2){
        int index = 0;
        double d2c[4][3] = {
            [0] = { 1, -1, -1},
            [1] = {-1, -1,  1},
            [2] = {-1,  1, -1},
            [3] = { 1,  1,  1}
        };
        memset(c, 0, sizeof(double[3][pg->order]));
        for(int s = 0;s < pg->order && !((index == 3) && e);s++){
            for(int i = 0;i < sg->order;i++){
                if(&pg->sops[s] != sg->sops[i]) continue;
                if(index == 3 && e) break;
                if(pg->sops[s].type == IDENTITY){
                    e = 1;
                    c[0][s] = c[1][s] = c[2][s] = c[3][s] = 1;
                    if(index == 3) break;
                } else {
                    c[0][s] = d2c[0][index];
                    c[1][s] = d2c[1][index];
                    c[2][s] = d2c[2][index];
                    c[3][s] = d2c[3][index];
                    index++;
                    if(index == 3 && e) break;
                }
            }
        }
        
    } else if (sg->type == MSYM_POINT_GROUP_TYPE_Dn && sg->n == 5){
        int sfound = 0, found = 0;
        msym_subgroup_t *ssg = NULL;
        for(int i = 0;i < 2; i++){
            if(sg->generators[i]->type == MSYM_POINT_GROUP_TYPE_Cn && sg->generators[i]->order == 2){
                ssg = sg->generators[i];
                break;
            }
        }
        if(NULL == ssg){
            ret = MSYM_INVALID_CHARACTER_TABLE;
            msymSetErrorDetails("Cannot find C2 subgroup of D5 %s");
            goto err;
        }
        cd[1] = cd[2] = 2;
        double d5c[2][5]= {
            [0] = {1, 1/PHI, -PHI},
            [1] = {1,  -PHI, 1/PHI}
        };
        
        double d2c[5] = {1,0,0};
        
        double c2c[5] = {1,-1};
        
        
        for(int s = 0;s < pg->order && (found < 10 || !sfound);s++){
            for(int i = 0;i < ssg->order && !sfound;i++){
                if(&pg->sops[s] != ssg->sops[i]) continue;
                if(pg->sops[s].type == PROPER_ROTATION && pg->sops[s].order == 2){
                    c[3][s] = c2c[0];
                    c[4][s] = c2c[1];
                    break;
                }
            }
            
            for(int i = 0;i < sg->order;i++){
                if(&pg->sops[s] != sg->sops[i]) continue;
                if(pg->sops[s].type == IDENTITY){
                    c[0][s] = c[3][s] = c[4][s] = 1;
                    c[1][s] = c[2][s] = 2;
                    found++;
                } else if (pg->sops[s].type == PROPER_ROTATION && pg->sops[s].order == 5) {
                    int index = (pg->sops[s].power >> 1) & 1; //2 and 3 -> index 1
                    c[0][s] = d5c[index][0];
                    c[1][s] = d5c[index][1];
                    c[2][s] = d5c[index][2];
                    found++;
                } else if(pg->sops[s].type == PROPER_ROTATION && pg->sops[s].order == 2){
                    c[0][s] = d2c[0];
                    c[1][s] = d2c[1];
                    c[2][s] = d2c[2];
                    found++;
                }
            }
        }
    }
    
    else {
        ret = MSYM_INVALID_CHARACTER_TABLE;
        msymSetErrorDetails("Cannot determine splitting field characters of subgroup %s",sg->name);
        goto err;
    }
    
err:
    return ret;
}


msym_error_t projectLinearlyIndependent(int dim, int vdim, double v[vdim][dim], int udim, double u[udim][dim], msym_thresholds_t *thresholds, double cmem[dim], double mem[dim][dim], double o[dim][dim], int *oirl){
    msym_error_t ret = MSYM_SUCCESS;
    
    memset(mem, 0, dim*sizeof(*mem));
    for(int vd = 0; vd < vdim;vd++){
        for(int ud = 0; ud < udim; ud++){
            double c = vldot(dim, v[vd], u[ud]);
            vlscale(c, dim, u[ud], cmem);
            vladd(dim, cmem, mem[vd], mem[vd]);
        }
    }
    
    int mdim = vdim > udim ? udim : vdim;

    int nirl = mgs2(dim, mdim, mem, o, *oirl, thresholds->orthogonalization/sqrt(dim));
    
    for(int i = *oirl; i < nirl;i++) vlnorm(dim, o[i]);
    
    *oirl = nirl;
    
    return ret;
}

msym_error_t generateSplittingOperation(msym_point_group_t *pg, msym_permutation_t perm[pg->order], int ld, double (*lrsops)[ld][ld], int sgl, const msym_subgroup_t *sg, const msym_subgroup_t **rsg, double rsop[perm->p_length*ld][perm->p_length*ld], msym_symmetry_operation_t **osop){
    msym_error_t ret = MSYM_SUCCESS;
    
    int pd = perm->p_length, dim = pd*ld;
    msym_symmetry_operation_t *sop = NULL;
    
    switch(pg->type){
        case MSYM_POINT_GROUP_TYPE_Cs  :
        case MSYM_POINT_GROUP_TYPE_Ci  :
        case MSYM_POINT_GROUP_TYPE_T   :
        case MSYM_POINT_GROUP_TYPE_Td  :
        case MSYM_POINT_GROUP_TYPE_Th  :
        case MSYM_POINT_GROUP_TYPE_O   :
        case MSYM_POINT_GROUP_TYPE_Oh  :
            break;
        case MSYM_POINT_GROUP_TYPE_Cnv :
        case MSYM_POINT_GROUP_TYPE_Dn  :
        case MSYM_POINT_GROUP_TYPE_Dnh :
        case MSYM_POINT_GROUP_TYPE_Dnd : {
            if(pg->n > 2){
                sop = pg->primary;
                if(NULL == sop){
                    ret = MSYM_SUBSPACE_ERROR;
                    msymSetErrorDetails("Cannot determine splitting operation for point group %s", pg->name);
                    goto err;
                }
            }
            
            break;
        }
        case MSYM_POINT_GROUP_TYPE_I   :
        case MSYM_POINT_GROUP_TYPE_Ih  : {
            for(int i = 0;i < sgl;i++){
                int f = 0;
                for(f = 0; f < pg->ct->d && &sg[i] != rsg[f];f++);

                if(f == pg->ct->d && MSYM_POINT_GROUP_TYPE_Dn == sg[i].type && sg[i].n == 5){
                    msym_symmetry_operation_t **sgsops = sg[i].sops;
                    for(int j = 0; j < sg[i].order; j++){
                        if(PROPER_ROTATION == sgsops[j]->type && 5 == sgsops[j]->order && 1 == sgsops[j]->power){
                            sop = sgsops[j];
                            break;
                        }
                    }
                    break;
                }
            }
            
            if(NULL == sop){
                ret = MSYM_SUBSPACE_ERROR;
                msymSetErrorDetails("Cannot determine splitting operation for point group %s", pg->name);
                goto err;
            }
            
            break;
        }
            
        case MSYM_POINT_GROUP_TYPE_Cn  :
        case MSYM_POINT_GROUP_TYPE_Cnh :
            if(2 == pg->n) break;
            // fallthrough
        case MSYM_POINT_GROUP_TYPE_Sn  :
        case MSYM_POINT_GROUP_TYPE_K   :
        case MSYM_POINT_GROUP_TYPE_Kh  :
            ret = MSYM_SUBSPACE_ERROR;
            msymSetErrorDetails("Point group %s has no splitting operation", pg->name);
            goto err;
    }
    
    *osop = sop;
    
    if(NULL != sop){
        int s = (int) (sop - pg->sops);
        memset(rsop, 0, dim*sizeof(*rsop));
        for(int pi = 0, po = 0;pi < pd;pi++, po += ld){
            int pr = perm[s].p[pi]*ld;
            for(int li = 0;li < ld;li++){
                int r = pr + li;
                for(int lj = 0;lj < ld;lj++){
                    rsop[r][po+lj] += lrsops[s][li][lj];
                }
            }
        }
    }
    

err:
    return ret;

}


msym_error_t determinePartnerFunctionsSearch(msym_point_group_t *pg, msym_permutation_t perm[pg->order], int ld, double (*lrsops)[ld][ld], int dim, int sd, int sspan, double (*sdss)[dim], int sdvi[5], double (*split)[dim], double (*mem)[dim], int *li, double (*pf)[dim]){
    msym_error_t ret = MSYM_SUCCESS;
    
    int pd = perm->p_length;
    
    //need at least 3 dimensions
    double *f = mem[0];
    double *proj = mem[1];
    
    if(sd == 1){
        memcpy(pf, sdss, sspan*sizeof(*pf));
        return ret;
    }
    
    if(dim < 2){
        ret = MSYM_SUBSPACE_ERROR;
        msymSetErrorDetails("Unexpected dimension %d < 2 when determining partner functions", dim);
        goto err;
    }
    
    memset(pf, 0, dim*sizeof(*pf));
    
    for(int i = 0;i < sspan;i++){
        int found[5] = {1,0,0,0,0};
        int s = 0;
        for(s = 0; s < pg->order;s++){
            int df = 0;
            for(int d = 0; d < sd; d++) df += found[d];
            if(df == sd) break;

            if(IDENTITY == pg->sops[s].type) continue;
            // build symmetry operation
            memset(split, 0, dim*sizeof(*split));
            for(int pi = 0, po = 0;pi < pd;pi++, po += ld){
                int pr = perm[s].p[pi]*ld;
                for(int li = 0;li < ld;li++){
                    int r = pr + li;
                    for(int lj = 0;lj < ld;lj++){
                        split[r][po+lj] += lrsops[s][li][lj];
                    }
                }
            }
            
            memcpy(pf[i*sd], sdss[i], dim*sizeof(*sdss[i]));
            mvlmul(dim, dim, split, sdss[i], f);
            
            for(int d = 1; d < sd; d++){
                if(found[d]) continue;
                int pfi = i*sd+d;
                for(int j = 0; j < sspan;j++){
                    int pi = sdvi[d] + j;
                    vlproj(dim, f, sdss[pi], proj);
                    vladd(dim, proj, pf[pfi], pf[pfi]);
                }
                if(vlabs(dim, pf[pfi]) >= PARTNER_THRESHOLD){
                    vlnorm(dim, pf[pfi]);
                    found[d] = 1;
                }
            }
        }
        
        if(s == pg->order){
            ret = MSYM_SUBSPACE_ERROR;
            msymSetErrorDetails("Could not find partner functions in %d dimensional space", sd);
            goto err;
        }
    }
err:
    return ret;
}




msym_error_t determinePartnerFunctions(msym_point_group_t *pg, int r, msym_permutation_t perm[pg->order], int ld, double (*lrsops)[ld][ld], int dim, int sd, int sspan, double (*sdss)[dim], int sdvi[5], double (*split)[dim], msym_symmetry_operation_t *splitop, double (*mem)[dim], int *li, double (*pf)[dim]){
    msym_error_t ret = MSYM_SUCCESS;
    
    
    if(sd/r == 1){
        memcpy(pf, sdss, r*sspan*sizeof(*pf));
        return ret;
    }
    
    if(NULL == splitop) return determinePartnerFunctionsSearch(pg, perm, ld, lrsops, dim, sd, sspan, sdss, sdvi, split, mem, li, pf);
    
    //need at least 3 dimensions
    double *f = mem[0];
    double *proj = mem[1];
    
    

    if(dim < 2){
        ret = MSYM_SUBSPACE_ERROR;
        msymSetErrorDetails("Unexpected dimension %d < 2 when determining partner functions", dim);
        goto err;
    }
    
    memset(pf, 0, dim*sizeof(*pf));
    
    for(int i = 0;i < sspan;i++){
        memcpy(pf[i*sd], sdss[i], dim*sizeof(*sdss[i]));
        mvlmul(dim, dim, split, sdss[i], f);
        
        for(int d = 1; d < sd; d++){
            int pfi = i*sd+d;
            for(int j = 0; j < sspan;j++){
                int pi = sdvi[d] + j;
                vlproj(dim, f, sdss[pi], proj);
                vladd(dim, proj, pf[pfi], pf[pfi]);
            }
            if(vlabs(dim, pf[pfi]) < PARTNER_THRESHOLD){
                ret = MSYM_SUBSPACE_ERROR;
                msymSetErrorDetails("Could not determine partner functions in %d dimensional space", sd);
                goto err;
            }
            vlnorm(dim, pf[pfi]);
        }
    }
err:
    return ret;
}



msym_error_t symmetrySpeciesComponents(msym_point_group_t *pg, int srsl, msym_subrepresentation_space_t *srs, int basisl, msym_basis_function_t *basis, double *wf, double *s){
    msym_error_t ret = MSYM_SUCCESS;
    
    if(srsl != pg->ct->d){
        ret = MSYM_SUBSPACE_ERROR;
        msymSetErrorDetails("Unexpected subspace length (expected %d got %d)",pg->ct->d, srsl);
        goto err;
    }
    
    for(int k = 0;k < srsl;k++){
        double kcomp = 0.0;
        for(int s = 0;s < srs[k].salcl;s++){
            msym_salc_t *salc = &srs[k].salc[s];
            double (*space)[salc->fl] = (double (*)[salc->fl]) salc->pf;
            for(int d = 0;d < salc->d;d++){
                double c = 0.0;
                for(int j = 0; j < salc->fl;j++){
                    c += wf[salc->f[j] - basis]*space[d][j];
                }
                kcomp += SQR(c);
            }
        }
        s[k] = sqrt(kcomp);
    }
    
err:
    
    return ret;
}



msym_error_t generateSubrepresentationSpaces(msym_point_group_t *pg, int sgl, const msym_subgroup_t sg[sgl], int esl, msym_equivalence_set_t *es, msym_permutation_t **perm, int basisl, msym_basis_function_t basis[basisl], msym_element_t *elements, msym_equivalence_set_t **esmap, msym_thresholds_t *thresholds, int *osrsl, msym_subrepresentation_space_t **osrs, msym_basis_function_t ***osrsbf, int **ospan){
    msym_error_t ret = MSYM_SUCCESS;
    msym_character_table_t *ct = pg->ct;
    int lmax = -1, nmax = 0, eslmax = 0;
    enum _msym_basis_type ftype = basis[0].type;
    for(int i = 0;i < basisl;i++){
        if(basis[i].type != ftype) {nmax = -1; break;}
        lmax = basis[i].f.rsh.l > lmax ? basis[i].f.rsh.l : lmax;
        nmax = basis[i].f.rsh.n > nmax ? basis[i].f.rsh.n : nmax;
    }
    
    for(int i = 0;i < esl;i++){
        eslmax = es[i].length > eslmax ? es[i].length : eslmax;
    }
    
    if(lmax < 0 || nmax < 1){
        if(nmax == -1) msymSetErrorDetails("Basis functions are not of the same type");
        else msymSetErrorDetails("Invalid sperical harmonics quantum numbers");
        ret = MSYM_INVALID_BASIS_FUNCTIONS;
        return ret;
    } else if (ftype != MSYM_BASIS_TYPE_REAL_SPHERICAL_HARMONIC){
        msymSetErrorDetails("Basis function type not supported");
        ret = MSYM_INVALID_BASIS_FUNCTIONS;
        return ret;
    }
    
    int projm = (2*lmax+1)*eslmax;
    
    double (*pmem)[projm][projm] = calloc(7, sizeof(*pmem));        // Memory for calculating projection operatorsle
    
    if(NULL == pmem){
        ret = MSYM_MEMORY_ERROR;
        msymSetErrorDetails("Could not allocate %ld bytes of memory for SALC generation", 7*sizeof(*pmem));
        return ret;
    }
    
    double (*bspan)[ct->d] = calloc(lmax+1, sizeof(*bspan));        // span of individual basis functions
    double (*pspan)[ct->d] = calloc(esl, sizeof(*pspan));               // span of permutation operators
    
    double *rspan = calloc(ct->d, sizeof(*rspan));                  // total direct product symmetrized basis
    double *dspan = calloc(ct->d, sizeof(*dspan));                  // decoposed total span of symmetrized basis (double)
    
    double (*sdssmem)[projm] = pmem[1];
    double (*dssmem)[projm] = pmem[2];
    double (*pssmem)[pg->order] = pmem[3];
    double (*ssmem)[projm] = pmem[4];
    double (*pfmem)[projm] = pmem[5];
    double (*splitmem)[projm] = pmem[6];
    double *cmem = calloc(pg->order*(2*lmax+1), sizeof(*cmem));    // Don't change this to elsmax*(2*lmax+1) needed for sops
    double (*(*sspmem)[5])[pg->order*(lmax+1)] = calloc(ct->d, sizeof(*sspmem));
    
    double (**psspmem)[pg->order] = calloc(ct->d, sizeof(*psspmem));
    double (*(*lssp)[ct->d])[2*lmax+1] = calloc(lmax+1, sizeof(*lssp));
    
    double *mspan = calloc(ct->d, sizeof(double));                  // span decomposition memory
    double (*sgc)[5][pg->order] = calloc(ct->d,sizeof(*sgc));
    
    const msym_subgroup_t **rsg = calloc(ct->d, sizeof(*rsg));
    
    int (*sgd)[5] = calloc(ct->d,sizeof(*sgd));
    int *ispan = calloc(ct->d, sizeof(*ispan));                               // decoposed total span of symmetrized basis (int)
    int (*iespan)[lmax+1][ct->d] = calloc(esl, sizeof(*iespan));
    int (*ipspan)[ct->d] = calloc(esl, sizeof(*ipspan));               // span of permutation operators
    int (*ibspan)[ct->d] = calloc(lmax+1, sizeof(*ibspan));
    int *isalc = calloc(ct->d, sizeof(*isalc));                               // number of added salcs to irrep
    int *esnmax = calloc(esl, sizeof(*esnmax));                                     // max n in eqset
    
    msym_basis_function_t *(*esbfmap)[pg->order][nmax+1][lmax+1][2*lmax+1] = calloc(esl,sizeof(*esbfmap));
    
    
    msym_basis_function_t *(*srsbf) = calloc(basisl, sizeof(*srsbf));
    
    int (*srsbfmap)[nmax+1][lmax+1] = calloc(esl,sizeof(*srsbfmap));
    
    rsh_representations_t *lts = calloc(lmax+1,sizeof(*lts)); // transformation matrices for rsh basis functions
    
    int (*les)[lmax+1] = calloc(esl, sizeof(*les));                      // number of l-type basis functions in each ES
    
    msym_basis_function_t dbf = {.type = ftype};
    double (*ctable)[ct->d] = ct->table;
    
    msym_subrepresentation_space_t *srs = calloc(ct->d, sizeof(*srs));
    
    /* determine number of l-type basis functions in each ES */
    for(int o = 0;o < basisl;o++){
        les[esmap[basis[o].element - elements] - es][basis[o].f.rsh.l] += basis[o].f.rsh.m == 0;
    }
    
    if(MSYM_SUCCESS != (ret = generateBasisRepresentations(pg->order+1, pg->order, pg->sops, lmax, lts))) goto err;
    
    
    for(int k = 0; k < ct->d;k++){
        if(ct->s[k].d > 1){
            if(MSYM_SUCCESS != (ret = findSplittingFieldSubgroup(pg, k, sgl, sg, thresholds, &rsg[k]))) goto err;
            if(MSYM_SUCCESS != (ret = getSplittingFieldCharacters(pg, rsg[k], sgc[k], sgd[k]))) goto err; //TODO: remove sgd
        }
    }
    
    for(int o = 0;o < basisl;o++){
        msym_basis_function_t *bf = &basis[o];
        int ei = (int)(bf->element - elements), esi = 0;
        msym_equivalence_set_t *e = esmap[ei];
        for(esi = 0;esi < e->length && e->elements[esi] != bf->element;esi++){}; //could improve perf with a map here
        if(esi >= e->length){
            ret = MSYM_INVALID_BASIS_FUNCTIONS;
            msymSetErrorDetails("Basis function does not map to any equivalence set");
            goto err;
        }
        if(bf->f.rsh.n > esnmax[e - es]) esnmax[e - es] = bf->f.rsh.n;
        esbfmap[e - es][esi][bf->f.rsh.n][bf->f.rsh.l][bf->f.rsh.m+bf->f.rsh.l] = bf;
    }
    
    int srsbfi = 0;
    
    for(int i = 0;i < esl;i++){
        for(int n = 0;n <= nmax;n++){
            for(int l = 0;l <= lmax;l++){
                srsbfmap[i][n][l] = srsbfi;
                for(int e = 0;e < es[i].length;e++){
                    msym_basis_function_t *lbf = esbfmap[i][e][n][l][0];
                    for(int m = -l;m <= l;m++){
                        msym_basis_function_t *bf = esbfmap[i][e][n][l][m+l];
                        if((NULL == bf && NULL != lbf) || (NULL == lbf && NULL != bf)) {
                            lbf = NULL == lbf ? bf : lbf;
                            ret = MSYM_INVALID_BASIS_FUNCTIONS;
                            msymSetErrorDetails("Found basis function %s but function where m = %d is missing on element %d of equivalence set %d",lbf->name,m,e,i);
                            goto err;
                        }
                        if(NULL != bf){
                            lbf = bf;
                            srsbf[srsbfi] = bf;
                            srsbfi++;
                        }
                    }
                }
                //when this is moved to a separate function it can be used and availability check (replace esbfmap)
                //if(srsbfmap[i][n][l] == srsbfi) srsbfmap[i][n][l] = -1;
            }
        }
    }
    
    if(srsbfi != basisl){
        ret = MSYM_SUBSPACE_ERROR;
        msymSetErrorDetails("Unexpected number of basis functions in subrepresentation map (expected %d, got %d)",basisl, srsbfi);
        goto err;
    }
    
    /* calculate span of irreducible representations for basis functions and permutations */
    for(int s = 0; s < pg->order;s++){
        for(int l = 0; l <= lmax;l++){
            dbf.f.rsh.l = l; //TODO: function type based
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
    
    /* scale reducible representations */
    for(int k = 0;k < ct->d;k++){
        double r = ct->s[k].r;
        if(r > 1){
            for(int l = 0; l <= lmax;l++) bspan[l][k] /= r;
            for(int i = 0;i < esl;i++) pspan[i][k] /= r;
        }
    }
    
    /* scale basis span, calculate basis function transforms and symmetry species vectors */
    for(int l = 0;l <= lmax;l++){
        int d = lts[l].d;
        double (*lproj)[d] = pmem[0];
        double (*lscal)[d] = pmem[1];
        
        vlscale(1.0/pg->order, ct->d, bspan[l], bspan[l]);
        
        double (*st)[d][d] = lts[l].t;
        memset(st[pg->order], 0, sizeof(double[d][d]));
        
        for(int k = 0, oirl = 0, nirl = 0;k < ct->d;k++, oirl = nirl){
            
            int vspan = ct->s[k].d*((int) round(bspan[l][k]));
            if(vspan == 0) continue;
            
            ibspan[l][k] = (int) round(bspan[l][k]);
            
            memset(lproj, 0, sizeof(double[d][d]));
            for(int s = 0;s < pg->order;s++){
                double c = ctable[k][pg->sops[s].cla];
                if(c != 0){
                    mlscale(ctable[k][pg->sops[s].cla], d, st[s], lscal);
                    mladd(d, lscal, lproj, lproj);
                }
            }
            
            nirl = mgs2(d, vspan,lproj, st[pg->order], oirl, thresholds->orthogonalization);
            
            if(nirl - oirl != vspan){
                debug_printTransform(d, d, st[pg->order]);
                ret = MSYM_SUBSPACE_ERROR;
                msymSetErrorDetails("Ortogonal subspace of dimension (%d) inconsistent with span (%d) in %s",nirl - oirl,vspan,ct->s[k].name);
                goto err;
                
            }
            
            lssp[l][k] = &st[pg->order][oirl];
            
        }
        
        for(int i = 0; i < d;i++) vlnorm(d, st[pg->order][i]);
        
    }
    
    /* scale permutation span and calculate total basis function span on each ES */
    for(int i = 0;i < esl;i++){
        vlscale(1.0/pg->order, ct->d, pspan[i], pspan[i]);
        for(int l = 0; l <= lmax;l++){
            if(les[i][l] == 0) continue;
            les[i][l] /= es[i].length;
            memset(dspan, 0, sizeof(double[ct->d]));
            for(int k = 0;k < ct->d;k++){
                for(int j = 0;j < ct->d && round(pspan[i][k]) > 0;j++){
                    directProduct(ct->d, ctable[k], ctable[j], mspan);
                    vlscale(pspan[i][k]*bspan[l][j], ct->d, mspan, mspan);
                    vladd(ct->d, mspan, dspan, dspan);
                }
            }
            
            decomposeRepresentation(ct, dspan, mspan);
            for(int k = 0;k < ct->d;k++){
                iespan[i][l][k] = (int)round(mspan[k]/ct->s[k].r);
                ispan[k] += les[i][l]*iespan[i][l][k];
            }
        }
        for(int k = 0;k < ct->d;k++) ipspan[i][k] = (int)round(pspan[i][k]);
    }

    for(int k = 0;k < ct->d;k++){
        srs[k].s = k;
        srs[k].salcl = ispan[k];
        srs[k].salc = calloc(srs[k].salcl, sizeof(msym_salc_t));
    }
    
    clean_debug_printf("decomposed %d\n", ct->d);
    for(int prk = 0;prk < ct->d;prk++){
        if(prk < ct->d - 1) clean_debug_printf("%d%s + ", ispan[prk], ct->s[prk].name);
        else clean_debug_printf("%d%s\n", ispan[prk], ct->s[prk].name);
    }
    
    int dbasisl = 0;
    
    for(int k = 0;k < ct->d;k++){
        dbasisl += ct->s[k].d*ispan[k];
    }
    
    if(dbasisl != basisl){
        ret = MSYM_SUBSPACE_ERROR;
        msymSetErrorDetails("Unexpected number of basis functions in decomposition (expected %d, got %d)",basisl,dbasisl);
        goto err;
    }
    
    
    for(int i = 0;i < esl;i++){
        for(int l = 0;l <= lmax;l++){

            if(les[i][l] == 0) continue;
            
            
            int ld = lts[l].d, esd = es[i].length, dim = esd*ld, li = 0;
            double (*lrsops)[ld][ld] = lts[l].t;
            double (*(*ssp)[5])[dim] = sspmem;
            double (*ss)[dim] = ssmem;
            double (**pssp)[esd] = psspmem;
            double (*pss)[esd] = pssmem;
            double (*split)[dim] = splitmem;
            msym_symmetry_operation_t *splitop = NULL;
            
            clean_debug_printf("e decomposed %d\n", ct->d);
            for(int prk = 0;prk < ct->d;prk++){
                if(prk < ct->d - 1) clean_debug_printf("%d%s + ", iespan[i][l][prk], ct->s[prk].name);
                else clean_debug_printf("%d%s\n", iespan[i][l][prk], ct->s[prk].name);
            }
            
            decomposeSubRepresentation(pg, rsg, sgc, iespan[i][l], sgd);
            if(MSYM_SUCCESS != (ret = generateSplittingOperation(pg, perm[i], ld, lrsops, sgl, sg, rsg, split, &splitop))) goto err;
            if(MSYM_SUCCESS != (ret = generateSubspaces(pg, perm[i], ld, lrsops, iespan[i][l], sgc, sgd, thresholds, cmem, pmem, ssp, ss))) goto err;
            // Must be done after generateSubspaces since memory overlaps with pss for icosahedral groups
            if(MSYM_SUCCESS != (ret = generatePermutationSubspaces(pg, perm[i], ipspan[i], thresholds, pmem, pssp, pss))) goto err;
            
            for(int pk = 0;pk < ct->d;pk++){
                int pvspan = ipspan[i][pk]*ct->s[pk].d;
                if(pvspan == 0) continue;
                for(int lk = 0;lk < ct->d;lk++){
                    int lvspan = ibspan[l][lk]*ct->s[lk].d, vspan = pvspan*lvspan;
                    if(lvspan == 0) continue;
                    
                    double (*dss)[dim] = dssmem;
                    
                    kron2(pvspan, esd, pssp[pk], lvspan, ld, lssp[l][lk], dss);
                    
                    directProduct(ct->d, ctable[pk], ctable[lk], rspan);
                    vlscale(pspan[i][pk]*bspan[l][lk], ct->d, rspan, rspan);
                    decomposeRepresentation(ct, rspan, mspan);
        
                    for(int dk = 0; dk < ct->d; dk++){
                        int sspan =  (int) round(mspan[dk]/ct->s[dk].r);
                        if(sspan == 0) continue;
                        
                        double (*sdss)[dim] = sdssmem;
                        
                        int sd = ct->s[dk].d;
                        
                        int sdvi[5] = {-1,-1,-1,-1,-1};
                        
                        for(int d = 0, oirl = 0, nirl = 0; d < sd;d++, oirl = nirl){
                            if(MSYM_SUCCESS != (ret = projectLinearlyIndependent(dim, vspan, dss, iespan[i][l][dk], ssp[dk][d], thresholds, cmem, pmem[0], sdss, &nirl))) goto err;
                            
                            int sdvl = nirl - oirl;
                            sdvi[d] = oirl;
                            
                            if(sdvl != sspan){
                                debug_printTransform(vspan, dim, dss);
                                ret = MSYM_SUBSPACE_ERROR;
                                msymSetErrorDetails("Linear projection subspace of dimension (%d) inconsistent with span (%d) in %s component %d",sdvl,sspan,ct->s[dk].name,d);
                                goto err;
                            }
                            
                            
                        }
                        
                        double (*pf)[dim] = pfmem;
                        
                        if(MSYM_SUCCESS != (ret = determinePartnerFunctions(pg, ct->s[dk].r, perm[i], ld, lrsops, dim, sd, sspan, sdss, sdvi, split, splitop, pmem[0], &li, pf))) goto err;
                        
                        for(int si = 0, pfi = 0; si < sspan;si++){
                            for(int n = l+1; n <= esnmax[i];n++){
                                if(esbfmap[i][0][n][l][0] == NULL) continue;
                                if(isalc[dk] >= ispan[dk]){
                                    ret = MSYM_SUBSPACE_ERROR;
                                    msymSetErrorDetails("Exceeded calculated number of SALCs %d >= %d",isalc[dk],ispan[dk]);
                                    goto err;
                                }

                                msym_salc_t *salc = &srs[dk].salc[isalc[dk]];
                                salc->d = sd;
                                // Add basis function and permutation irrep here
                                double (*salcpf)[dim] = calloc(salc->d,sizeof(*salcpf));
                                memcpy(salcpf, pf[pfi], sd*sizeof(*salcpf));
                                
                                salc->pf = (double*) salcpf;
                                salc->fl = dim;
                                salc->f = &srsbf[srsbfmap[i][n][l]];
                                
                                isalc[dk]++;
                            }
                            
                            pfi += sd;
                        }
                    }
                }
            }
        }
    }
    
    debug_printSubspace(ct,ct->d,srs);
    *ospan = ispan;
    *osrsl = ct->d;
    *osrs = srs;
    *osrsbf = srsbf;
    

    free(pmem);
    free(bspan);
    free(pspan);
    free(rspan);
    free(dspan);
    free(mspan);
    free(cmem);
    free(iespan);
    free(ipspan);
    free(ibspan);
    free(sspmem);
    free(psspmem);
    free(lssp);
    free(rsg);
    free(sgc);
    free(sgd);
    free(isalc);
    free(esnmax);
    free(esbfmap);
    for(int l = 0;l <= lmax;l++){
        free(lts[l].t);
    }
    free(lts);
    free(les);
    free(srsbfmap);
    
    return ret;
    
err:
    
    free(pmem);
    free(bspan);
    free(pspan);
    free(rspan);
    free(dspan);
    free(mspan);
    free(cmem);
    free(iespan);
    free(ipspan);
    free(ibspan);
    free(sspmem);
    free(psspmem);
    free(lssp);
    free(rsg);
    free(sgc);
    free(sgd);
    free(ispan);
    free(isalc);
    free(esnmax);
    free(esbfmap);
    for(int l = 0;l <= lmax;l++){
        free(lts[l].t);
    }
    free(lts);
    free(les);
    free(srsbfmap);
    for(int k = 0;k < ct->d;k++){
        for(int i = 0;i < srs[k].salcl;i++){
            free(srs[k].salc[i].pf);
        }
        free(srs[k].salc);
    }
    free(srs);
    free(srsbf);
    
    return ret;
}



msym_error_t generateSubrepresentationSpacesLowMem(msym_point_group_t *pg, int sgl, const msym_subgroup_t sg[sgl], int esl, msym_equivalence_set_t *es, msym_permutation_t **perm, int basisl, msym_basis_function_t basis[basisl], msym_element_t *elements, msym_equivalence_set_t **esmap, msym_thresholds_t *thresholds, int *osrsl, msym_subrepresentation_space_t **osrs, msym_basis_function_t ***osrsbf, int **ospan){
    msym_error_t ret = MSYM_SUCCESS;
    msym_character_table_t *ct = pg->ct;
    int lmax = -1, nmax = 0;
    enum _msym_basis_type ftype = basis[0].type;
    for(int i = 0;i < basisl;i++){
        if(basis[i].type != ftype) {nmax = -1; break;}
        lmax = basis[i].f.rsh.l > lmax ? basis[i].f.rsh.l : lmax;
        nmax = basis[i].f.rsh.n > nmax ? basis[i].f.rsh.n : nmax;
    }
    
    if(lmax < 0 || nmax < 1){
        if(nmax == -1) msymSetErrorDetails("Basis functions are not of the same type");
        else msymSetErrorDetails("Invalid sperical harmonics quantum numbers");
        ret = MSYM_INVALID_BASIS_FUNCTIONS;
        return ret;
    } else if (ftype != MSYM_BASIS_TYPE_REAL_SPHERICAL_HARMONIC){
        msymSetErrorDetails("Basis function type not supported");
        ret = MSYM_INVALID_BASIS_FUNCTIONS;
        return ret;
    }
    
    int projm = (2*lmax+1)*pg->order;
    
    double (*bspan)[ct->d] = calloc(lmax+1, sizeof(*bspan));        // span of individual basis functions
    double (*pspan)[ct->d] = calloc(esl, sizeof(*pspan));               // span of permutation operators
    
    double (*lspan)[ct->d] = calloc(esl, sizeof(*lspan));           // total span of basis function on each ES
    double *rspan = calloc(ct->d, sizeof(*rspan));                  // total direct product symmetrized basis
    double *dspan = calloc(ct->d, sizeof(*dspan));                  // decoposed total span of symmetrized basis (double)
    
    double *mspan = calloc(ct->d, sizeof(double));                  // span decomposition memory
    double (*mproj)[projm] = calloc(projm, sizeof(*mproj));         // projection operator memory
    double (*mscal)[projm] = calloc(projm, sizeof(*mscal));         // projection scaling memory
    double (*mpih)[projm] = calloc(projm, sizeof(*mpih));           // icosahedral projection memory
    double (*mperm)[pg->order] = calloc(pg->order, sizeof(*mperm)); // permutation memory
    double (*morth)[pg->order] = calloc(pg->order, sizeof(*morth)); // permutation orthoginalization memory
    double (*mbasis)[projm] = calloc(basisl, sizeof(*mbasis));      // basis function coefficients
    double (*mdec)[projm] = calloc(basisl, sizeof(*mdec));          // directo product decomposition memory
    double (*sgc)[5][pg->order] = calloc(ct->d,sizeof(*sgc));
    
    const msym_subgroup_t **rsg = calloc(ct->d, sizeof(*rsg));
    
    double *mdcomp = NULL;
    double *mdproj = NULL;
    int *mdfound = NULL;
    
    int (*sgd)[5] = calloc(ct->d,sizeof(*sgd));
    int *ispan = calloc(ct->d, sizeof(int));                               // decoposed total span of symmetrized basis (int)
    int *isalc = calloc(ct->d, sizeof(int));                               // number of added salcs to irrep
    int *esnmax = calloc(esl, sizeof(int));                                     // max n in eqset
    
    msym_basis_function_t *(*esbfmap)[pg->order][nmax+1][lmax+1][2*lmax+1] = calloc(esl,sizeof(*esbfmap));
    
    
    msym_basis_function_t *(*srsbf) = calloc(basisl, sizeof(*srsbf));
    
    int (*srsbfmap)[nmax+1][lmax+1] = calloc(esl,sizeof(*srsbfmap));
    
    rsh_representations_t *lts = calloc(lmax+1,sizeof(*lts)); // transformation matrices for rsh basis functions
    
    int (*les)[lmax+1] = calloc(esl, sizeof(*les));                      // number of l-type basis functions in each ES
    
    msym_basis_function_t dbf = {.type = ftype};
    double (*ctable)[ct->d] = ct->table;
    
    msym_subrepresentation_space_t *srs = calloc(ct->d, sizeof(*srs));
    
    /* determine number of l-type basis functions in each ES */
    for(int o = 0;o < basisl;o++){
        les[esmap[basis[o].element - elements] - es][basis[o].f.rsh.l] += basis[o].f.rsh.m == 0;
    }
    
    if(MSYM_SUCCESS != (ret = generateBasisRepresentations(pg->order+1, pg->order, pg->sops, lmax, lts))) goto err;
    
    
    for(int k = 0; k < ct->d;k++){
        if(ct->s[k].d > 1){
            if(MSYM_SUCCESS != (ret = findSplittingFieldSubgroup(pg, k, sgl, sg, thresholds, &rsg[k]))) goto err;
            if(MSYM_SUCCESS != (ret = getSplittingFieldCharacters(pg, rsg[k], sgc[k], sgd[k]))) goto err;
        }
    }
    
    
    for(int o = 0;o < basisl;o++){
        msym_basis_function_t *bf = &basis[o];
        int ei = (int)(bf->element - elements), esi = 0;
        msym_equivalence_set_t *e = esmap[ei];
        for(esi = 0;esi < e->length && e->elements[esi] != bf->element;esi++){}; //could improve perf with a map here
        if(esi >= e->length){
            ret = MSYM_INVALID_BASIS_FUNCTIONS;
            msymSetErrorDetails("Basis function does not map to any equivalence set");
            goto err;
        }
        if(bf->f.rsh.n > esnmax[e - es]) esnmax[e - es] = bf->f.rsh.n;
        esbfmap[e - es][esi][bf->f.rsh.n][bf->f.rsh.l][bf->f.rsh.m+bf->f.rsh.l] = bf;
    }
    
    int srsbfi = 0;
    
    for(int i = 0;i < esl;i++){
        for(int n = 0;n <= nmax;n++){
            for(int l = 0;l <= lmax;l++){
                srsbfmap[i][n][l] = srsbfi;
                for(int e = 0;e < es[i].length;e++){
                    msym_basis_function_t *lbf = esbfmap[i][e][n][l][0];
                    for(int m = -l;m <= l;m++){
                        msym_basis_function_t *bf = esbfmap[i][e][n][l][m+l];
                        if((NULL == bf && NULL != lbf) || (NULL == lbf && NULL != bf)) {
                            lbf = NULL == lbf ? bf : lbf;
                            ret = MSYM_INVALID_BASIS_FUNCTIONS;
                            msymSetErrorDetails("Found basis function %s but function where m = %d is missing on element %d of equivalence set %d",lbf->name,m,e,i);
                            goto err;
                        }
                        if(NULL != bf){
                            lbf = bf;
                            srsbf[srsbfi] = bf;
                            srsbfi++;
                        }
                    }
                }
                //when this is moved to a separate function it can be used and availability check (replace esbfmap)
                //if(srsbfmap[i][n][l] == srsbfi) srsbfmap[i][n][l] = -1;
            }
        }
    }
    
    if(srsbfi != basisl){
        ret = MSYM_SUBSPACE_ERROR;
        msymSetErrorDetails("Unexpected number of basis functions in subrepresentation map (expected %d, got %d)",basisl, srsbfi);
        goto err;
    }
    
    /* calculate span of irreducible representations for basis functions and permutations */
    for(int s = 0; s < pg->order;s++){
        for(int l = 0; l <= lmax;l++){
            dbf.f.rsh.l = l; //TODO: function type based
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
    
    /* scale reducible representations */
    for(int k = 0;k < ct->d;k++){
        double r = ct->s[k].r;
        if(r > 1){
            for(int l = 0; l <= lmax;l++) bspan[l][k] /= r;
            for(int i = 0;i < esl;i++) pspan[i][k] /= r;
        }
    }
    
    /* scale basis span, calculate basis function transforms and symmetry species vectors */
    for(int l = 0;l <= lmax;l++){
        int d = lts[l].d;
        double (*lproj)[d] = mproj;
        double (*lscal)[d] = mscal;
        
        vlscale(1.0/pg->order, ct->d, bspan[l], bspan[l]);
        
        double (*st)[d][d] = lts[l].t;
        memset(st[pg->order], 0, sizeof(double[d][d]));
        
        for(int k = 0, oirl = 0, nirl = 0;k < ct->d;k++, oirl = nirl){
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
                debug_printTransform(d, d, st[pg->order]);
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
        ispan[k] = (int)round(dspan[k]/ct->s[k].r);
        ddim_max = ddim_max > ispan[k] ? ddim_max : ispan[k];
        srs[k].s = k;
        srs[k].salcl = ispan[k];
        srs[k].salc = calloc(srs[k].salcl, sizeof(msym_salc_t));
    }
    
    mdcomp = malloc(sizeof(double[5][ddim_max][ddim_max]));
    mdproj = malloc(sizeof(double[projm][projm]));
    mdfound = malloc(sizeof(int[5][ddim_max]));
    
    
    clean_debug_printf("decomposed %d\n", ct->d);
    for(int prk = 0;prk < ct->d;prk++){
        if(prk < ct->d - 1) clean_debug_printf("%d%s + ", ispan[prk], ct->s[prk].name);
        else clean_debug_printf("%d%s\n", ispan[prk], ct->s[prk].name);        
    }
    
    int dbasisl = 0;
    
    for(int k = 0;k < ct->d;k++){
        dbasisl += ct->s[k].d*ispan[k];
    }
    
    if(dbasisl != basisl){
        ret = MSYM_SUBSPACE_ERROR;
        msymSetErrorDetails("Unexpected number of basis functions in decomposition (expected %d, got %d)",basisl,dbasisl);
        goto err;
    }
    
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
                    
                    for(int dk = 0; dk < ct->d; dk++) mspan[dk] /= ct->s[dk].r;
                    
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
                            //icosahedral symmetry map (dependant on getSplittingFieldCharacters)
                            int ihdim[] = {0,1,1,2,2};
                            int ihsub[] = {0,3,4,3,4};
                            double (*ihproj)[dd] = mpih;
                            double (*dproj)[dd] = mproj; // projection operator for direct product
                            double (*dscal)[dd] = mscal;
                            
                            memset(dproj, 0, sizeof(double[dd][dd]));
                            memset(mdec, 0, sizeof(double[vspan*lvspan][dd]));
                            for(int dim = 0, doirl = 0, dnirl = 0;dim < ct->s[sk].d;dim++, doirl = dnirl){
                                int pdim = ct->s[sk].d == 5 ? ihdim[dim] : dim; //icosahedral symmetry has split character table
                                memset(dproj, 0, sizeof(double[dd][dd]));
                                for(int s = 0;s < pg->order;s++){
                                    if(sgc[sk][pdim][s] == 0) continue;
                                    // we really should have precomputed this by now, it's the third time,
                                    // or a better approach would be a an optimized function
                                    permutationMatrix(&perm[i][s], mperm);
                                    kron(d, mperm, ld, lst[s], dd, dscal);
                                    mlscale(sgc[sk][pdim][s], dd, dscal, dscal);
                                    mladd(dd, dscal, dproj, dproj);
                                }
                                if(sgd[sk][pdim] == 2){ // only icosahedral
                                    int c2dim = ihsub[dim];
                                    memset(ihproj, 0, sizeof(double[dd][dd]));
                                    //getting rediculous
                                    for(int s = 0;s < pg->order;s++){
                                        if(sgc[sk][c2dim][s] == 0) continue;
                                        permutationMatrix(&perm[i][s], mperm);
                                        kron(d, mperm, ld, lst[s], dd, dscal);
                                        mlscale(sgc[sk][c2dim][s], dd, dscal, dscal);
                                        mladd(dd, dscal, ihproj, ihproj);
                                    }
                                    mmlmul(dd, dd, ihproj, dd, dproj, dscal); //avoid the mul malloc
                                    memcpy(dproj, dscal, sizeof(double[dd][dd]));
                                }
                                
                                
                                
                                mlscale(1.0/rsg[sk]->order, dd, dproj, dproj);
                                
                                //printf("constructed subgroup projection operator %s\n",cs_name[dim]);
                                //debug_printTransform(dd, dd, dproj);
                                
                                mmtlmul(dd,dd,dproj,svspan,&sbasis[si],dscal);
                                
                                memset(dproj, 0, sizeof(double[dd][dd]));
                                mltranspose(dd, svspan, dscal, dproj);
                                dnirl = mgs(dd, dproj, sdec, doirl, thresholds->orthogonalization/basisl);
                                if(dnirl - doirl != round(mspan[sk])){
                                    ret = MSYM_SUBSPACE_ERROR;
                                    msymSetErrorDetails("Multi-dimensional subspace decomposition of %s dimension (%d) inconsistent with representaion (%d) in subgroup irrep %d",ct->s[sk].name, dnirl - doirl,(int) round(mspan[sk]),dim);
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
                            for(int s = 0;s < pg->order ;s++){
                                permutationMatrix(&perm[i][s], mperm);
                                kron(d, mperm, ld, lst[s], dd, dscal);
                                for(int dim = 0;dim < mdim;dim++){
                                    for(int sg = 1;sg < ct->s[sk].d;sg++){
                                        if(found[sg][dim]) continue;
                                        mmtlmul(dd, dd, dscal, mdim, &sdec[sg*mdim], mt);
                                        mmlmul(1,dd,&sdec[dim], mdim, mt, &g[sg][dim]);
                                        if(vlabs(mdim, g[sg][dim]) > thresholds->zero) found[sg][dim] = 1;
                                    }
                                }
                            }
                            
                            for(int cdim = 0;cdim < mdim;cdim++){
                                vlnorm2(dd, sdec[cdim], sbasis[si+ct->s[sk].d*cdim]);
                            }
                            for(int dim = 1; dim < ct->s[sk].d;dim++){
                                mmlmul(mdim, mdim, g[dim], dd, &sdec[dim*mdim], dproj);
                                
                                for(int cdim = 0;cdim < mdim;cdim++){
                                    vlnorm2(dd, dproj[cdim], sbasis[si+dim+ct->s[sk].d*cdim]);
                                }
                            }
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
                                msym_salc_t *salc = &srs[sk].salc[isalc[sk]];
                                salc->d = ct->s[sk].d;
                                
                                double (*pf)[dd] = calloc(salc->d,sizeof(double[dd]));
                                
                                for(int dim = 0; dim < salc->d;dim++){
                                    vlnorm2(dd, sbasis[si+dim+ir], pf[dim]);
                                }
                                salc->pf = (double*) pf;
                                salc->fl = dd;
                                salc->f = &srsbf[srsbfmap[i][n][l]];
                                
                                isalc[sk]++;
                            }
                        }
                        si += svspan;
                    }
                }
            }
        }
    }
    
    debug_printSubspace(ct,ct->d,srs);
    *ospan = ispan;
    *osrsl = ct->d;
    *osrs = srs;
    *osrsbf = srsbf;
    
    free(bspan);
    free(pspan);
    free(lspan);
    free(rspan);
    free(dspan);
    free(mspan);
    free(mproj);
    free(mscal);
    free(mpih);
    free(mperm);
    free(morth);
    free(mbasis);
    free(mdec);
    free(rsg);
    free(sgc);
    free(sgd);
    free(mpih);
    free(mdcomp);
    free(mdproj);
    free(mdfound);
    free(isalc);
    free(esnmax);
    free(esbfmap);
    for(int l = 0;l <= lmax;l++){
        free(lts[l].t);
    }
    free(lts);
    free(les);
    free(srsbfmap);
    
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
    free(rsg);
    free(sgc);
    free(sgd);
    free(mpih);
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
    free(srsbfmap);
    for(int k = 0;k < ct->d;k++){
        for(int i = 0;i < srs[k].salcl;i++){
            free(srs[k].salc[i].pf);
        }
        free(srs[k].salc);
    }
    free(srs);
    free(srsbf);
    
    return ret;
}

void freeSubrepresentationSpaces(int srsl, msym_subrepresentation_space_t *srs){
    for(int i = 0;i < srsl && NULL != srs;i++){
        for(int j = 0;j < srs[i].salcl;j++){
            free(srs[i].salc[j].pf);
        }
        free(srs[i].salc);
    }
    free(srs);
}


