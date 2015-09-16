//
//  vibration.c
//  libmsym
//
//  Created by Marcus Johansson on 26/05/15.
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
#include "vibration.h"
#include "linalg.h"
#include "symop.h"
#include "permutation.h"
#include "point_group.h"
#include "subspace.h"

//void printVibrationsSubspace(CharacterTable *ct, msym_subspace_t *ss);


/*
msym_error_t generateDisplacementSubspaces(msym_point_group_t *pg, int esl, msym_equivalence_set_t *es, msym_permutation_t **perm, msym_thresholds_t *thresholds, int *subspacel, msym_subspace_t **subspace, int **pspan){
    printCharacterTable(pg->ct);
    msym_error_t ret = MSYM_SUCCESS;
    msym_subspace_t *iss = calloc(pg->ct->l, sizeof(msym_subspace_t));
    //double (*dspan) = malloc(sizeof(double[pg->ct->l]));
    int ssl = 0;
    int *issl = NULL;
    int (*ispan)[pg->ct->l] = calloc(esl,sizeof(int[pg->ct->l]));
    int (*aspan) = calloc(pg->ct->l,sizeof(int));
    double (*mkron)[3*pg->order] = malloc(sizeof(double[3*pg->order][3*pg->order]));
    double (*mperm)[pg->order] = malloc(sizeof(double[pg->order][pg->order]));
    double (*mproj)[3*pg->order][3*pg->order] = calloc(pg->ct->l+1,sizeof(double[3*pg->order][3*pg->order]));
    
    // Just precalculate span, makes things easier even though it may be slightly more calculations we can avoid alot of reallocs
    for(int k = 0;k < pg->ct->l;k++){
        iss[ssl].type = MASS_WEIGHTED_COORDINATES;
        iss[ssl].irrep = k;
        for(int i = 0;i < esl;i++){
            double dspan = 0;
            //memset(dspan,0,sizeof(double[pg->ct->l]));
            for(int s = 0; s < pg->order;s++){
                int uma = 0;
                for(int j = 0; j < perm[i][s].c_length;j++) uma += perm[i][s].c[j].l == 1;
                dspan += uma*pg->ct->irrep[k].v[pg->sops[s].cla]*symmetryOperationCartesianCharacter(&pg->sops[s]);
            }
            int ssvl = round(pg->ct->irrep[k].d*dspan/pg->order);
            iss[ssl].subspacel += (ssvl > 0);
            ispan[i][k] += ssvl;
            aspan[k] += ssvl;
        }
        
        if(iss[ssl].subspacel){
            iss[ssl].subspace = calloc(iss[ssl].subspacel, sizeof(msym_subspace_t));
            ssl++;
        }
        
        printf("span %s = %d\n",pg->ct->irrep[k].name,aspan[k]);
        for(int i = 0; i < esl;i++) printf("\tspan[%d] %s = %d\n",i,pg->ct->irrep[k].name,ispan[i][k]);
        
    }
    
    iss = realloc(iss, sizeof(msym_subspace_t[ssl]));
    issl = calloc(ssl, sizeof(int));
                  
    for(int i = 0;i < esl;i++){
        int d = 3*es[i].length;
        for(int s = 0;s < pg->order;s++){
            double m[3][3];
            permutationMatrix(&perm[i][s], mperm);
            symmetryOperationMatrix(&pg->sops[s], m);
            kron(perm[i][s].p_length,mperm,3,m,d,mkron);
            for(int k = 0;k < pg->ct->l;k++){
                mlscale(pg->ct->irrep[k].v[pg->sops[s].cla], d, mkron, mproj[pg->ct->l]); //mlproj[pg->ct->l] is a buffer
                mladd(d, mproj[pg->ct->l], mproj[k], mproj[k]);
                
            }
        }
        
        memset(mproj[pg->ct->l],0,sizeof(double[d][d]));
        
        int nirrepl = 0;
        printf("issl = %d\n",ssl);
        for(int is = 0;is < ssl;is++){
            int k = iss[is].irrep;
            int lirrepl = nirrepl;
            
            
            mlscale(((double) pg->ct->irrep[k].d)/pg->order, d, mproj[k], mproj[k]);
            nirrepl = mgs(d, mproj[k], mproj[pg->ct->l], nirrepl, thresholds->orthogonalization/d); //Use mkron as storage
            
            int nv = nirrepl - lirrepl;
            
            if(nv != ispan[i][k]){
                ret = MSYM_SUBSPACE_ERROR;
                msymSetErrorDetails("Ortogonal subspace of dimension (%d) inconsistent with span (%d) in %s",nv,ispan[i][k],pg->ct->irrep[k].name);
                goto err;
                
            }
            //printf("span[%d][%d] = %d subspace length = %d\n",i,k,ispan[i][k],iss[is].subspacel);
            if(ispan[i][k] > 0){
                if(issl[is] >= iss[is].subspacel){
                    ret = MSYM_SUBSPACE_ERROR;
                    msymSetErrorDetails("Number of subspaces (%d) larger than expected (%d) in %s",issl[is]+1,iss[is].subspacel,pg->ct->irrep[k].name);
                    goto err;
                }
                msym_subspace_t *ss = &iss[is].subspace[issl[is]];
                double (*dproj)[d] = mproj[pg->ct->l];
                printf("adding %d vectors to subspace %d in %s\n",nv,issl[is],pg->ct->irrep[k].name);
                ss->type = MASS_WEIGHTED_COORDINATES;
                ss->irrep = k;
                ss->d = ispan[i][k];
                ss->basisl = d;
                double (*space)[d] = malloc(sizeof(double[ss->d][d]));
                ss->space = (double*) space;
                //printTransform(d, d, dproj);
                
                for(int dim = 0; dim < ss->d;dim++){
                    vlnorm2(d, dproj[lirrepl+dim], space[dim]);
                }
                
                ss->basis.q = malloc(sizeof(msym_displacement_t[d]));
                double bv[3][3];
                mleye(3, bv);
                for(int j = 0;j < d;j++){
                    ss->basis.q[j].element = es[i].elements[j/3];
                    vcopy(bv[j % 3], ss->basis.q[j].v);
                }
                
                issl[is]++;
                
            }
        }
    }
    msym_subspace_t tss = {.subspacel = ssl, .subspace = iss, .d = 0, .basisl = 0, .space = NULL};
    printSubspace(pg->ct, &tss);
    
    *subspace = iss;
    *subspacel = ssl;
    *pspan = aspan;
    
    free(ispan);
    free(mproj);
    free(mperm);
    free(mkron);
    free(issl);
    
    return ret;
    
err:
    free(aspan);
    free(ispan);
    free(mproj);
    free(mperm);
    free(mkron);
    free(issl);

    for(int k = 0;k < ssl;k++){
        freeSubspace(&iss[k]);
    }
    free(iss);
    return ret;
    
}*/

/*
msym_error_t getDisplacementSubspaceCoefficients(msym_subspace_t *ss, int basisl, msym_orbital_t basis[basisl], int *offset, double c[basisl][basisl]){
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
            if(MSYM_SUCCESS != (ret = getDisplacementSubspaceCoefficients(&ss->subspace[i],basisl,basis,&index,c))) goto err;
        }
    }
    
    *offset = index;
    
    return ret;
err:
    return ret;
}*/

