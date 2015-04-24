//
//  permutation.c
//  Symmetry
//
//  Created by Marcus Johansson on 01/02/15.
//  Copyright (c) 2015 Marcus Johansson. 
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#include <stdlib.h>
#include <string.h>

#include "msym.h"
#include "permutation.h"
#include "linalg.h"


msym_error_t setPermutationCycles(msym_permutation_t *perm);

void freePermutationData(msym_permutation_t *perm){
    if(perm != NULL){
        free(perm->c);
        free(perm->p);
    }
}

msym_error_t findPermutation(msym_symmetry_operation_t *sop, int l, double (*v[l])[3], msym_thresholds_t *t, msym_permutation_t *perm){
    msym_error_t ret = MSYM_SUCCESS;
    double m[3][3];
    symmetryOperationMatrix(sop, m);
    
    perm->p = malloc(sizeof(int[l]));
    memset(perm->p, -1, sizeof(int[l]));
    perm->p_length = l;
    
    for(int i = 0; i < l;i++){
        int j;
        double r[3];
        mvmul(*v[i], m, r);
        for(j = 0;j < l;j++){
            if(vequal(r, *v[j],t->permutation)){
                perm->p[i] = j;
                break;
            }
        }
        if(j == l) {
            char buf[16];
            symmetryOperationName(sop, sizeof(buf), buf);
            msymSetErrorDetails("Unable to determine permutation for symmetry operation %s",buf);
            ret = MSYM_PERMUTATION_ERROR;
            goto err;
        }
    }
    if(MSYM_SUCCESS != (ret = setPermutationCycles(perm))) goto err;
    
    return ret;
    
err:
    free(perm->p);
    return ret;
}


//This is not functioning properly

msym_error_t findSymmetryOperationPermutations(int l, msym_symmetry_operation_t sops[l], msym_thresholds_t *t, msym_permutation_t **rperm){
    
    msym_error_t ret = MSYM_SUCCESS;
    //Don't block allocate this, it's a pain to keep track of the pointers
    msym_permutation_t *permutations = malloc(sizeof(msym_permutation_t[l]));
    
    for(int i = 0; i < l;i++){
        permutations[i].p = malloc(sizeof(int[l]));
        memset(permutations[i].p, -1, sizeof(int[l]));
        permutations[i].p_length = l;
    }
    
    double (*msops)[3][3] = malloc(sizeof(double[l][3][3]));
    
    for(int i = 0; i < l;i++){
        symmetryOperationMatrix(&sops[i], msops[i]);
    }
    
    for(int i = 0; i < l;i++){
        if((sops[i].type == PROPER_ROTATION && sops[i].order == 0) || sops[i].type == IDENTITY){
            for(int j = 0;j < l;j++) permutations[i].p[j] = j;
        } else {
            for(int j = 0; j < l;j++){
                int k;
                double rsop[3][3];
                mmmul(msops[i], msops[j], rsop);
                for(k = 0;k < l;k++){
                    if(mequal(rsop,msops[k],t->permutation)){
                        permutations[i].p[j] = k;
                        break;
                    }
                }
                if(k == l){
                    char buf1[16];
                    char buf2[16];
                    symmetryOperationName(&sops[i], sizeof(buf1), buf1);
                    symmetryOperationName(&sops[j], sizeof(buf2), buf2);
                    msymSetErrorDetails("Unable to determine permutation for symmetry operation %s and %s",buf1, buf2);
                    ret = MSYM_PERMUTATION_ERROR;
                    goto err;
                }
            }
        }
    }
    
    for(int i = 0; i < l;i++){
        if(MSYM_SUCCESS != (ret = setPermutationCycles(&permutations[i]))) goto err;
    }
    
    free(msops);
    *rperm = permutations;
    return ret;
    
err:
    free(msops);
    for(int i = 0; i < l;i++){
        free(permutations[i].p);
    }
    free(permutations);
    *rperm = NULL;
    return ret;
    
}

msym_error_t setPermutationCycles(msym_permutation_t *perm){
    msym_error_t ret = MSYM_SUCCESS;
    int l = perm->p_length;
    int *icycle = malloc(sizeof(int[l]));
    int *pcycle = malloc(sizeof(int[l]));
    int *lcycle = malloc(sizeof(int[l]));
    
    int cl = 0;
    memset(icycle, -1,sizeof(int[l]));
    memset(lcycle,  0,sizeof(int[l]));
    
    perm->c = NULL;
    perm->c_length = 0;
    
    for(int i = 0; i < l;i++){
        if(icycle[i] >= 0) continue;
        lcycle[cl] = 1;
        pcycle[cl] = i;
        icycle[i] = cl;
        for(int next = perm->p[i], loop = 0; next != i;next = perm->p[next]){
            if(loop++ > l) {
                msymSetErrorDetails("Encountered loop when determining permutation cycle");
                ret = MSYM_PERMUTATION_ERROR;
                goto err;
            }
            icycle[next] = cl;
            lcycle[cl]++;
        };
        cl++;
    }
    perm->c_length = cl;
    perm->c = malloc(sizeof(msym_permutation_cycle_t[cl]));
    for(int c = 0; c < cl;c++){
        perm->c[c].l = lcycle[c];
        perm->c[c].s = pcycle[c];
    }
    
err:
    free(icycle);
    free(pcycle);
    free(lcycle);
    return ret;
}


//We need these as doubles later, so might as well, even though we could represent these with ALOT less memory
void permutationMatrix(msym_permutation_t *perm, double m[perm->p_length][perm->p_length]){
    memset(m, 0, sizeof(double[perm->p_length][perm->p_length]));
    for(int i = 0;i < perm->p_length;i++){
        //m[i][perm->p[i]] = 1.0;
        m[perm->p[i]][i] = 1.0;
    }
}

void printPermutation(msym_permutation_t *perm){
    /*int l = perm->p_length;
    printf("(");
    for(int j = 0; j < l; j++){
        printf(j == l -1 ? "%d" : "%d\t",j);
    }
    printf(")\n(");
    for(int j = 0; j < l; j++){
        printf(j == l -1 ? "%d" : "%d\t",perm->p[j]);
    }
    printf(")\n");*/
    
    for(msym_permutation_cycle_t* c = perm->c; c < (perm->c + perm->c_length);c++){
        printf("(");
        for(int next = c->s, j = 0;j < c->l;j++){
            printf(j == c->l -1 ? "%d" : "%d ",next);
            next = perm->p[next];
        }
        printf(")");
    }
    
    printf("\n");
}