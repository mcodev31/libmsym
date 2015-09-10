//
//  point_group.c
//  Symmetry
//
//  Created by Marcus Johansson on 14/04/14.
//  Copyright (c) 2014 Marcus Johansson. 
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//


#include <math.h>
#include <float.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "symmetry.h"
#include "linalg.h"
#include "point_group.h"
#include "permutation.h"

#define PHI 1.618033988749894848204586834
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288419716939937510582
#endif

#ifndef M_PI_2
#define M_PI_2 (M_PI/2)
#endif


#define CLASSIFICATION_THRESHOLD 0.01



msym_error_t generatePointGroup(msym_point_group_type_t type, int n, msym_symmetry_operation_t *primary, int sopsl, msym_symmetry_operation_t sops[sopsl], msym_thresholds_t *thresholds, msym_point_group_t **opg);
msym_error_t reorientAxes(msym_symmetry_operation_t *primary, int sopsl, msym_symmetry_operation_t sops[sopsl], msym_thresholds_t *thresholds);
msym_error_t transformAxes(msym_point_group_type_t type, int n, msym_symmetry_operation_t *primary, int sopsl, msym_symmetry_operation_t sops[sopsl], msym_thresholds_t *thresholds, double transform[3][3]);
msym_error_t transformPrimary(msym_symmetry_operation_t *primary, int sopsl, msym_symmetry_operation_t sops[sopsl], msym_thresholds_t *thresholds, double transform[3][3]);
msym_error_t transformSecondary(msym_point_group_type_t type, msym_symmetry_operation_t *primary, int sopsl, msym_symmetry_operation_t sops[sopsl], msym_thresholds_t *thresholds, double transform[3][3]);
msym_error_t getPointGroupOrder(msym_point_group_type_t type, int n, int *order);
msym_error_t getPointGroupName(msym_point_group_type_t type, int n, size_t max, char name[max]);

msym_error_t findSecondaryAxisSigma(msym_symmetry_operation_t *primary, int sopsl, msym_symmetry_operation_t sops[sopsl], msym_thresholds_t *thresholds, double r[3]);
msym_error_t findSecondaryAxisC2(msym_symmetry_operation_t *primary, int sopsl, msym_symmetry_operation_t sops[sopsl], msym_thresholds_t *thresholds, double r[3]);

msym_error_t findSecondaryAxisC4(msym_symmetry_operation_t *primary, int sopsl, msym_symmetry_operation_t sops[sopsl], msym_thresholds_t *thresholds, double r[3]);

msym_error_t findSecondaryAxisC2C5(msym_symmetry_operation_t *primary, int sopsl, msym_symmetry_operation_t sops[sopsl], msym_thresholds_t *thresholds, double r[3]);



msym_error_t generateSymmetryOperations(msym_point_group_type_t type, int n, int order, msym_symmetry_operation_t **osops);


//msym_error_t generateSymmetryOperationsImpliedS(msym_point_group_t *pg, msym_thresholds_t *thresholds);
//msym_error_t generateSymmetryOperationsImpliedCPow(msym_point_group_t *pg, msym_thresholds_t *thresholds);
//msym_error_t generateSymmetryOperationsImpliedSPow(msym_point_group_t *pg, msym_thresholds_t *thresholds);
msym_error_t generateSymmetryOperationsImpliedRot(int sopsl, msym_symmetry_operation_t sops[sopsl], int order, msym_thresholds_t *thresholds, int *osopsl);

void generateSymmetryOperationsCiOld(msym_point_group_t *pg);
void generateSymmetryOperationsCsOld(msym_point_group_t *pg);
void generateSymmetryOperationsCnOld(msym_point_group_t *pg);
void generateSymmetryOperationsCnhOld(msym_point_group_t *pg);
void generateSymmetryOperationsCnvOld(msym_point_group_t *pg);
void generateSymmetryOperationsDnOld(msym_point_group_t *pg);
void generateSymmetryOperationsDnhOld(msym_point_group_t *pg);
void generateSymmetryOperationsDndOld(msym_point_group_t *pg);
void generateSymmetryOperationsS2nOld(msym_point_group_t *pg);
void generateSymmetryOperationsTOld(msym_point_group_t *pg);
void generateSymmetryOperationsTdOld(msym_point_group_t *pg);
void generateSymmetryOperationsThOld(msym_point_group_t *pg);
void generateSymmetryOperationsOOld(msym_point_group_t *pg);
void generateSymmetryOperationsOhOld(msym_point_group_t *pg);
void generateSymmetryOperationsIOld(msym_point_group_t *pg);
void generateSymmetryOperationsIhOld(msym_point_group_t *pg);

msym_error_t generateSymmetryOperationsDn(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsDnd(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsDnh(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsSn(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsCn(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsCnv(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsCnh(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsT(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsTd(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsO(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsOh(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsI(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsIh(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);

msym_error_t generateSymmetryOperationsOctahedral(int l, msym_symmetry_operation_t sops[l], int c2l, msym_symmetry_operation_t c2[c2l], int c3l, msym_symmetry_operation_t c3[c3l], int c4l, msym_symmetry_operation_t c4[c4l], int *pk);
msym_error_t generateSymmetryOperationsIcosahedral(int l, msym_symmetry_operation_t sops[l], int c2l, msym_symmetry_operation_t c2[c2l], int c3l, msym_symmetry_operation_t c3[c3l], int c5l, msym_symmetry_operation_t c5[c5l], int *pk);
msym_error_t generateReflectionPlanes(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateC2Axes(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);

msym_error_t generateSymmetryOperationsUnknown(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);


msym_error_t pointGroupFromName(char *name, msym_point_group_t *pg);

int classifySymmetryOperations(msym_point_group_t *pg);
void sortSymmetryOperations(msym_point_group_t *pg, int classes);

void print_transform(double M[3][3], double axis[3]);



msym_error_t generatePointGroupFromName(char *name, msym_thresholds_t *thresholds, msym_point_group_t **opg){
    msym_error_t ret = MSYM_SUCCESS;
    msym_point_group_t *pg = calloc(1,sizeof(msym_point_group_t));
    if(MSYM_SUCCESS != (ret = pointGroupFromName(name,pg))) goto err;
    if(MSYM_SUCCESS != (ret = generateSymmetryOperations(pg->type, pg->n, pg->order, &pg->sops))) goto err;
    //int classes = classifySymmetryOperations(pg);
    //sortSymmetryOperations(pg,classes);
    
    if((pg->type == POINT_GROUP_Cnv && pg->n == 0) || (pg->type == POINT_GROUP_Dnh && pg->n == 0)){
        pg->perm = NULL;
    } else {
        if(MSYM_SUCCESS != (ret = findSymmetryOperationPermutations(pg->order,pg->sops, thresholds, &pg->perm))) goto err;
    }
    for(msym_symmetry_operation_t *s = pg->sops;s < (pg->sops + pg->order);s++){
        if(pg->primary == NULL || (s->type == PROPER_ROTATION && s->order > pg->primary->order)) pg->primary = s;
    }
    
    mleye(3,pg->transform);
    
    *opg = pg;
    return ret;
    
err:
    *opg = NULL;
    free(pg);
    return ret;
}


msym_error_t pointGroupFromName(char *name, msym_point_group_t *pg){
    msym_error_t ret = MSYM_SUCCESS;
    int n = 0, gi = 0, ri = 0;;
    char g = 0, r = 0;
    
    int map[7][6];
    
    
    struct _pg_map {
        int i;
        msym_point_group_type_t type;
    } pg_map[] = {
        {1,  POINT_GROUP_Cn},
        {2,  POINT_GROUP_Cnv},
        {3,  POINT_GROUP_Cnh},
        {4,  POINT_GROUP_Ci},
        {5,  POINT_GROUP_Cs},
        {6,  POINT_GROUP_Dn},
        {7,  POINT_GROUP_Dnh},
        {8,  POINT_GROUP_Dnd},
        {9,  POINT_GROUP_Sn},
        {10, POINT_GROUP_T},
        {11, POINT_GROUP_Th},
        {12, POINT_GROUP_Td},
        {13, POINT_GROUP_O},
        {14, POINT_GROUP_Oh},
        {15, POINT_GROUP_I},
        {16, POINT_GROUP_Ih},
        {17, POINT_GROUP_K},
        {18, POINT_GROUP_Kh}
        
    };
    
    memset(map,0,sizeof(int[7][6]));
    
    map[0][0] =  1;
    map[0][1] =  2;
    map[0][2] =  3;
    map[0][4] =  4;
    map[0][5] =  5;
    map[1][0] =  6;
    map[1][2] =  7;
    map[1][3] =  8;
    map[2][0] =  9;
    map[3][0] =  10;
    map[3][2] =  11;
    map[3][3] =  12;
    map[4][0] =  13;
    map[4][2] =  14;
    map[5][0] =  15;
    map[5][2] =  16;
    map[6][0] =  17;
    map[6][2] =  18;
    
    if(sscanf(name,"%c%d%c",&g,&n,&r) < 2 && sscanf(name,"%c%c",&g,&r) < 1){
        ret = MSYM_INVALID_POINT_GROUP;
        msymSetErrorDetails("Invalid point group name %s",name);
        goto err;
    }
    
    if(n < 0) {
        ret = MSYM_INVALID_POINT_GROUP;
        msymSetErrorDetails("Invalid point group order %d",n);
        goto err;
    }
    
    switch(g){
        case 'C' : gi = 0; break;
        case 'D' : gi = 1; break;
        case 'S' : {
            if(n < 4 || n % 2 != 0){
                ret = MSYM_INVALID_POINT_GROUP;
                msymSetErrorDetails("Improper rotation order (%d) must be even",n);
                goto err;
            }
            gi = 2;
            break;
        }
        case 'T' : gi = 3; break;
        case 'O' : gi = 4; break;
        case 'I' : gi = 5; break;
        case 'K' : gi = 6; break;
        default :
            ret = MSYM_INVALID_POINT_GROUP;
            msymSetErrorDetails("Invalid point group type %c",g);
            goto err;
    }
    
    switch(r){
        case 0   : ri = 0; break;
        case 'v' : ri = 1; break;
        case 'h' : ri = 2; break;
        case 'd' : ri = 3; break;
        case 'i' : ri = 4; break;
        case 's' : ri = 5; break;
        default :
            ret = MSYM_INVALID_POINT_GROUP;
            msymSetErrorDetails("Invalid point group subtype %c",r);
            goto err;
    }
    
    int fi, fil = sizeof(pg_map)/sizeof(pg_map[0]);
    for(fi = 0;fi < fil;fi++){
        if(pg_map[fi].i == map[gi][ri]) break;
    }
    
    if(fi == fil){
        ret = MSYM_INVALID_POINT_GROUP;
        msymSetErrorDetails("Cannot find point group %s",name);
        goto err;
    }
    
    pg->type = pg_map[fi].type;
    pg->n = n;
    if(MSYM_SUCCESS != (ret = getPointGroupOrder(pg->type, pg->n, &pg->order))) goto err;
    if(MSYM_SUCCESS != (ret = getPointGroupName(pg->type, pg->n, sizeof(pg->name)/sizeof(char), pg->name))) goto err;
    
err:
    return ret;
    
}


msym_error_t getPointGroupName(msym_point_group_type_t type, int n, size_t max, char name[max]){
    msym_error_t ret = MSYM_SUCCESS;
    switch(type) {
        case POINT_GROUP_Ci  : snprintf(name,max,"Ci"); break;
        case POINT_GROUP_Cs  : snprintf(name,max,"Cs"); break;
        case POINT_GROUP_Cn  : snprintf(name,max,"C%d",n); break;
        case POINT_GROUP_Cnh : snprintf(name,max,"C%dh",n); break;
        case POINT_GROUP_Cnv : snprintf(name,max,"C%dv",n); break;
        case POINT_GROUP_Dn  : snprintf(name,max,"D%d",n); break;
        case POINT_GROUP_Dnh : snprintf(name,max,"D%dh",n); break;
        case POINT_GROUP_Dnd : snprintf(name,max,"D%dd",n); break;
        case POINT_GROUP_Sn  : snprintf(name,max,"S%d",n); break;
        case POINT_GROUP_T   : snprintf(name,max,"T"); break;
        case POINT_GROUP_Td  : snprintf(name,max,"Td"); break;
        case POINT_GROUP_Th  : snprintf(name,max,"Th"); break;
        case POINT_GROUP_O   : snprintf(name,max,"O"); break;
        case POINT_GROUP_Oh  : snprintf(name,max,"Oh"); break;
        case POINT_GROUP_I   : snprintf(name,max,"I"); break;
        case POINT_GROUP_Ih  : snprintf(name,max,"Ih"); break;
        case POINT_GROUP_K   : snprintf(name,max,"K"); break;
        case POINT_GROUP_Kh  : snprintf(name,max,"Kh"); break;
        default :
            msymSetErrorDetails("Unknown point group when determining name");
            ret = MSYM_POINT_GROUP_ERROR;
            goto err;
    }
err:
    return ret;
}

//init point group, copy all point groups so we can free the original list later
/*msym_error_t createPointGroup(msym_thresholds_t *thresholds,int n, msym_point_group_type_t type, msym_symmetry_operation_t *primary, msym_symmetry_operation_t *sops, unsigned int sopsl, msym_point_group_t **rpg){
    msym_error_t ret = MSYM_SUCCESS;
    
    msym_point_group_t pg = {.n = n, .type = type, .sops = sops, .sopsl = sopsl, .primary = primary, .ct = NULL};

    if(MSYM_SUCCESS != (ret = setPointGroupOrder(&pg))) goto err;
    if(MSYM_SUCCESS != (ret = setPointGroupName(sizeof(pg.name)/sizeof(char),n,type,pg.name))) goto err;
    if(MSYM_SUCCESS != (ret = symmetrizePointGroup(&pg, rpg, thresholds))) goto err;
    
    if(((*rpg)->type == POINT_GROUP_Cnv && (*rpg)->n == 0) || ((*rpg)->type == POINT_GROUP_Dnh && (*rpg)->n == 0)){
        (*rpg)->perm = NULL;
    } else {
        if(MSYM_SUCCESS != (ret = findSymmetryOperationPermutations((*rpg)->sopsl,(*rpg)->sops, thresholds, &(*rpg)->perm))) goto err;
    }
    
err:
    return ret;
    
}*/

msym_error_t generatePointGroup(msym_point_group_type_t type, int n, msym_symmetry_operation_t *primary, int sopsl, msym_symmetry_operation_t sops[sopsl], msym_thresholds_t *thresholds, msym_point_group_t **opg){
    msym_error_t ret = MSYM_SUCCESS;
    msym_point_group_t *pg = calloc(1,sizeof(msym_point_group_t));
    pg->type = type;
    pg->n = n;
    if(MSYM_SUCCESS != (ret = getPointGroupName(type,n,sizeof(pg->name)/sizeof(char),pg->name))) goto err;
    if(MSYM_SUCCESS != (ret = getPointGroupOrder(type, n, &pg->order))) goto err;
    if(pg->order < sopsl){
        ret = MSYM_POINT_GROUP_ERROR;
        msymSetErrorDetails("More symmetry operations than order of point group (%s). Order: %d Number of operations: %d",pg->name, pg->order,sopsl);
        goto err;
    }
    
    
    //Allocate more memory (allocate sops here, change generate)
    //if(MSYM_SUCCESS != (ret = generateSymmetryOperationsImpliedRot(sopsl, sops, pg->order, thresholds, &sopsl))) goto err;
    
    if(MSYM_SUCCESS != (ret = transformAxes(type, n, primary, sopsl, sops, thresholds, pg->transform))) goto err;
    
    if(MSYM_SUCCESS != (ret = generateSymmetryOperations(type, n, pg->order, &pg->sops))) goto err;
    
    if((pg->type == POINT_GROUP_Cnv && pg->n == 0) || (pg->type == POINT_GROUP_Dnh && pg->n == 0)){
        pg->perm = NULL;
    } else {
        if(MSYM_SUCCESS != (ret = findSymmetryOperationPermutations(pg->order,pg->sops, thresholds, &pg->perm))) goto err;
    }

    double T[3][3];
    minv(pg->transform, T);
    
    for(int i = 0; i < pg->order;i++){
        mvmul(pg->sops[i].v,T,pg->sops[i].v);
        if(pg->sops[i].type == PROPER_ROTATION){
            if(pg->primary == NULL || pg->sops[i].order > pg->primary->order) pg->primary = &(pg->sops[i]);
        }
    }
    
    for(int i = 0;i < pg->order;i++){
        printSymmetryOperation(&pg->sops[i]);
    }
    
    *opg = pg;
    return ret;
    
err:
    if(pg) free(pg->sops);
    free(pg);
    return ret;
}

msym_error_t getPointGroupOrder(msym_point_group_type_t type, int n, int *order){
    
    msym_error_t ret = MSYM_SUCCESS;
    switch (type) {
        case (POINT_GROUP_Cs)  :
        case (POINT_GROUP_Ci)  : *order = 2; break;
        case (POINT_GROUP_Cn)  :
        case (POINT_GROUP_Sn)  : *order = n; break;
        case (POINT_GROUP_Cnh) :
        case (POINT_GROUP_Dn)  : *order = 2*n; break;
        case (POINT_GROUP_Cnv) : *order = (n == 0 ? 2 : 2*n); break;
        case (POINT_GROUP_Dnh) : *order = (n == 0 ? 4 : 4*n); break;
        case (POINT_GROUP_Dnd) : *order = 4*n; break;
        case (POINT_GROUP_T)   : *order = 12; break;
        case (POINT_GROUP_Td)  :
        case (POINT_GROUP_Th)  :
        case (POINT_GROUP_O)   : *order = 24; break;
        case (POINT_GROUP_Oh)  : *order = 48; break;
        case (POINT_GROUP_I)   : *order = 60; break;
        case (POINT_GROUP_Ih)  : *order = 120; break;
        case (POINT_GROUP_K)   :
        case (POINT_GROUP_Kh)  : *order = 0;
        default                :
            msymSetErrorDetails("Point group has no primary axis for reorientation");
            goto err;
    }    
err:
    return ret;
}

msym_error_t findPointGroup(int sopsl, msym_symmetry_operation_t *sops, msym_thresholds_t *thresholds, msym_point_group_t **pg){
    msym_error_t ret = MSYM_SUCCESS;
    int inversion = 0, sigma = 0, nC[6] = {0,0,0,0,0,0}, linear = 0;
    msym_symmetry_operation_t *primary = NULL;
    msym_symmetry_operation_t *s = NULL;
    *pg = NULL;
    
    for(int i = 0;i < sopsl;i++){
        if(sops[i].type == PROPER_ROTATION && sops[i].order == 0){
            linear = 1;
            break;
        } else if (sops[i].type == PROPER_ROTATION && sops[i].order > 2){
            break;
        }
    }
    if(!linear){
        for(int i = 0;i < sopsl;i++){
            switch(sops[i].type){
                case PROPER_ROTATION :
                    if(primary == NULL || sops[i].order > primary->order) primary = &(sops[i]);
                    if(sops[i].order <= 5) nC[sops[i].order]++;
                    break;
                case INVERSION :
                    inversion = 1;
                    break;
                case REFLECTION :
                    sigma = 1;
                    break;
                case IMPROPER_ROTATION :
                    if(s == NULL || sops[i].order > s->order) s = &(sops[i]);
                    break;
                default :
                    break;
                    
            }
        }
        
        if(nC[3] >= 2) {
            if(nC[5] >= 2){
                if(inversion){
                    ret = generatePointGroup(POINT_GROUP_Ih, primary->order, primary, sopsl, sops, thresholds, pg);
                } else {
                    ret = generatePointGroup(POINT_GROUP_I, primary->order, primary, sopsl, sops, thresholds, pg);
                }
            } else if (nC[4] >= 2) {
                if(inversion){
                    ret = generatePointGroup(POINT_GROUP_Oh, primary->order, primary, sopsl, sops, thresholds, pg);
                } else {
                    ret = generatePointGroup(POINT_GROUP_O, primary->order, primary, sopsl, sops, thresholds, pg);
                }
                
            } else if (sigma){
                if(inversion){
                    ret = generatePointGroup(POINT_GROUP_Th, primary->order, primary, sopsl, sops, thresholds, pg);
                } else {
                    ret = generatePointGroup(POINT_GROUP_Td, primary->order, primary, sopsl, sops, thresholds, pg);
                }
            } else {
                ret = generatePointGroup(POINT_GROUP_T, primary->order, primary, sopsl, sops, thresholds, pg);
            }
            
        } else if (primary == NULL){
            if(sigma){
                ret = generatePointGroup(POINT_GROUP_Cs, 1, primary, sopsl, sops, thresholds, pg);
            } else if(inversion){
                ret = generatePointGroup(POINT_GROUP_Ci, 1, primary, sopsl, sops, thresholds, pg);
            } else {
                ret = generatePointGroup(POINT_GROUP_Cn, 1, primary, 0, NULL, thresholds, pg);
            }
        } else {
            int nC2 = 0;
            int sigma_h = 0;
            int nsigma_v = 0;
            
            //Special case for D2d
            if(primary->order == 2 && s != NULL && !vparallel(primary->v, s->v, thresholds->angle)){
                for(int i = 0; i < sopsl;i++){
                    if(sops[i].type == PROPER_ROTATION && sops[i].order == 2 && vparallel(sops[i].v, s->v, thresholds->angle)){
                        primary = &sops[i];
                        break;
                    }
                }
            }
            
            for(int i = 0; i < sopsl;i++){
                nC2 += sops[i].type == PROPER_ROTATION && sops[i].order == 2 && vperpendicular(primary->v,sops[i].v, thresholds->angle);
                sigma_h = sigma_h || (sops[i].type == REFLECTION && vparallel(primary->v,sops[i].v,thresholds->angle));
                nsigma_v += (sops[i].type == REFLECTION && vperpendicular(primary->v,sops[i].v,thresholds->angle));
            }
            if(nC2){ //actually nC2 == primary->order but less is acceptable here since we can generate the rest
                if(sigma_h){
                    ret = generatePointGroup(POINT_GROUP_Dnh, primary->order, primary, sopsl, sops, thresholds, pg);
                } else if (nsigma_v){ //actually nsigma_v == primary->order but less is acceptable here since we can generate the rest
                    ret = generatePointGroup(POINT_GROUP_Dnd, primary->order, primary, sopsl, sops, thresholds, pg);
                } else {
                    ret = generatePointGroup(POINT_GROUP_Dn, primary->order, primary, sopsl, sops, thresholds, pg);
                }
                
            } else if (sigma_h) {
                ret = generatePointGroup(POINT_GROUP_Cnh, primary->order, primary, sopsl, sops, thresholds, pg);
            } else if(nsigma_v) { //actually nsigma_v == primary->order but less is acceptable here since we can generate the rest
                ret = generatePointGroup(POINT_GROUP_Cnv, primary->order, primary, sopsl, sops, thresholds, pg);
            } else if(s != NULL){
                ret = generatePointGroup(POINT_GROUP_Sn, s->order, primary, sopsl, sops, thresholds, pg);
            } else {
                ret = generatePointGroup(POINT_GROUP_Cn, primary->order, primary, sopsl, sops, thresholds, pg);
            }
        }
    } else {
        for(int i = 0; i < sopsl;i++){
            inversion = inversion || sops[i].type == INVERSION;
            
            if(sops[i].type == PROPER_ROTATION && (sops[i].order == 0 || primary == NULL)){ //TODO The primary == NULL is needed because of a bug!
                primary = &(sops[i]);
            }
        }
        
        if(inversion){
            ret = generatePointGroup(POINT_GROUP_Dnh, primary->order, primary, sopsl, sops, thresholds, pg);
        } else {
            ret = generatePointGroup(POINT_GROUP_Cnv, primary->order, primary, sopsl, sops, thresholds, pg);
        }
    }
    
    return ret;
    
err:
    return ret;
    
}


msym_error_t findSubgroup(msym_subgroup_t *subgroup, msym_thresholds_t *thresholds){
    msym_error_t ret = MSYM_SUCCESS;
    int inversion = 0, sigma = 0, nC[6] = {0,0,0,0,0,0}, linear = 0;
    msym_symmetry_operation_t *primary = NULL;
    msym_symmetry_operation_t *s = NULL, *sop = NULL;;
    
    
    for(int i = 0;i < subgroup->sopsl;i++){
        if(subgroup->sops[i]->type == PROPER_ROTATION && subgroup->sops[i]->order == 0){
            linear = 1;
            break;
        } else if (subgroup->sops[i]->type == PROPER_ROTATION && subgroup->sops[i]->order > 2){
            break;
        }
    }
    if(!linear){
        for(int i = 0;i < subgroup->sopsl;i++){
            sop = subgroup->sops[i];
            if(sop->power > 1) continue;
            switch(subgroup->sops[i]->type){
                case PROPER_ROTATION :
                    if(primary == NULL || sop->order > primary->order) primary = sop;
                    if(sop->order <= 5) nC[sop->order]++;
                    break;
                case INVERSION :
                    inversion = 1;
                    break;
                case REFLECTION :
                    sigma = 1;
                    break;
                case IMPROPER_ROTATION :
                    if(s == NULL || sop->order > s->order) s = sop;
                    break;
                default :
                    break;
                    
            }
        }
        if(nC[3] >= 2) {
            if(nC[5] >= 2){
                if(inversion){
                    subgroup->type = POINT_GROUP_Ih;
                    subgroup->n = primary->order;
                } else {
                    subgroup->type = POINT_GROUP_I;
                    subgroup->n = primary->order;
                }
            } else if (nC[4] >= 2) {
                if(inversion){
                    subgroup->type = POINT_GROUP_Oh;
                    subgroup->n = primary->order;
                } else {
                    subgroup->type = POINT_GROUP_O;
                    subgroup->n = primary->order;
                }
                
            } else if (sigma){
                if(inversion){
                    subgroup->type = POINT_GROUP_Th;
                    subgroup->n = primary->order;
                } else {
                    subgroup->type = POINT_GROUP_Td;
                    subgroup->n = primary->order;
                }
            } else {
                subgroup->type = POINT_GROUP_T;
                subgroup->n = primary->order;
            }
            
        } else if (primary == NULL){
            if(sigma){
                subgroup->type = POINT_GROUP_Cs;
                subgroup->n = 1;
            } else if(inversion){
                subgroup->type = POINT_GROUP_Ci;
                subgroup->n = 1;
            } else {
                subgroup->type = POINT_GROUP_Cn;
                subgroup->n = 1;
            }
        } else {
            int nC2 = 0;
            int sigma_h = 0;
            int nsigma_v = 0;
            
            if(primary->order == 2 && s != NULL && !vparallel(primary->v, s->v, thresholds->angle)){
                for(int i = 0; i < subgroup->sopsl;i++){
                    sop = subgroup->sops[i];
                    if(sop->power > 1) continue;
                    if(sop->type == PROPER_ROTATION && sop->order == 2 && vparallel(sop->v, s->v, thresholds->angle)){
                        primary = sop;
                        break;
                    }
                }
            }
            
            for(int i = 0; i < subgroup->sopsl;i++){
                sop = subgroup->sops[i];
                if(sop->power > 1) continue;
                nC2 += sop->type == PROPER_ROTATION && sop->order == 2 && vperpendicular(primary->v,sop->v, thresholds->angle);
                sigma_h = sigma_h || (sop->type == REFLECTION && vparallel(primary->v,sop->v,thresholds->angle));
                nsigma_v += (sop->type == REFLECTION && vperpendicular(primary->v,sop->v,thresholds->angle));
            }
            if(nC2 == primary->order){
                if(sigma_h){
                    subgroup->type = POINT_GROUP_Dnh;
                    subgroup->n = primary->order;
                } else if (nsigma_v == primary->order){
                    subgroup->type = POINT_GROUP_Dnd;
                    subgroup->n = primary->order;
                } else {
                    subgroup->type = POINT_GROUP_Dn;
                    subgroup->n = primary->order;
                }
                
            } else if (sigma_h) {
                subgroup->type = POINT_GROUP_Cnh;
                subgroup->n = primary->order;
            } else if(nsigma_v == primary->order) {
                subgroup->type = POINT_GROUP_Cnv;
                subgroup->n = primary->order;
            } else if(s != NULL){
                subgroup->type = POINT_GROUP_Sn;
                subgroup->n = s->order;
            } else {
                subgroup->type = POINT_GROUP_Cn;
                subgroup->n = primary->order;
            }
        }
    } else {
        for(int i = 0; i < subgroup->sopsl;i++){
            inversion = inversion || subgroup->sops[i]->type == INVERSION;
            
            if(subgroup->sops[i]->type == PROPER_ROTATION && (subgroup->sops[i]->order == 0 || primary == NULL)){
                primary = subgroup->sops[i];
            }
        }
        
        if(inversion){
            subgroup->type = POINT_GROUP_Dnh;
            subgroup->n = primary->order;
        } else {
            subgroup->type = POINT_GROUP_Cnv;
            subgroup->n = primary->order;
        }
    }
    
    subgroup->primary = primary;
    if(MSYM_SUCCESS != (ret = getPointGroupName(subgroup->type, subgroup->n, sizeof(subgroup->name)/sizeof(char), subgroup->name))) goto err;
    return ret;
    
    
err:
    return ret;
    
}

msym_error_t transformAxes(msym_point_group_type_t type, int n, msym_symmetry_operation_t *primary, int sopsl, msym_symmetry_operation_t sops[sopsl], msym_thresholds_t *thresholds, double transform[3][3]){
    msym_error_t ret = MSYM_SUCCESS;
    switch (type){
        case (POINT_GROUP_Ci)  :
        case (POINT_GROUP_K)   :
        case (POINT_GROUP_Kh)  :
            break;
        case (POINT_GROUP_Cs)  :
            for(primary = sops; primary < (primary + sopsl) && primary->type != REFLECTION; primary++){};
        case (POINT_GROUP_Cn)  :
        case (POINT_GROUP_Cnh) :
        case (POINT_GROUP_Sn) :
            if(MSYM_SUCCESS != (ret = reorientAxes(primary, sopsl, sops, thresholds))) goto err;
            if(MSYM_SUCCESS != (ret = transformPrimary(primary, sopsl, sops, thresholds, transform))) goto err;
            break;
        case (POINT_GROUP_Cnv) :
        case (POINT_GROUP_Dnh) :
            if(MSYM_SUCCESS != (ret = reorientAxes(primary, sopsl, sops, thresholds))) goto err;
            if(MSYM_SUCCESS != (ret = transformPrimary(primary, sopsl, sops, thresholds, transform))) goto err;
            if(n > 0){
                if(MSYM_SUCCESS != (ret = transformSecondary(type, primary, sopsl, sops, thresholds, transform))) goto err;
            }
            break;
        case (POINT_GROUP_Dn)  :
        case (POINT_GROUP_Dnd) :
        case (POINT_GROUP_O)   :
        case (POINT_GROUP_Oh)  :
            if(MSYM_SUCCESS != (ret = reorientAxes(primary, sopsl, sops, thresholds))) goto err;
            if(MSYM_SUCCESS != (ret = transformPrimary(primary, sopsl, sops, thresholds, transform))) goto err;
            if(MSYM_SUCCESS != (ret = transformSecondary(type, primary, sopsl, sops, thresholds, transform))) goto err;
            break;
        case (POINT_GROUP_T)   :
        case (POINT_GROUP_Td)  :
        case (POINT_GROUP_Th)  :
        case (POINT_GROUP_I)   :
        case (POINT_GROUP_Ih)  : {
            msym_symmetry_operation_t *dprimary = NULL;
            msym_symmetry_operation_t *sop;
            for(sop = sops; sop < (sops + sopsl); sop++){
                double z = -2.0;
                if(sop->type == PROPER_ROTATION && sop->order == 2){
                    double v[3];
                    vnorm2(sop->v,v);
                    if(v[2] > z || dprimary == NULL) {
                        dprimary = sop;
                        z = sop->v[2];
                    }
                }
            }
            if(dprimary == NULL) {
                msymSetErrorDetails("Cannot find C2 axis for point group symmetrization of polyhedral point group");
                ret = MSYM_POINT_GROUP_ERROR;
                goto err;
            }
            if(MSYM_SUCCESS != (ret = reorientAxes(dprimary, sopsl, sops, thresholds))) goto err;
            if(MSYM_SUCCESS != (ret = transformPrimary(dprimary, sopsl, sops, thresholds, transform))) goto err;
            if(MSYM_SUCCESS != (ret = transformSecondary(type, dprimary, sopsl, sops, thresholds, transform))) goto err;
            break;
        }
            
    }
    
err:
    return ret;
}


msym_error_t pointGroupFromSubgroup(msym_subgroup_t *sg, msym_thresholds_t *thresholds, msym_point_group_t **opg){
    msym_error_t ret = MSYM_SUCCESS;
    *opg = calloc(1,sizeof(msym_point_group_t));
    msym_point_group_t *pg = *opg;
    pg->type = sg->type;
    pg->primary = sg->primary;
    pg->n = sg->n;
    pg->sops = malloc(sizeof(msym_symmetry_operation_t[sg->sopsl]));
    memcpy(pg->name,sg->name,sizeof(pg->name));
    
    if(MSYM_SUCCESS != (ret = getPointGroupOrder(pg->type, pg->n, &pg->order))) goto err;
    if(pg->order != sg->sopsl){
        msymSetErrorDetails("Point group order %d does not equal nuber of operations in subgroup %d for point gropu %s",pg->order,sg->sopsl,pg->name);
        ret = MSYM_POINT_GROUP_ERROR;
        goto err;
    }
    
    for(int i = 0;i < sg->sopsl;i++){
        if(sg->primary == sg->sops[i]) pg->primary = &pg->sops[i];
        memcpy(&pg->sops[i], sg->sops[i], sizeof(msym_symmetry_operation_t));
    }
    
    mleye(3, pg->transform);
    
    if(MSYM_SUCCESS != (ret = transformAxes(pg->type, pg->n, pg->primary, pg->order, pg->sops, thresholds, pg->transform))) goto err;
    
    /* Unfortunately we need to regenerate these as they need a specific
     * class ordering for orbital symmetrization */
    free(pg->sops);
    pg->sops = NULL;
    pg->primary = NULL;
    
    if(MSYM_SUCCESS != (ret = generateSymmetryOperations(pg->type, pg->n, pg->order, &pg->sops))) goto err;
    
    if((pg->type == POINT_GROUP_Cnv && pg->n == 0) || (pg->type == POINT_GROUP_Dnh && pg->n == 0)){
        pg->perm = NULL;
    } else {
        if(MSYM_SUCCESS != (ret = findSymmetryOperationPermutations(pg->order,pg->sops, thresholds, &pg->perm))) goto err;
    }
    
    double T[3][3];
    minv(pg->transform, T);
    
    for(int i = 0; i < pg->order;i++){
        mvmul(pg->sops[i].v,T,pg->sops[i].v);
    }
    
    return ret;
err:
    *opg = NULL;
    free(pg->sops);
    free(pg);
    return ret;
}

/*
msym_error_t symmetrizePointGroup(msym_point_group_t *ipg, msym_point_group_t **opg, msym_thresholds_t *thresholds){
    msym_error_t ret = MSYM_SUCCESS;
    
    *opg = calloc(1,sizeof(msym_point_group_t));
    msym_point_group_t *pg = *opg;
    if(MSYM_SUCCESS != (ret = copyPointGroup(ipg, pg))) goto err;
    if(MSYM_SUCCESS != (ret = generateSymmetryOperationsImpliedRot(pg, thresholds))) goto err;
    if(MSYM_SUCCESS != (ret = transformAxes(pg, thresholds))) goto err;
    
    free(pg->sops);
    pg->sops = NULL;
    pg->sopsl = 0;
    pg->primary = NULL;
    
    if(MSYM_SUCCESS != (ret = generateSymmetryOperations(pg,thresholds))) goto err;
    int classes = classifySymmetryOperations(pg);
    sortSymmetryOperations(pg,classes);
    
    //for(int i = 0;i < pg->sopsl;i++) printSymmetryOperation(&pg->sops[i]);
    
    double T[3][3];
    minv(pg->transform, T);
        
    for(int i = 0; i < pg->sopsl;i++){
        mvmul(pg->sops[i].v,T,pg->sops[i].v);
        if(pg->sops[i].type == PROPER_ROTATION){
            if(pg->primary == NULL || pg->sops[i].order > pg->primary->order) pg->primary = &(pg->sops[i]);
        }
    }
    
    return ret;
err:
    free(pg->sops);
    free(pg);
    *opg = NULL;
    return ret;
    
    
}*/


/* Point primary axes above xy plane and all symops above the primary plane.
 * Not really needed but we may not have to deal with transforms that flip things over
 * This could be improved with some thresholds, since some things lie in the plane
 * and get flipped. Should try to align with x as well
 */
msym_error_t reorientAxes(msym_symmetry_operation_t *primary, int sopsl, msym_symmetry_operation_t sops[sopsl], msym_thresholds_t *thresholds
){
    double x[3] = {1.0,0.0,0.0}, y[3] = {0.0,1.0,0.0}, z[3] = {0.0, 0.0, 1.0};
    
    if(primary == NULL) goto err;
    
    if(vdot(primary->v,z) < 0.0) vinv(primary->v);
    
    for(msym_symmetry_operation_t *sop = sops; sop < (sops + sopsl); sop++){
        if(vperpendicular(sop->v, z, thresholds->angle)){
            double proj = vdot(sop->v, y)/vabs(sop->v);
            if(fabs(fabs(proj)-1.0) < thresholds->angle && proj < 0.0) { //along y axis
                vinv(sop->v);
            } else if (vdot(sop->v,x) < 0.0){ //in xy-plane not in y axis, reorient to positive x
                vinv(sop->v);
            }
            
        } else if(vdot(primary->v,sop->v) < 0.0) vinv(sop->v); //not in xy-plane reorient to positive along primary axiz
    }
    
    return MSYM_SUCCESS;
    
err:
    msymSetErrorDetails("Point group has no primary axis for reorientation");
    return MSYM_POINT_GROUP_ERROR;
}




msym_error_t transformPrimary(msym_symmetry_operation_t *primary, int sopsl, msym_symmetry_operation_t sops[sopsl], msym_thresholds_t *thresholds, double transform[3][3]){
    msym_error_t ret = MSYM_SUCCESS;
    if(primary != NULL){
        double z[3] = {0.0, 0.0, 1.0};
        malign(primary->v,z,transform);
        for(msym_symmetry_operation_t *sop = sops; sop < (sops + sopsl); sop++){
            mvmul(sop->v,transform,sop->v);
        }
        vcopy(z,primary->v); // get rid of small errors
    } else {
        msymSetErrorDetails("Point group has no primary axis for transformation");
        ret = MSYM_POINT_GROUP_ERROR;
    }
    return ret;
}

msym_error_t transformSecondary(msym_point_group_type_t type, msym_symmetry_operation_t *primary, int sopsl, msym_symmetry_operation_t sops[sopsl], msym_thresholds_t *thresholds, double transform[3][3]){
    msym_error_t ret = MSYM_SUCCESS;
    double axis[3], x[3] = {1.0,0.0,0.0};
    
    switch(type){
        case (POINT_GROUP_Cnv) :
            if(MSYM_SUCCESS != (ret = findSecondaryAxisSigma(primary, sopsl, sops, thresholds, axis))) goto err;
            break;
        case (POINT_GROUP_O)   :
        case (POINT_GROUP_Oh)  :
            if(MSYM_SUCCESS != (ret = findSecondaryAxisC4(primary, sopsl, sops, thresholds, axis))) goto err;
            break;
        case (POINT_GROUP_Dn)  :
        case (POINT_GROUP_Dnh) :
        case (POINT_GROUP_Dnd) :
        case (POINT_GROUP_T)   :
        case (POINT_GROUP_Td)  :
        case (POINT_GROUP_Th)  :
            if(MSYM_SUCCESS != (ret = findSecondaryAxisC2(primary, sopsl, sops, thresholds, axis))) goto err;
            break;
        case (POINT_GROUP_I)   :
        case (POINT_GROUP_Ih)  :
            if(MSYM_SUCCESS != (ret = findSecondaryAxisC2C5(primary, sopsl, sops, thresholds, axis))) goto err;
            break;
        default :
            msymSetErrorDetails("Unknown point group when determining secondary axis");
            ret = MSYM_POINT_GROUP_ERROR;
            goto err;
    }
    
    double m[3][3];
    
    malign(axis, x, m);
    
    for(msym_symmetry_operation_t *sop = sops; sop < (sops + sopsl); sop++){
        mvmul(sop->v,m,sop->v);
    }
    mmmul(m,transform,transform);
    
err:
    return ret;
}



/* For point groups where we use a perpendicular reflection plane to indicate direction.
   We use a vector where the the xy-plane and reflection plane cross
*/
msym_error_t findSecondaryAxisSigma(msym_symmetry_operation_t *primary, int sopsl, msym_symmetry_operation_t sops[sopsl], msym_thresholds_t *thresholds, double r[3]){
    msym_error_t ret = MSYM_SUCCESS;
    msym_symmetry_operation_t *sop = NULL;
    for(sop = sops; sop < (sops + sopsl); sop++){
        if(sop->type == REFLECTION){
            vcross(sop->v, primary->v, r);
            vnorm(r);
            break;
        }
    }
    if(sop == (sops + sopsl)){
        msymSetErrorDetails("Can't find secondary reflection plane when symmetrizing point group");
        ret = MSYM_POINT_GROUP_ERROR;
        goto err;
    }
err:
    return ret;
}

/* For point groups where we use a perpendicular C2 axis to indicate direction.
   Adjusted to make sure it's perfectly in the xy-plane.
 */
msym_error_t findSecondaryAxisC2(msym_symmetry_operation_t *primary, int sopsl, msym_symmetry_operation_t sops[sopsl], msym_thresholds_t *thresholds, double r[3]){
    msym_error_t ret = MSYM_SUCCESS;
    msym_symmetry_operation_t *sop = NULL;
    for(msym_symmetry_operation_t *sop = sops; sop < (sops + sopsl); sop++){
        //Choose a C2 perpendicular to the primary axis, it'll make things relatively easy
        if(sop != primary && sop->type == PROPER_ROTATION && sop->order == 2 && vperpendicular(sop->v, primary->v,thresholds->angle)){
            vproj_plane(sop->v, primary->v, r);
            vnorm(r);
            break;
        }
    }
    if(sop == (sops + sopsl)){
        msymSetErrorDetails("Can't find secondary C2 axis when symmetrizing point group");
        ret = MSYM_POINT_GROUP_ERROR;
        goto err;
    }
err:
    return ret;
}


msym_error_t findSecondaryAxisC2C5(msym_symmetry_operation_t *primary, int sopsl, msym_symmetry_operation_t sops[sopsl], msym_thresholds_t *thresholds, double r[3]){
    msym_error_t ret = MSYM_SUCCESS;
    msym_symmetry_operation_t *c2[2], *c5 = NULL;
    int c2i = 0;

    for(msym_symmetry_operation_t *sop = sops; sop < (sops + sopsl) && (c5 == NULL || c2i < 2); sop++){
        if(vperpendicular(sop->v, primary->v,thresholds->angle)){
            if(sop->type == PROPER_ROTATION && sop->order == 2){
                //printf("Found perpendicular C2\n");
                c2[c2i++] = sop;
            } else if (sop->type == PROPER_ROTATION && sop->order == 5){
                //printf("Found perpendicular C5\n");
                c5 = sop;
            }
        }
    }
    
    if(c5 == NULL || c2i < 2) {
        msymSetErrorDetails("Can't find secondary C2 axis when symmetrizing point group: (%s %s)",c5 == NULL ? "C5 axis missing" : "", c2i < 2 ? "C2 axis missing" : "");
        ret = MSYM_POINT_GROUP_ERROR;
        goto err;
    }
    
    if(fabs(vdot(c5->v, c2[0]->v)) > fabs(vdot(c5->v, c2[1]->v))){
        vproj_plane(c2[0]->v, primary->v, r);
    } else {
        vproj_plane(c2[1]->v, primary->v, r);
    }
    
err:
    return ret;
}



// Same as C2
msym_error_t findSecondaryAxisC4(msym_symmetry_operation_t *primary, int sopsl, msym_symmetry_operation_t sops[sopsl], msym_thresholds_t *thresholds, double r[3]){
    msym_error_t ret = MSYM_SUCCESS;
    msym_symmetry_operation_t *sop = NULL;
    for(sop = sops; sop < (sops + sopsl); sop++){
        //Choose a C2 perpendicular to the primary axis, it'll make things relatively easy
        if(sop != primary && sop->type == PROPER_ROTATION && sop->order == 4 && vperpendicular(sop->v, primary->v,thresholds->angle)){
            vproj_plane(sop->v, primary->v, r);
            vnorm(r);
            break;
        }
    }
    if(sop == (sops + sopsl)){
        msymSetErrorDetails("Can't find secondary C4 axis when symmetrizing point group");
        ret = MSYM_POINT_GROUP_ERROR;
        goto err;
    }
err:
    return ret;
}



/*
msym_error_t generateSymmetryOperationsOld(msym_point_group_t *pg, msym_thresholds_t *thresholds){
    msym_error_t ret = MSYM_SUCCESS;
    double origo[3] = {0.0,0.0,0.0};

    pg->sops = malloc(sizeof(msym_symmetry_operation_t[pg->order+1]));
    vcopy(origo,pg->sops[0].v);
    
    pg->sops[0].type = IDENTITY;
    pg->sops[0].order = 0;
    pg->sopsl = 1;
    
    switch (pg->type){
        case (POINT_GROUP_Ci)  :
            generateSymmetryOperationsCi(pg);
            break;
        case (POINT_GROUP_Cs)  :
            generateSymmetryOperationsCs(pg);
            break;
        case (POINT_GROUP_Cn)  :
            generateSymmetryOperationsCn(pg);
            break;
        case (POINT_GROUP_Cnh) :
            generateSymmetryOperationsCnh(pg);
            break;
        case (POINT_GROUP_Sn) :
            generateSymmetryOperationsS2n(pg);
            break;
        case (POINT_GROUP_Cnv) :
            generateSymmetryOperationsCnv(pg);
            break;
        case (POINT_GROUP_Dn)  :
            generateSymmetryOperationsDn(pg);
            break;
        case (POINT_GROUP_Dnh) :
            generateSymmetryOperationsDnh(pg);
            break;
        case (POINT_GROUP_Dnd) :
            generateSymmetryOperationsDnd(pg);
            break;
        case (POINT_GROUP_T)   :
            generateSymmetryOperationsT(pg);
            break;
        case (POINT_GROUP_Td)  :
            generateSymmetryOperationsTd(pg);
            break;
        case (POINT_GROUP_Th)  :
            generateSymmetryOperationsTh(pg);
            break;
        case (POINT_GROUP_O)   :
            generateSymmetryOperationsO(pg);
            break;
        case (POINT_GROUP_Oh)  :
            generateSymmetryOperationsOh(pg);
            break;
        case (POINT_GROUP_I)   :
            generateSymmetryOperationsI(pg);
            break;
        case (POINT_GROUP_Ih)  :
            generateSymmetryOperationsIh(pg);
            break;
        case (POINT_GROUP_K)   :
        case (POINT_GROUP_Kh)  :
            pg->sops = NULL;
            pg->sopsl = 0;
            break;
        default :
            ret = MSYM_POINT_GROUP_ERROR;
            msymSetErrorDetails("Unknown point group when generating symmetry operations");
            goto err;
    }
    
    if(MSYM_SUCCESS != (ret = generateSymmetryOperationsImpliedCPow(pg,thresholds))) goto err;
    if(MSYM_SUCCESS != (ret = generateSymmetryOperationsImpliedS(pg,thresholds))) goto err;
    if(MSYM_SUCCESS != (ret = generateSymmetryOperationsImpliedSPow(pg,thresholds))) goto err;
    if(MSYM_SUCCESS != (ret = generateSymmetryOperationsImpliedRot(pg,thresholds))) goto err;
    
    if(pg->sopsl != pg->order){
        ret = MSYM_POINT_GROUP_ERROR;
        msymSetErrorDetails("Number of generated operations (%d) not equal to point group order (%d)",pg->sopsl,pg->order);
        goto err;
    }
    
    pg->sops = realloc(pg->sops,sizeof(msym_symmetry_operation_t[pg->order]));
    
    return ret;
    
err:
    free(pg->sops);
    pg->sops = NULL;
    return MSYM_POINT_GROUP_ERROR;
    
}*/
/*
msym_error_t generateSymmetryOperationsImpliedS(msym_point_group_t *pg, msym_thresholds_t *thresholds){
    double origo[3] = {0.0,0.0,0.0};
    int n = pg->sopsl;
    for(msym_symmetry_operation_t *sopi = pg->sops; sopi < (pg->sops + n); sopi++){
        if(sopi->type == REFLECTION){
            for(msym_symmetry_operation_t *sopj = pg->sops; sopj < (pg->sops + n) && pg->sopsl < pg->order; sopj++){
                if(sopj->type == PROPER_ROTATION && sopj->order == 2 && sopj->power == 1 && vparallel(sopi->v, sopj->v,thresholds->angle)){
                    pg->sops[pg->sopsl].type = INVERSION;
                    pg->sops[pg->sopsl].order = 0;
                    pg->sops[pg->sopsl].power = 1;
                    vcopy(origo, pg->sops[pg->sopsl].v);
                    pg->sopsl += !findSymmetryOperation(&(pg->sops[pg->sopsl]), pg->sops, pg->sopsl,thresholds);
                    if(pg->sopsl > pg->order) goto err;
                } else if (sopj->type == PROPER_ROTATION && sopj->power == 1 && sopj->order > 0 && vparallel(sopi->v, sopj->v,thresholds->angle)){
                    
                    copySymmetryOperation(&(pg->sops[pg->sopsl]), sopj);
                    pg->sops[pg->sopsl].type = IMPROPER_ROTATION;
                    pg->sopsl += !findSymmetryOperation(&(pg->sops[pg->sopsl]), pg->sops, pg->sopsl,thresholds);
                    if(pg->sopsl > pg->order) goto err;
                }
            }
        }
    }
    return MSYM_SUCCESS;
err:
    msymSetErrorDetails("Generation of implied symmetry operations by reflection resulted in more operations than point group order");
    return MSYM_POINT_GROUP_ERROR;
}

msym_error_t generateSymmetryOperationsImpliedCPow(msym_point_group_t *pg, msym_thresholds_t *thresholds){
    int n = pg->sopsl;
    for(msym_symmetry_operation_t *sop = pg->sops; sop < (pg->sops + n); sop++){
        if(sop->type == PROPER_ROTATION){
            for(int pow = 2; pow < sop->order && pg->sopsl < pg->order; pow++){
                symopPow(sop, pow, &(pg->sops[pg->sopsl]));
                pg->sopsl += !findSymmetryOperation(&(pg->sops[pg->sopsl]), pg->sops, pg->sopsl,thresholds);
                if(pg->sopsl > pg->order) goto err;
            }
        }
    }
    return MSYM_SUCCESS;
err:
    msymSetErrorDetails("Generation of implied proper rotations resulted in more operations than point group order");
    return MSYM_POINT_GROUP_ERROR;
}

msym_error_t generateSymmetryOperationsImpliedSPow(msym_point_group_t *pg, msym_thresholds_t *thresholds){
    int n = pg->sopsl;
    for(msym_symmetry_operation_t *sop = pg->sops; sop < (pg->sops + n); sop++){
        if(sop->type == IMPROPER_ROTATION){
            int mpow = sop->order % 2 == 1 ? 2*sop->order : sop->order;
            for(int j = 2; j < mpow; j++){
                symopPow(sop, j, &(pg->sops[pg->sopsl]));
                pg->sopsl += !findSymmetryOperation(&(pg->sops[pg->sopsl]), pg->sops, pg->sopsl,thresholds);
                if(pg->sopsl > pg->order) goto err;
            }
        }
    }
    return MSYM_SUCCESS;
err:
    msymSetErrorDetails("Generation of implied improper operations resulted in more operations than point group order");
    return MSYM_POINT_GROUP_ERROR;
}*/


msym_error_t generateSymmetryOperationsImpliedRot(int sopsl, msym_symmetry_operation_t sops[sopsl], int order, msym_thresholds_t *thresholds, int *osopsl){
    
    int isopsl = sopsl;
    for(msym_symmetry_operation_t *sopi = sops; sopi < (sops + sopsl) && isopsl < order; sopi++){
        if(sopi->type == PROPER_ROTATION){
            for(msym_symmetry_operation_t *sopj = sops; sopj < (sops + sopsl); sopj++){
                int istype = (sopj->type == REFLECTION || sopj->type == IMPROPER_ROTATION || (sopj->type == PROPER_ROTATION));
                if(sopi != sopj && istype && !vparallel(sopi->v, sopj->v,thresholds->angle)){
                    copySymmetryOperation(&(sops[isopsl]), sopj);
                    applySymmetryOperation(sopi,sops[isopsl].v,sops[isopsl].v);
                    isopsl += !findSymmetryOperation(&(sops[isopsl]), sops, isopsl,thresholds);
                    if(isopsl > order) goto err;
                }
            }
        }
    }
    *osopsl = isopsl;
    return MSYM_SUCCESS;
err:
    msymSetErrorDetails("Generation of implied symmetry operations by rotation resulted in more operations than point group order");
    return MSYM_POINT_GROUP_ERROR;
}

/*
void generateSymmetryOperationsCiOld(msym_point_group_t *pg){
    double origo[3] = {0.0,0.0,0.0};
    int n = pg->sopsl;
    
    vcopy(origo,pg->sops[n].v);
    pg->sops[n].type = INVERSION;
    pg->sops[n].order = 0;
    pg->sops[n].power = 1;
    n++;
    
    pg->sopsl = n;
}

void generateSymmetryOperationsCsOld(msym_point_group_t *pg){
    double z[3] = {0.0,0.0,1.0};
    int n = pg->sopsl;
    
    vcopy(z,pg->sops[n].v);
    pg->sops[n].type = REFLECTION;
    pg->sops[n].order = 0;
    pg->sops[n].power = 1;
    n++;
    
    pg->sopsl = n;
}

void generateSymmetryOperationsCnOld(msym_point_group_t *pg){
    double z[3] = {0.0,0.0,1.0};
    int n = pg->sopsl;
    
    //Only need to generate the Cn the other will come through pow
    vcopy(z,pg->sops[n].v);
    pg->sops[n].type = PROPER_ROTATION;
    pg->sops[n].order = pg->n;
    pg->sops[n].power = 1;
    
    n++;
    pg->sopsl = n;
}

void generateSymmetryOperationsCnhOld(msym_point_group_t *pg){
    generateSymmetryOperationsCnOld(pg);
    generateSymmetryOperationsCsOld(pg);

}

void generateSymmetryOperationsCnvOld(msym_point_group_t *pg){
    double y[3] = {0.0,1.0,0.0}, z[3] = {0.0,0.0,1.0};
    int n;
    
    generateSymmetryOperationsCnOld(pg);
    
    if(pg->n > 0){
        n = pg->sopsl;
        
        vcopy(y,pg->sops[n].v);
        pg->sops[n].type = REFLECTION;
        pg->sops[n].order = 0;
        pg->sops[n].power = 1;
        copySymmetryOperation(&(pg->sops[n+1]), &(pg->sops[n]));
        vrotate(M_PI/pg->n, pg->sops[n].v, z, pg->sops[n+1].v);
        n += 2;
        
        pg->sopsl = n;
    }
    
}

void generateSymmetryOperationsDnOld(msym_point_group_t *pg){
    double x[3] = {1.0,0.0,0.0}, z[3] = {0.0,0.0,1.0};
    int n;
    
    generateSymmetryOperationsCn(pg);
    
    n = pg->sopsl;
    
    vcopy(x,pg->sops[n].v);
    pg->sops[n].type = PROPER_ROTATION;
    pg->sops[n].order = 2;
    pg->sops[n].power = 1;
    copySymmetryOperation(&(pg->sops[n+1]), &(pg->sops[n]));
    vrotate(M_PI/pg->n, pg->sops[n].v, z, pg->sops[n+1].v);
    n += 2;
    
    pg->sopsl = n;
    
}

void generateSymmetryOperationsDnhOld(msym_point_group_t *pg){
    double x[3] = {1.0,0.0,0.0}, y[3] = {0.0,1.0,0.0}, z[3] = {0.0,0.0,1.0};
    int n;
    
    generateSymmetryOperationsCnhOld(pg);
    if(pg->n > 0){
        n = pg->sopsl;
        
        vcopy(x,pg->sops[n].v);
        pg->sops[n].type = PROPER_ROTATION;
        pg->sops[n].order = 2;
        pg->sops[n].power = 1;
        copySymmetryOperation(&(pg->sops[n+1]), &(pg->sops[n]));
        vrotate(M_PI/pg->n, pg->sops[n].v, z, pg->sops[n+1].v);
        n += 2;
        
        vcopy(y,pg->sops[n].v);
        pg->sops[n].type = REFLECTION;
        pg->sops[n].order = 0;
        pg->sops[n].power = 1;
        copySymmetryOperation(&(pg->sops[n+1]), &(pg->sops[n]));
        vrotate(M_PI/pg->n, pg->sops[n].v, z, pg->sops[n+1].v);
        n += 2;
        
        pg->sopsl = n;
    } else {
        n = pg->sopsl;
        
        pg->sops[n].type = INVERSION;
        pg->sops[n].order = 0;
        pg->sops[n].power = 1;
        
        pg->sopsl++;
    }
}

void generateSymmetryOperationsDndOld(msym_point_group_t *pg){
    double z[3] = {0.0,0.0,1.0}, x[3] = {1.0,0.0,0.0};
    int n;
    
    generateSymmetryOperationsDnOld(pg);
    
    n = pg->sopsl;
    
    vrotate(M_PI/(2*pg->n),x,z,pg->sops[n].v);
    vcrossnorm(pg->sops[n].v,z,pg->sops[n].v);
    //vcopy(x,pg->sops[n].v);
    pg->sops[n].type = REFLECTION;
    pg->sops[n].order = 0;
    pg->sops[n].power = 1;
    copySymmetryOperation(&(pg->sops[n+1]), &(pg->sops[n]));
    vrotate(M_PI/pg->n, pg->sops[n].v, z, pg->sops[n+1].v);
    n += 2;
    
    vcopy(z,pg->sops[n].v);
    pg->sops[n].type = IMPROPER_ROTATION;
    pg->sops[n].order = 2*pg->n;
    pg->sops[n].power = 1;
    n++;
    
    pg->sopsl = n;
    
}

void generateSymmetryOperationsS2nOld(msym_point_group_t *pg){
    double z[3] = {0.0,0.0,1.0};
    int n = pg->sopsl;
    
    vcopy(z,pg->sops[n].v);
    pg->sops[n].type = IMPROPER_ROTATION;
    pg->sops[n].order = pg->n;
    pg->sops[n].power = 1;
    n++;
    
    pg->sopsl = n;
}

void generateSymmetryOperationsTOld(msym_point_group_t *pg){
    double v[4][3] = { {1.0,1.0,1.0}, {-1.0,1.0,1.0}, {1.0,-1.0,1.0}, {-1.0,-1.0,1.0} };
    
    int n;
    
    pg->n = 2;
    generateSymmetryOperationsDnOld(pg);
    pg->n = 3;
    
    n = pg->sopsl;
    
    for (int i = 0; i < 4; n++,i++){
        vnorm(v[i]);
        vcopy(v[i], pg->sops[n].v);
        pg->sops[n].type = PROPER_ROTATION;
        pg->sops[n].order = 3;
        pg->sops[n].power = 1;
    }
    
    pg->sopsl = n;
    
}

void generateSymmetryOperationsTdOld(msym_point_group_t *pg){
    double v[3][3] = { {1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0} };
    double xy[3] = {1.0,1.0,0.0};
    
    int n;
    
    generateSymmetryOperationsTOld(pg);
    
    n = pg->sopsl;
    
    for (int i = 0; i < 3; n++,i++){
        vnorm(v[i]);
        vcopy(v[i], pg->sops[n].v);
        pg->sops[n].type = IMPROPER_ROTATION;
        pg->sops[n].order = 4;
        pg->sops[n].power = 1;
    }
    
    vnorm(xy);
    vcopy(xy, pg->sops[n].v);
    pg->sops[n].type = REFLECTION;
    pg->sops[n].order = 0;
    pg->sops[n].power = 1;
    n++;
    
    pg->sopsl = n;
    
}

void generateSymmetryOperationsThOld(msym_point_group_t *pg){
    double v[4][3] = { {1.0,1.0,1.0}, {-1.0,1.0,1.0}, {1.0,-1.0,1.0}, {-1.0,-1.0,1.0} };
    
    int n;
    
    pg->n = 2;
    generateSymmetryOperationsDnhOld(pg);
    pg->n = 3;
    
    n = pg->sopsl;
    
    for (int i = 0; i < 4; n++,i++){
        vnorm(v[i]);
        vcopy(v[i], pg->sops[n].v);
        pg->sops[n].type = IMPROPER_ROTATION;
        pg->sops[n].order = 6;
        pg->sops[n].power = 1;
    }
    
    pg->sopsl = n;
}

void generateSymmetryOperationsOOld(msym_point_group_t *pg){
    double v[3][3] = { {1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0} }, xy[3] = {1.0,1.0,0.0};
    int n;
    
    pg->n = 4;
    generateSymmetryOperationsTOld(pg);
    
    n = pg->sopsl;
    
    vnorm(xy);
    vcopy(xy, pg->sops[n].v);
    pg->sops[n].type = PROPER_ROTATION;
    pg->sops[n].order = 2;
    pg->sops[n].power = 1;
    n++;
    
    for (int i = 0; i < 3; n++,i++){
        vnorm(v[i]);
        vcopy(v[i], pg->sops[n].v);
        pg->sops[n].type = PROPER_ROTATION;
        pg->sops[n].order = 4;
        pg->sops[n].power = 1;
    }
    
    pg->sopsl = n;
    
    
}

void generateSymmetryOperationsOhOld(msym_point_group_t *pg){
    double v[3][3] = { {1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0} }, xy[3] = {1.0,1.0,0.0};
    int n;
    
    pg->n = 4;
    generateSymmetryOperationsThOld(pg);
    
    n = pg->sopsl;
    
    vnorm(xy);
    vcopy(xy, pg->sops[n].v);
    pg->sops[n].type = PROPER_ROTATION;
    pg->sops[n].order = 2;
    pg->sops[n].power = 1;
    n++;
    
    vcopy(xy, pg->sops[n].v);
    pg->sops[n].type = REFLECTION;
    pg->sops[n].order = 0;
    pg->sops[n].power = 1;
    n++;
    
    
    for (int i = 0; i < 3; n++,i++){
        vnorm(v[i]);
        vcopy(v[i], pg->sops[n].v);
        pg->sops[n].type = PROPER_ROTATION;
        pg->sops[n].order = 4;
        pg->sops[n].power = 1;
    }
    
    pg->sopsl = n;
}

void generateSymmetryOperationsIOld(msym_point_group_t *pg){
    double v[6][3] = { {PHI,1.0,0.0}, {-PHI,1.0,0.0}, {0.0,PHI,1.0}, {0.0,-PHI,1.0}, {1.0,0.0,PHI}, {1.0,0.0,-PHI} };
    
    int n;
    
    pg->n = 4;
    generateSymmetryOperationsTOld(pg);
    
    n = pg->sopsl;
    
    for (int i = 0; i < 6; n++,i++){
        vnorm(v[i]);
        vcopy(v[i], pg->sops[n].v);
        pg->sops[n].type = PROPER_ROTATION;
        pg->sops[n].order = 5;
        pg->sops[n].power = 1;
    }
    
    pg->sopsl = n;
    
}

void generateSymmetryOperationsIhOld(msym_point_group_t *pg){
    double v[6][3] = { {PHI,1.0,0.0}, {-PHI,1.0,0.0}, {0.0,PHI,1.0}, {0.0,-PHI,1.0}, {1.0,0.0,PHI}, {1.0,0.0,-PHI} };
    
    int n;
    
    pg->n = 4;
    generateSymmetryOperationsThOld(pg);
    
    n = pg->sopsl;
    
    for (int i = 0; i < 6; n++,i++){
        vnorm(v[i]);
        vcopy(v[i], pg->sops[n].v);
        pg->sops[n].type = PROPER_ROTATION;
        pg->sops[n].order = 5;
        pg->sops[n].power = 1;
    }
    
    for (int i = 0; i < 6; n++,i++){
        vnorm(v[i]);
        vcopy(v[i], pg->sops[n].v);
        pg->sops[n].type = IMPROPER_ROTATION;
        pg->sops[n].order = 10;
        pg->sops[n].power = 1;
    }
    
    pg->sopsl = n;
}


*/

/*

msym_error_t generatePointGroupTest(int n){
    msym_error_t ret = MSYM_SUCCESS;
    int l = 1, cla = 1;
    msym_symmetry_operation_t *sops = calloc(n << 2, sizeof(msym_symmetry_operation_t));
    if(MSYM_SUCCESS != (ret = generatePointGroupDnh(n,n << 2,sops,&l,&cla))) goto err;
    //if(MSYM_SUCCESS != (ret = generateSymmetryOperationsSn(n,n,sops,&l,&cla))) goto err;
    printf("Generated %d operations and %d classes n = %d\n",l,cla,n);
    for(int i = 0;i < l;i++){
        printSymmetryOperation(&sops[i]);
    }
    
    msym_character_table_t ct;
    
    if(MSYM_SUCCESS != (ret = new_characterTableCn(POINT_GROUP_Dnh, n, l, sops, &ct))) goto err;
    printf("Generated %d representations\n",ct.d);
    
    return ret;
err:
    printf("Error\n");
    return ret;
}*/


msym_error_t generateSymmetryOperations(msym_point_group_type_t type, int n, int order, msym_symmetry_operation_t **osops){
    msym_error_t ret = MSYM_SUCCESS;
    msym_symmetry_operation_t *sops = calloc(order, sizeof(msym_symmetry_operation_t));
    sops[0].cla = 0;
    sops[0].type = IDENTITY;
    sops[0].power = 1;
    sops[0].order = 1;
    sops[0].orientation = NONE;
    int cla = 1, gsopsl = 1;
    
    const struct _fmap {
        msym_point_group_type_t type;
        msym_error_t (*f)(int, int, msym_symmetry_operation_t *, int *, int *);
    } fmap[18] = {
        
        [ 0] = {POINT_GROUP_Ci,  generateSymmetryOperationsUnknown},
        [ 1] = {POINT_GROUP_Cs,  generateSymmetryOperationsUnknown},
        [ 2] = {POINT_GROUP_Cn,  generateSymmetryOperationsCn},
        [ 3] = {POINT_GROUP_Cnh, generateSymmetryOperationsCnh},
        [ 4] = {POINT_GROUP_Cnv, generateSymmetryOperationsCnv},
        [ 5] = {POINT_GROUP_Dn,  generateSymmetryOperationsDn},
        [ 6] = {POINT_GROUP_Dnh, generateSymmetryOperationsDnh},
        [ 7] = {POINT_GROUP_Dnd, generateSymmetryOperationsDnd},
        [ 8] = {POINT_GROUP_Sn,  generateSymmetryOperationsSn},
        [ 9] = {POINT_GROUP_T,   generateSymmetryOperationsT},
        [10] = {POINT_GROUP_Td,  generateSymmetryOperationsTd},
        [11] = {POINT_GROUP_Th,  generateSymmetryOperationsUnknown},
        [12] = {POINT_GROUP_O,   generateSymmetryOperationsO},
        [13] = {POINT_GROUP_Oh,  generateSymmetryOperationsOh},
        [14] = {POINT_GROUP_I,   generateSymmetryOperationsI},
        [15] = {POINT_GROUP_Ih,  generateSymmetryOperationsIh},
        [16] = {POINT_GROUP_K,   generateSymmetryOperationsUnknown},
        [17] = {POINT_GROUP_Kh,  generateSymmetryOperationsUnknown}
    };
    
    int fi, fil = sizeof(fmap)/sizeof(fmap[0]);
    for(fi = 0; fi < fil;fi++){
        if(fmap[fi].type == type) {
            if(MSYM_SUCCESS != (ret = fmap[fi].f(n,order,sops,&gsopsl,&cla))) goto err;
            break;
        }
    }
    
    if(fi == fil){
        msymSetErrorDetails("Unknown point group when generating symmetry operations");
        ret = MSYM_POINT_GROUP_ERROR;
        goto err;
    }
    
    if(gsopsl != order){
        msymSetErrorDetails("Generated incorrect number of symmetry operations %d != %d",gsopsl,order);
        ret = MSYM_INVALID_POINT_GROUP;
        goto err;
    }
    
    for(int i = 0;i < order;i++) printSymmetryOperation(&sops[i]);
    
    *osops = sops;
    return ret;
err:
    free(sops);
    return ret;

}

msym_error_t generateSymmetryOperationsUnknown(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msymSetErrorDetails("Generating symmetry operations for unknown point group");
    return MSYM_POINT_GROUP_ERROR;
}
msym_error_t generateSymmetryOperationsSn(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla, m = (n << (n & 1));
    double z[3] = {0.0,0.0,1.0};
    msym_symmetry_operation_t sn = {.type = IMPROPER_ROTATION, .order = n, .power = 1, .orientation = HORIZONTAL};
    vcopy(z,sn.v);
    if(k + m - 1 > l){ret = MSYM_POINT_GROUP_ERROR; msymSetErrorDetails("Too many operations when generating S%d symmetry operations",n); goto err;}
    for(int i = 1;i <= m >> 1;i++){
        int index = k + ((i-1) << 1);
        symopPow(&sn, i, &sops[index]);
        sops[index].cla = cla + i - 1;
        printf("i = %d m = %d index = %d ",i,m,index);
        printSymmetryOperation(&sops[index]);
    }
/*
    int ri = k + (((m >> 1)-1) << 1);
    sops[ri].cla = cla + (m >> 1) - 1;
    sops[ri].power = 1;
    sops[ri].p.orientation = HORIZONTAL;
    sops[ri].type = REFLECTION;
    vcopy(z,sops[ri].v);
    printf("replacing symmetry operation %d\n",ri);
    printSymmetryOperation(&sops[ri]);*/
    
    for(int i = 1;i < m >> 1;i++){
        int index = k + 1 + ((i-1) << 1);
        symopPow(&sn, m-i, &sops[index]);
        sops[index].cla = cla + i - 1;
        printf("i = %d m = %d index = %d ",i,m,index);
        printSymmetryOperation(&sops[index]);
        
    }
    k += m - 1;
    cla += m >> 1;
    
    printf("------ Sn %d operations %d classes------\n",k-*pk, cla-*pcla);
    *pk = k; *pcla = cla;
    
    return ret;
err:
    return ret;
}

msym_error_t generateSymmetryOperationsCn(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla;
    double z[3] = {0.0,0.0,1.0};
    msym_symmetry_operation_t cn = {.type = PROPER_ROTATION, .order = n, .power = 1, .orientation = HORIZONTAL};
    if(k + n - 1 > l){ret = MSYM_POINT_GROUP_ERROR; msymSetErrorDetails("Too many operations when generating C%d symmetry operations",n); goto err;}
    vcopy(z,cn.v);
    
    for(int i = 1;i <= (n >> 1);i++){
        int index = k + (i << 1) - 2;
        symopPow(&cn, i, &sops[index]);
        sops[index].cla = cla + (index >> 1);
        printf("i = %d index = %d ",i,index);
        printSymmetryOperation(&sops[index]);
    }
    
    for(int i = 1;i < (n >> 1) + (n & 1);i++){
        int index = k + (i << 1) - 1;
        symopPow(&cn, n-i, &sops[index]);
        sops[index].cla = cla + ((index - 1) >> 1);
        printf("i = %d index = %d ",i,index);
        printSymmetryOperation(&sops[index]);
    }
    
    k += n - 1;
    cla += n >> 1;
    
    printf("------ Cn %d operations %d classes------\n",k-*pk, cla-*pcla);
    *pk = k; *pcla = cla;
    
    return ret;
err:
    return ret;
}

msym_error_t generateSymmetryOperationsCnh(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla, s = 0;
    double z[3] = {0.0,0.0,1.0};
    msym_symmetry_operation_t cn = {.type = PROPER_ROTATION, .order = n, .power = 1, .orientation = HORIZONTAL};
    msym_symmetry_operation_t sn = {.type = IMPROPER_ROTATION, .order = n, .power = 1, .orientation = HORIZONTAL};
    if(k + (n << 1) - 1 > l){ret = MSYM_POINT_GROUP_ERROR; msymSetErrorDetails("Too many operations when generating C%dh symmetry operations",n); goto err;}
    vcopy(z,cn.v); vcopy(z,sn.v);
    printf("------ Cnh begin ------\n");
    for(s = n;s % 2 == 0;s = s >> 1){
        cn.order = s;
        printf("-------- 1s = %d ---------\n",s);
        for(int i = 1;i <= s >> 1;i += 2){
            int index = k + ((i >> 1) << 1);
            symopPow(&cn, i, &sops[index]);
            sops[index].cla = cla + (i >> 1);
            printf("i = %d s = %d index = %d ",i,s,index);
            printSymmetryOperation(&sops[index]);
        }
        printf("-------- 2s = %d ---------\n",s);
        for(int i = 1;i < s >> 1;i += 2){
            int index = k + 1 + ((i >> 1) << 1);
            symopPow(&cn, s-i, &sops[index]);
            sops[index].cla = cla + (i >> 1);
            printf("i = %d s = %d index = %d ",i,s,index);
            printSymmetryOperation(&sops[index]);
        }
        
        k += (s >> 1);
        cla += (s >> 2) + ((s >> 1) & 1);
        
        printf("-------- 3s = %d ---------\n",s);
        sn.order = s;
        for(int i = 1;i <= s >> 1;i += 2){
            int index = k + ((i >> 1) << 1);
            symopPow(&sn, i, &sops[index]);
            sops[index].cla = cla + (i >> 1);
            printf("i = %d s = %d index = %d ",i,s,index);
            printSymmetryOperation(&sops[index]);
        }
        printf("-------- 4s = %d ---------\n",s);
        for(int i = 1;i < s >> 1;i += 2){
            int index = k + 1 + ((i >> 1) << 1);
            symopPow(&sn, s-i, &sops[index]);
            sops[index].cla = cla + (i >> 1);
            printf("i = %d s = %d index = %d ",i,s,index);
            printSymmetryOperation(&sops[index]);
        }
        printf("-------- 5s = %d ---------\n",s);
        k += (s >> 1);
        cla += (s >> 2) + ((s >> 1) & 1);
        
    }
    
    printf("------ Cnh end ------\n");
    
    if(MSYM_SUCCESS != (ret = generateSymmetryOperationsSn(s,l,sops,&k,&cla))) goto err;
    //k += (s << (s & 1)) - 1;
    //cla = sops[k-1].cla + 1;
    
    
    printf("------ Cnh %d operations %d classes------\n",k-*pk, cla-*pcla);
    *pk = k; *pcla = cla;
    
    return ret;
err:
    return ret;
    
}

msym_error_t generateSymmetryOperationsCnv(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla;
    if(k + (n << 1) - 1 > l){ret = MSYM_POINT_GROUP_ERROR; msymSetErrorDetails("Too many operations when generating C%dv symmetry operations",n); goto err;}
    if(MSYM_SUCCESS != (ret = generateSymmetryOperationsCn(n,l,sops,&k,&cla))) goto err;
    //k += n-1;
    //cla = sops[k-1].cla + 1;
    if(MSYM_SUCCESS != (ret = generateReflectionPlanes(n,l,sops,&k,&cla))) goto err;
    //k += n;
    //cla = sops[k-1].cla + 1;
    printf("------ Cnv %d operations %d classes------\n",k-*pk, cla-*pcla);
    *pk = k; *pcla = cla;
    
    return ret;
err:
    return ret;
}

msym_error_t generateSymmetryOperationsDn(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla;
    if(k + (n << 1) - 1 > l){ret = MSYM_POINT_GROUP_ERROR; msymSetErrorDetails("Too many operations when generating D%d symmetry operations",n); goto err;}
    //k += (n << 1) - 1;
    //cla = sops[k-1].cla + 1;
    if(MSYM_SUCCESS != (ret = generateSymmetryOperationsCn(n,l,sops,&k,&cla))) goto err;
    //k += n;
    //cla = sops[k-1].cla + 1;
    if(MSYM_SUCCESS != (ret = generateC2Axes(n,l,sops,&k,&cla))) goto err;
    //k += n;
    //cla = sops[k-1].cla + 1;
    
    printf("\n------ Dn %d operations %d classes------\n",k-*pk, cla-*pcla);
    *pk = k; *pcla = cla;
    
    return ret;
err:
    return ret;
}

msym_error_t generateSymmetryOperationsDnh(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla;
    if(k + (n << 2) - 1 > l){ret = MSYM_POINT_GROUP_ERROR; msymSetErrorDetails("Too many operations when generating D%dh symmetry operations",n); goto err;}
    if(MSYM_SUCCESS != (ret = generateSymmetryOperationsCnh(n,l,sops,&k,&cla))) goto err;
    if(MSYM_SUCCESS != (ret = generateReflectionPlanes(n,l,sops,&k,&cla))) goto err;
    if(MSYM_SUCCESS != (ret = generateC2Axes(n,l,sops,&k,&cla))) goto err;
    printf("\n------ Dnh %d operations %d classes------\n",k-*pk, cla-*pcla);
    *pk = k; *pcla = cla;
    
    return ret;
err:
    return ret;
}

msym_error_t generateSymmetryOperationsDnd(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla;
    double x[3] = {1.0,0.0,0.0}, y[3] = {0.0,1.0,0.0}, z[3] = {0.0,0.0,1.0};
    msym_symmetry_operation_t sigma = {.type = REFLECTION, .orientation = DIHEDRAL, .power = 1};
    msym_symmetry_operation_t c2 = {.type = PROPER_ROTATION, .order = 2, .power = 1, .orientation = VERTICAL};
    
    if(k + (n << 2) - 1 > l){ret = MSYM_POINT_GROUP_ERROR; msymSetErrorDetails("Too many operations when generating D%dd symmetry operations",n); goto err;}
    
    vcopy(x,c2.v); vrotate(M_PI_2/n, y, z, sigma.v);
    
    if(MSYM_SUCCESS != (ret = generateSymmetryOperationsSn((n << 1),l,sops, &k, &cla))) goto err;
    
    for(int i = 0;i < n;i++){
        memcpy(&(sops[k+i]), &sigma, sizeof(msym_symmetry_operation_t));
        vrotate(i*M_PI/n, sigma.v, z, sops[k+i].v);
        sops[k+i].cla = cla;
    }
    k += n;
    cla += 1;
    
    for(int i = 0;i < n;i++){
        memcpy(&(sops[k+i]), &c2, sizeof(msym_symmetry_operation_t));
        vrotate(i*M_PI/n, c2.v, z, sops[k+i].v);
        sops[k+i].cla = cla;
    }
    k += n;
    cla += 1;
    
    
    printf("\n------ Dnd %d operations %d classes------\n",k-*pk, cla-*pcla);
    *pk = k; *pcla = cla;
    
    return ret;
err:
    return ret;
}

msym_error_t generateSymmetryOperationsT(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla;
    double v2[3][3] = { {1,0,0}, {0,1,0}, {0,0,1} };
    double v3[4][3] = { {1,1,1}, {-1,1,1}, {1,-1,1}, {-1,-1,1} };
    
    if(k + 11 > l){ret = MSYM_POINT_GROUP_ERROR; msymSetErrorDetails("Too many operations when generating T operations %d >= %d",k + 11,l); goto err;}
    
    for(int i = 0; i < 3; k++,i++){
        vnorm2(v2[i],sops[k].v);
        sops[k].type = PROPER_ROTATION;
        sops[k].order = 2;
        sops[k].power = 1;
        sops[k].orientation = NONE;
        sops[k].cla = cla;
    }
    
    cla += 1;
    
    for(int i = 0; i < 8; k++,i++){
        vnorm2(v3[i % 4],sops[k].v);
        sops[k].type = PROPER_ROTATION;
        sops[k].order = 3;
        sops[k].power = 1 + (i/4);
        sops[k].orientation = NONE;
        sops[k].cla = cla;
    }
    
    cla += 1;
    
    printf("\n------ T %d operations %d classes------\n",k-*pk, cla-*pcla);
    
    *pk = k; *pcla = cla;
    
    return ret;
err:
    return ret;

}

msym_error_t generateSymmetryOperationsTd(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla;
    const double v2[3][3] = { {1,0,0}, {0,1,0}, {0,0,1} };
    const double vs[6][3] = {
        { 1, 0, 1},
        { 0, 1, 1},
        {-1, 0, 1},
        { 0,-1, 1},
        { 1, 1, 0},
        { 1,-1, 0}
    };
    
    
    
    if(MSYM_SUCCESS != (ret = generateSymmetryOperationsT(n,l,sops, &k, &cla))) goto err;
    if(k + 12 > l){ret = MSYM_POINT_GROUP_ERROR; msymSetErrorDetails("Too many operations when generating Td operations %d >= %d",k + 12, l); goto err;}
    
    msym_symmetry_operation_t sigma = {.type = REFLECTION, .order = 1, .power = 1, .cla = cla, .orientation = NONE};

    for(int i = 0; i < 6; i++){
        memcpy(&sops[k+i], &sigma, sizeof(msym_symmetry_operation_t));
        vnorm2(vs[i],sops[k+i].v);
    }
    
    k += 6;
    cla += 1;
    
    msym_symmetry_operation_t s4 = {.type = IMPROPER_ROTATION, .order = 4, .power = 1, .cla = cla, .orientation = NONE};
    
    for(int i = 0; i < 3; i++){
        memcpy(&sops[k+i], &s4, sizeof(msym_symmetry_operation_t));
        memcpy(&sops[k+i+3], &s4, sizeof(msym_symmetry_operation_t));
        vnorm2(v2[i],sops[k+i].v);
        vnorm2(v2[i],sops[k+i+3].v);
        sops[k+i+3].power = 3;
    }
    
    k += 6;
    cla += 1;
    
    printf("\n------ Td %d operations %d classes------\n",k-*pk, cla-*pcla);
    
    *pk = k; *pcla = cla;
    
    return ret;
err:
    return ret;

}

msym_error_t generateSymmetryOperationsO(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla;
    
    msym_symmetry_operation_t c2[1] = {
        [0] = {.type = PROPER_ROTATION, .order = 2, .power = 1, .cla = cla, .orientation = VERTICAL}
    };
    
    msym_symmetry_operation_t c3[2] = {
        [0] = {.type = PROPER_ROTATION, .order = 3, .power = 1, .cla = cla+1, .orientation = NONE},
        [1] = {.type = PROPER_ROTATION, .order = 3, .power = 2, .cla = cla+1, .orientation = NONE}
    };
    
    msym_symmetry_operation_t c4[4] = {
        [0] = {.type = PROPER_ROTATION, .order = 2, .power = 1, .cla = cla+2, .orientation = NONE},
        [1] = {.type = PROPER_ROTATION, .order = 4, .power = 1, .cla = cla+3, .orientation = NONE},
        [2] = {.type = PROPER_ROTATION, .order = 4, .power = 3, .cla = cla+3, .orientation = NONE}
    };
    
    if(MSYM_SUCCESS != (ret = generateSymmetryOperationsOctahedral(l,sops, 1, c2, 2, c3, 3, c4, &k))) goto err;
    
    cla += 4;
    
    printf("------ O %d operations %d classes------\n",k-*pk, cla-*pcla);
    *pk = k; *pcla = cla;
    
    return ret;
err:
    return ret;
}

msym_error_t generateSymmetryOperationsOh(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla;
    
    msym_symmetry_operation_t c2[2] = {
        [0] = {.type = PROPER_ROTATION, .order = 2, .power = 1, .cla = cla, .orientation = VERTICAL},
        [1] = {.type = REFLECTION, .order = 1, .power = 1, .cla = cla+1, .orientation = DIHEDRAL}
    };
    
    msym_symmetry_operation_t c3[4] = {
        [0] = {.type = PROPER_ROTATION, .order = 3, .power = 1, .cla = cla+2, .orientation = NONE},
        [1] = {.type = PROPER_ROTATION, .order = 3, .power = 2, .cla = cla+2, .orientation = NONE},
        [2] = {.type = IMPROPER_ROTATION, .order = 6, .power = 1, .cla = cla+3, .orientation = NONE},
        [3] = {.type = IMPROPER_ROTATION, .order = 6, .power = 5, .cla = cla+3, .orientation = NONE}
    };
    
    msym_symmetry_operation_t c4[6] = {
        [0] = {.type = PROPER_ROTATION, .order = 2, .power = 1, .cla = cla+4, .orientation = NONE},
        [1] = {.type = PROPER_ROTATION, .order = 4, .power = 1, .cla = cla+5, .orientation = NONE},
        [2] = {.type = PROPER_ROTATION, .order = 4, .power = 3, .cla = cla+5, .orientation = NONE},
        [3] = {.type = IMPROPER_ROTATION, .order = 4, .power = 1, .cla = cla+6, .orientation = NONE},
        [4] = {.type = IMPROPER_ROTATION, .order = 4, .power = 3, .cla = cla+6, .orientation = NONE},
        [5] = {.type = REFLECTION, .order = 1, .power = 1, .cla = cla+7, .orientation = HORIZONTAL}
    };
    
    if(MSYM_SUCCESS != (ret = generateSymmetryOperationsOctahedral(l,sops, 2, c2, 4, c3, 6, c4, &k))) goto err;
    
    if(k + 1 > l){ret = MSYM_POINT_GROUP_ERROR; msymSetErrorDetails("Too many operations when generating Oh operations %d >= %d",k + 12, l); goto err;}
    
    sops[k].type = INVERSION;
    sops[k].power = 1;
    sops[k].order = 1;
    sops[k].orientation = NONE;
    sops[k].cla = cla+8;
    k++;
    
    cla += 9;
    
    printf("------ O %d operations %d classes------\n",k-*pk, cla-*pcla);
    *pk = k; *pcla = cla;
    
    return ret;
err:
    return ret;
}

msym_error_t generateSymmetryOperationsI(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla;
    
    msym_symmetry_operation_t c2[1] = {
        [0] = {.type = PROPER_ROTATION, .order = 2, .power = 1, .cla = cla, .orientation = NONE}
    };
    
    msym_symmetry_operation_t c3[2] = {
        [0] = {.type = PROPER_ROTATION, .order = 3, .power = 1, .cla = cla+1, .orientation = NONE},
        [1] = {.type = PROPER_ROTATION, .order = 3, .power = 2, .cla = cla+1, .orientation = NONE}
    };
    
    msym_symmetry_operation_t c4[4] = {
        [0] = {.type = PROPER_ROTATION, .order = 5, .power = 1, .cla = cla+2, .orientation = NONE},
        [1] = {.type = PROPER_ROTATION, .order = 5, .power = 4, .cla = cla+2, .orientation = NONE},
        [2] = {.type = PROPER_ROTATION, .order = 5, .power = 2, .cla = cla+3, .orientation = NONE},
        [3] = {.type = PROPER_ROTATION, .order = 5, .power = 3, .cla = cla+3, .orientation = NONE}
    };
    
    if(MSYM_SUCCESS != (ret = generateSymmetryOperationsIcosahedral(l,sops, 1, c2, 2, c3, 4, c4, &k))) goto err;
    
    cla += 4;
    
    printf("------ I %d operations %d classes------\n",k-*pk, cla-*pcla);
    *pk = k; *pcla = cla;
    
    return ret;
err:
    return ret;
}


msym_error_t generateSymmetryOperationsIh(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla;
    
    msym_symmetry_operation_t c2[2] = {
        [0] = {.type = PROPER_ROTATION, .order = 2, .power = 1, .cla = cla, .orientation = NONE},
        [1] = {.type = REFLECTION, .order = 1, .power = 1, .cla = cla+1, .orientation = NONE}
    };
    
    msym_symmetry_operation_t c3[4] = {
        [0] = {.type = PROPER_ROTATION, .order = 3, .power = 1, .cla = cla+2, .orientation = NONE},
        [1] = {.type = PROPER_ROTATION, .order = 3, .power = 2, .cla = cla+2, .orientation = NONE},
        [2] = {.type = IMPROPER_ROTATION, .order = 6, .power = 1, .cla = cla+3, .orientation = NONE},
        [3] = {.type = IMPROPER_ROTATION, .order = 6, .power = 5, .cla = cla+3, .orientation = NONE}
    };
    
    msym_symmetry_operation_t c5[8] = {
        [0] = {.type = PROPER_ROTATION, .order = 5, .power = 1, .cla = cla+4, .orientation = NONE},
        [1] = {.type = PROPER_ROTATION, .order = 5, .power = 4, .cla = cla+4, .orientation = NONE},
        [2] = {.type = PROPER_ROTATION, .order = 5, .power = 2, .cla = cla+5, .orientation = NONE},
        [3] = {.type = PROPER_ROTATION, .order = 5, .power = 3, .cla = cla+5, .orientation = NONE},
        [4] = {.type = IMPROPER_ROTATION, .order = 10, .power = 1, .cla = cla+6, .orientation = NONE},
        [5] = {.type = IMPROPER_ROTATION, .order = 10, .power = 9, .cla = cla+6, .orientation = NONE},
        [6] = {.type = IMPROPER_ROTATION, .order = 10, .power = 3, .cla = cla+7, .orientation = NONE},
        [7] = {.type = IMPROPER_ROTATION, .order = 10, .power = 7, .cla = cla+7, .orientation = NONE},
    };
    
    if(MSYM_SUCCESS != (ret = generateSymmetryOperationsIcosahedral(l,sops, 2, c2, 4, c3, 8, c5, &k))) goto err;
    
    if(k + 1 > l){ret = MSYM_POINT_GROUP_ERROR; msymSetErrorDetails("Too many operations when generating Ih operations %d >= %d",k + 12, l); goto err;}
    
    sops[k].type = INVERSION;
    sops[k].power = 1;
    sops[k].order = 1;
    sops[k].orientation = NONE;
    sops[k].cla = cla+8;
    k++;
    
    cla += 9;
    
    printf("------ Ih %d operations %d classes------\n",k-*pk, cla-*pcla);
    *pk = k; *pcla = cla;
    
    return ret;
err:
    return ret;
}

msym_error_t generateSymmetryOperationsTetrahedral(int l, msym_symmetry_operation_t sops[l], int c2l, msym_symmetry_operation_t c2[c2l], int csl, msym_symmetry_operation_t cs[csl], int c3l, msym_symmetry_operation_t c3[c3l], int *pk){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk;
    
    if(k + c2l*3 + csl*6 + c3l*4 > l){ret = MSYM_POINT_GROUP_ERROR; msymSetErrorDetails("Too many operations when generating tetrahedral operations %d >= %d",k + c2l*3 + csl*6 + c3l*4, l); goto err;}

    const double v2[3][3] = {
        { 1, 0, 0},
        { 0, 1, 0},
        { 0, 0, 1},
    };
    
    const double vs[6][3] = {
        { 1, 0, 1},
        { 0, 1, 1},
        {-1, 0, 1},
        { 0,-1, 1},
        { 1, 1, 0},
        { 1,-1, 0}
    };
    
    const double v3[4][3] = {
        {1,1,1},
        {-1,1,1},
        {1,-1,1},
        {-1,-1,1}
    };
    

    
    for(int i = 0;i < c2l;i++){
        for(int j = 0; j < 3;j++){
            memcpy(&sops[k], &c2[i], sizeof(msym_symmetry_operation_t));
            vnorm2(v2[j],sops[k].v);
            k++;
        }
    }
    
    for(int i = 0;i < csl;i++){
        for(int j = 0; j < 6;j++){
            memcpy(&sops[k], &cs[i], sizeof(msym_symmetry_operation_t));
            vnorm2(vs[j],sops[k].v);
            k++;
        }
    }
    
    for(int i = 0;i < c3l;i++){
        for(int j = 0; j < 4;j++){
            memcpy(&sops[k], &c3[i], sizeof(msym_symmetry_operation_t));
            vnorm2(v3[j],sops[k].v);
            k++;
        }
    }
    
    printf("------ Tetra %d operations 0 classes------\n",k-*pk);
    *pk = k;
    
    return ret;
err:
    return ret;
    
}

msym_error_t generateSymmetryOperationsOctahedral(int l, msym_symmetry_operation_t sops[l], int c2l, msym_symmetry_operation_t c2[c2l], int c3l, msym_symmetry_operation_t c3[c3l], int c4l, msym_symmetry_operation_t c4[c4l], int *pk){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk;
    
    if(k + c2l*6 + c3l*4 + c4l*3 > l){ret = MSYM_POINT_GROUP_ERROR; msymSetErrorDetails("Too many operations when generating octahedral operations %d >= %d",k + c2l*6 + c3l*4 + c4l*3, l); goto err;}
    
    const double v2[6][3] = {
        { 1, 0, 1},
        { 0, 1, 1},
        {-1, 0, 1},
        { 0,-1, 1},
        { 1, 1, 0},
        { 1,-1, 0}
    };
    
    const double v3[4][3] = {
        {1,1,1},
        {-1,1,1},
        {1,-1,1},
        {-1,-1,1}
    };
    
    const double v4[3][3] = {
        { 1, 0, 0},
        { 0, 1, 0},
        { 0, 0, 1},
    };
    
    for(int i = 0;i < c2l;i++){
        for(int j = 0; j < 6;j++){
            memcpy(&sops[k], &c2[i], sizeof(msym_symmetry_operation_t));
            vnorm2(v2[j],sops[k].v);
            k++;
        }
    }
    
    for(int i = 0;i < c3l;i++){
        for(int j = 0; j < 4;j++){
            memcpy(&sops[k], &c3[i], sizeof(msym_symmetry_operation_t));
            vnorm2(v3[j],sops[k].v);
            k++;
        }
    }
    
    for(int i = 0;i < c4l;i++){
        for(int j = 0; j < 3;j++){
            memcpy(&sops[k], &c4[i], sizeof(msym_symmetry_operation_t));
            vnorm2(v4[j],sops[k].v);
            k++;
        }
    }
    
    printf("------ Octa %d operations 0 classes------\n",k-*pk);
    *pk = k;
    
    return ret;
err:
    return ret;
    
}


msym_error_t generateSymmetryOperationsIcosahedral(int l, msym_symmetry_operation_t sops[l], int c2l, msym_symmetry_operation_t c2[c2l], int c3l, msym_symmetry_operation_t c3[c3l], int c5l, msym_symmetry_operation_t c5[c5l], int *pk){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk;
    
    if(k + c2l*15 + c3l*10 + c5l*6 > l){ret = MSYM_POINT_GROUP_ERROR; msymSetErrorDetails("Too many operations when generating icosahedral operations %d >= %d",k + c2l*15 + c3l*10 + c5l*6, l); goto err;}
    
    const double v2[15][3] = {
        {0,0,1},
        {1,0,0},
        {0,1,0},
        {PHI,-(PHI+1),1},
        {(PHI+1),1,-PHI},
        {1,PHI,(PHI+1)},
        {PHI,(PHI+1),1},
        {(PHI+1),-1,-PHI},
        {-1,PHI,-(PHI+1)},
        {(PHI+1),1,PHI},
        {1,PHI,-(PHI+1)},
        {-PHI,(PHI+1),1},
        {1,-PHI,-(PHI+1)},
        {PHI,(PHI+1),-1},
        {(PHI+1),-1,PHI}
    };
    
    const double v3[10][3] = {
        {1,1,1},
        {-1,1,1},
        {1,-1,1},
        {-1,-1,1},
        {(PHI+1),0,1},
        {0,-1,(PHI+1)},
        {-1,-(PHI+1),0},
        {-1,(PHI+1),0},
        {0,1,(PHI+1)},
        {(PHI+1),0,-1}
    };
    
    const double v5[6][3] = {
        {PHI,1,0},
        {-PHI,1,0},
        {0,PHI,1},
        {0,-PHI,1},
        {1,0,PHI},
        {1,0,-PHI}
    };
    
    for(int i = 0;i < c2l;i++){
        for(int j = 0; j < 15;j++){
            memcpy(&sops[k], &c2[i], sizeof(msym_symmetry_operation_t));
            vnorm2(v2[j],sops[k].v);
            k++;
        }
    }
    
    for(int i = 0;i < c3l;i++){
        for(int j = 0; j < 10;j++){
            memcpy(&sops[k], &c3[i], sizeof(msym_symmetry_operation_t));
            vnorm2(v3[j],sops[k].v);
            k++;
        }
    }
    
    for(int i = 0;i < c5l;i++){
        for(int j = 0; j < 6;j++){
            memcpy(&sops[k], &c5[i], sizeof(msym_symmetry_operation_t));
            vnorm2(v5[j],sops[k].v);
            k++;
        }
    }    
    
    printf("------ Icosa %d operations 0 classes------\n",k-*pk);
    *pk = k;
    
    return ret;
err:
    return ret;
    
}

msym_error_t generateReflectionPlanes(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla;
    double z[3] = {0.0,0.0,1.0}, y[3] = {0.0,1.0,0.0};
    msym_symmetry_operation_t sigma = {.type = REFLECTION, .power = 1};
    if(k + n > l){ret = MSYM_POINT_GROUP_ERROR; msymSetErrorDetails("Too many operations when generating reflection planes"); goto err;}
    vcopy(y,sigma.v);
    enum _msym_symmetry_operation_orientation orientation[2] = {VERTICAL, DIHEDRAL};
    for(int i = 0;i < n;i++){
        int e = 1 & ~n, ie = ((i & e)), index = k + (i >> e) + (ie ? (n >> 1) : 0); //power = 1 - (ie << 1),
        memcpy(&(sops[index]), &sigma, sizeof(msym_symmetry_operation_t));
        vrotate(i*M_PI/n, sigma.v, z, sops[index].v);
        //sops[index].power = power; // Used the finding of eigenvalues for character tables
        sops[index].orientation = orientation[ie];
        sops[index].cla = cla + ie;
    }
    
    k += n;
    cla += 1 << (~n & 1); //1 or 2 added classes
    printf("------ R %d operations %d classes------\n",k-*pk, cla-*pcla);
    *pk = k; *pcla = cla;
    
    return ret;
err:
    return ret;
}

msym_error_t generateC2Axes(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla;
    double z[3] = {0.0,0.0,1.0}, x[3] = {1.0,0.0,0.0};
    msym_symmetry_operation_t c2 = {.type = PROPER_ROTATION, .order = 2, .power = 1};
    if(k + n > l){ret = MSYM_POINT_GROUP_ERROR; msymSetErrorDetails("Too many operations when generating C2 axes"); goto err;}
    vcopy(x,c2.v);
    enum _msym_symmetry_operation_orientation orientation[2] = {VERTICAL, DIHEDRAL};
    for(int i = 0;i < n;i++){
        int e = 1 & ~n, ie = ((i & e)), index = k + (i >> e) + (ie ? (n >> 1) : 0);
        memcpy(&(sops[index]), &c2, sizeof(msym_symmetry_operation_t));
        vrotate(i*M_PI/n, c2.v, z, sops[index].v);
        //sops[index].power = power; // Used the finding of eigenvalues for character tables
        sops[index].orientation = orientation[ie];
        sops[index].cla = cla + ie;
    }
    
    k += n;
    cla += 1 << (~n & 1); //1 or 2 added classes
    printf("------ C2 %d operations %d classes------\n",k-*pk, cla-*pcla);
    *pk = k; *pcla = cla;
    
    return ret;
err:
    return ret;
}


int classifySymmetryOperations(msym_point_group_t *pg){
    int c = 1;
    double (*mop)[3][3] = malloc(sizeof(double[pg->order][3][3]));
    double (*imop)[3][3] = malloc(sizeof(double[pg->order][3][3]));
    
    //There may be a better way to do this
    for(int i = 0; i < pg->order;i++){
        if(pg->sops[i].type == IDENTITY){
            pg->sops[i].cla = 0;
        } else {
            pg->sops[i].cla = -1;
        }
        msym_symmetry_operation_t isop;
        invertSymmetryOperation(&(pg->sops[i]), &isop);
        symmetryOperationMatrix(&(pg->sops[i]), mop[i]);
        symmetryOperationMatrix(&isop, imop[i]);
    }
    
    for(int i = 0; i < pg->order;i++){
        if(pg->sops[i].cla < 0){
            pg->sops[i].cla = c;
            for(int j = 0; j < pg->order;j++){
                double m[3][3];
                mmmul(mop[i], imop[j], m);
                mmmul(mop[j],m,m);
                for(int k = 0; k < pg->order;k++){
                    if(mequal(mop[k],m,CLASSIFICATION_THRESHOLD)){ //Don't need to be dynamic, this is done on generated point groups (always the same)
                        pg->sops[k].cla = c;
                    }
                }
            }
            c++;
        }
    }
    
    free(mop);
    free(imop);
    
    return c;
    
    
}

//cant be bothereed writing an efficient sorting alg for this
void sortSymmetryOperations(msym_point_group_t *pg, int classes){
    msym_symmetry_operation_t *tmp = malloc(pg->order*sizeof(msym_symmetry_operation_t));
    int n = 0;
    
    for(int c = 0; c < classes;c++){
        for(int i = 0; i < pg->order;i++){
            if(pg->sops[i].cla == c){
                copySymmetryOperation(&tmp[n], &pg->sops[i]);
                n++;
            }
        }
    }
    memcpy(pg->sops, tmp,pg->order*sizeof(msym_symmetry_operation_t));

    free(tmp);
}

int numberOfSubgroups(msym_point_group_t *pg){
    
    int n = pg->n;
    int size = 0, ndiv = n >= 2, sdiv = 0, nodd = 0, sodd = 0, neven = 0, seven = 0;
    
    switch (pg->type) {
        case POINT_GROUP_Kh  : size = -1; break;
        case POINT_GROUP_K   : size = -1; break;
        case POINT_GROUP_Ih  : size = 162; break;
        case POINT_GROUP_I   : size = 57; break;
        case POINT_GROUP_Oh  : size = 96; break;
        case POINT_GROUP_O   : size = 28; break;
        case POINT_GROUP_Th  : size = 24; break;
        case POINT_GROUP_Td  : size = 28; break;
        case POINT_GROUP_T   : size = 9; break;
        case POINT_GROUP_Ci  : size = 0; break;
        case POINT_GROUP_Cs  : size = 0; break;
        default: {
            for(int i = 2; i < n;i++){
                if(n % i == 0){ ndiv++; sdiv += i; }
            }
            for(int i = 3; i < n;i += 2){
                if(n % i == 0){ nodd++; sodd += i;}
            }
            for(int i = 4; i <= n;i += 2){
                if(n % i == 0){ neven++; seven += (n << 1)/i; }
            }
            switch (pg->type) {
                case POINT_GROUP_Cnh : {
                    size = 2*ndiv;
                    if(n % 2 == 0){
                        int n2 = n >> 1;
                        for(int i = 2; i < n2;i++){
                            if(n2 % i == 0){ size++;}
                        }
                        size += 1 + (n2 >= 2);
                    }
                    break;
                }
                case POINT_GROUP_Dn  :
                case POINT_GROUP_Cnv : size = n + ndiv + sdiv; break;
                case POINT_GROUP_Cn  : size = ndiv - 1; break;
                case POINT_GROUP_Dnh : {
                    if(n % 2 == 0) size = 4*n + 2*ndiv + 3*sdiv + 4 + neven + seven;
                    else size = 3*(n+sdiv+1) + 2*ndiv;
                    break;
                }
                case POINT_GROUP_Dnd : {
                    if(n % 2 == 0) size = 2*n + 3 + ndiv + 2*sdiv + nodd + sodd;
                    else size = 3*(n+sdiv+1) + 2*ndiv;
                    break;
                }
                case POINT_GROUP_Sn : size = ndiv - 1; break;
                default : break;
            }
        }
    }
    
    return size;
}

msym_error_t findCharacterTable(msym_point_group_t *pg){
    printf("WARNING skipping old ct\n");
    return MSYM_SUCCESS;
    
    const struct _fmap {
        msym_point_group_type_t type;
        msym_error_t (*f)(int, CharacterTable*);
        
    } fmap[18] = {
        
        [ 0] = {POINT_GROUP_Ci, characterTableUnknown},
        [ 1] = {POINT_GROUP_Cs, characterTableUnknown},
        [ 2] = {POINT_GROUP_Cn, characterTableUnknown},
        [ 3] = {POINT_GROUP_Cnh,characterTableCnh},
        [ 4] = {POINT_GROUP_Cnv,characterTableCnv},
        [ 5] = {POINT_GROUP_Dn, characterTableUnknown},
        [ 6] = {POINT_GROUP_Dnh,characterTableDnh},
        [ 7] = {POINT_GROUP_Dnd,characterTableUnknown},
        [ 8] = {POINT_GROUP_Sn,characterTableUnknown},
        [ 9] = {POINT_GROUP_T,  characterTableUnknown},
        [10] = {POINT_GROUP_Td, characterTableTd},
        [11] = {POINT_GROUP_Th, characterTableUnknown},
        [12] = {POINT_GROUP_O,  characterTableUnknown},
        [13] = {POINT_GROUP_Oh, characterTableUnknown},
        [14] = {POINT_GROUP_I,  characterTableUnknown},
        [15] = {POINT_GROUP_Ih, characterTableIh},
        [16] = {POINT_GROUP_K,  characterTableUnknown},
        [17] = {POINT_GROUP_Kh, characterTableUnknown}
    };
    
    msym_error_t ret = MSYM_SUCCESS;
    
    CharacterTable *ct = malloc(sizeof(CharacterTable));
    
    int fi, fil = sizeof(fmap)/sizeof(fmap[0]);
    for(fi = 0; fi < fil;fi++){
        if(fmap[fi].type == pg->type) {
            if(MSYM_SUCCESS != (ret = fmap[fi].f(pg->n,ct))) goto err;
            break;
        }
    }
    
    if(fi == fil){
        msymSetErrorDetails("Unknown point group when finding character table");
        ret = MSYM_POINT_GROUP_ERROR;
        goto err;
    }
    
    ct = realloc(ct, sizeof(CharacterTable)+sizeof(int[ct->l])+ct->l*sizeof(*ct->name));

    ct->classc = (int*)(ct + 1);
    ct->name = (char (*)[6]) ((int *)ct->classc + ct->l);
        
    memset(ct->classc, 0, sizeof(int[ct->l]));
    memset(ct->name, 0, ct->l*sizeof(*(ct->name)));
    for(int i = 0; i < pg->order;i++){
        ct->classc[pg->sops[i].cla]++;
        symmetryOperationShortName(&pg->sops[i], sizeof(*(ct->name)), ct->name[pg->sops[i].cla]);
    }
    pg->ct = ct;
    return ret;
err:
    free(ct);
    return ret;
}

void printPointGroup(msym_point_group_t *pg){
    char buf[64];
    if(pg == NULL){
        printf("No point group\n");
        return;
    }
    printf("PointGroup %s (%d,%d)\nPrimary:\n",pg->name, pg->order, pg->order);
    if(pg->primary != NULL) {
        symmetryOperationName(pg->primary, 64, buf);
        printf("%s\n",buf);
    } else {
        printf("No primary rotation axis\n");
    }
    for(int i = 0; i < pg->order;i++){
        symmetryOperationName(&pg->sops[i], 64, buf);
        printf("\t%s\n",buf);
    }
}


void print_transform(double M[3][3], double axis[3]){
    
    fprintf(stderr,"M = \n");
    fprintf(stderr,"[[%lf, %lf, %lf], ",M[0][0],M[0][1],M[0][2]);
    fprintf(stderr,"[%lf, %lf, %lf], ",M[1][0],M[1][1],M[1][2]);
    fprintf(stderr,"[%lf, %lf, %lf]]\n",M[2][0],M[2][1],M[2][2]);
    
    double v[3];
    mvmul(axis,M,v);
    fprintf(stderr,"After transform:\n");
    fprintf(stderr,"[%lf, %lf, %lf]\n",v[0],v[1],v[2]);
    
}
