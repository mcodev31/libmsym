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

#include "debug.h"

#define PHI 1.618033988749894848204586834
#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288419716939937510582
#endif

#ifndef M_PI_2
#define M_PI_2 (M_PI/2)
#endif


#define CLASSIFICATION_THRESHOLD 0.01

msym_error_t determinePointGroup(int sopsl, msym_symmetry_operation_t *sops, msym_thresholds_t *thresholds, msym_point_group_t *pg);

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


msym_error_t generateSymmetryOperationsImpliedRot(int sopsl, msym_symmetry_operation_t sops[sopsl], int order, msym_thresholds_t *thresholds, int *osopsl);


msym_error_t generateSymmetryOperationsCs(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsCi(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsDn(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsDnd(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsDnh(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsSn(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsCn(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsCnv(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsCnh(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsT(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsTd(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsTh(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsO(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsOh(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsI(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateSymmetryOperationsIh(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);

msym_error_t generateSymmetryOperationsTetrahedral(int l, msym_symmetry_operation_t sops[l], int c2l, msym_symmetry_operation_t c2[c2l], int csl, msym_symmetry_operation_t cs[csl], int c3l, msym_symmetry_operation_t c3[c3l], int *pk);
msym_error_t generateSymmetryOperationsOctahedral(int l, msym_symmetry_operation_t sops[l], int c2l, msym_symmetry_operation_t c2[c2l], int c3l, msym_symmetry_operation_t c3[c3l], int c4l, msym_symmetry_operation_t c4[c4l], int *pk);
msym_error_t generateSymmetryOperationsIcosahedral(int l, msym_symmetry_operation_t sops[l], int c2l, msym_symmetry_operation_t c2[c2l], int c3l, msym_symmetry_operation_t c3[c3l], int c5l, msym_symmetry_operation_t c5[c5l], int *pk);
msym_error_t generateReflectionPlanes(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);
msym_error_t generateC2Axes(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);

msym_error_t generateSymmetryOperationsUnknown(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla);


msym_error_t pointGroupFromName(const char *name, msym_point_group_t *pg);
msym_error_t pointGroupFromType(msym_point_group_type_t type, int n, msym_point_group_t *pg);
msym_error_t generatePointGroupFromStruct(msym_point_group_t *pg, double transform[3][3], msym_thresholds_t *thresholds);

int classifySymmetryOperations(msym_point_group_t *pg);
void sortSymmetryOperations(msym_point_group_t *pg, int classes);

msym_error_t generatePointGroupFromType(msym_point_group_type_t type, int n, double transform[3][3], msym_thresholds_t *thresholds, msym_point_group_t **opg){
    msym_error_t ret = MSYM_SUCCESS;
    msym_point_group_t *pg = calloc(1,sizeof(msym_point_group_t));
    if(MSYM_SUCCESS != (ret = pointGroupFromType(type,n,pg))) goto err;
    if(MSYM_SUCCESS != (ret = generatePointGroupFromStruct(pg, transform, thresholds))) goto err;
    *opg = pg;
    return ret;
err:
    free(pg);
    return ret;
}

msym_error_t generatePointGroupFromName(const char *name, double transform[3][3], msym_thresholds_t *thresholds, msym_point_group_t **opg){
    msym_error_t ret = MSYM_SUCCESS;
    msym_point_group_t *pg = calloc(1,sizeof(msym_point_group_t));
    if(MSYM_SUCCESS != (ret = pointGroupFromName(name,pg))) goto err;
    if(MSYM_SUCCESS != (ret = generatePointGroupFromStruct(pg, transform, thresholds))) goto err;
    *opg = pg;
    return ret;
err:
    free(pg);
    return ret;
}

msym_error_t generatePointGroupFromStruct(msym_point_group_t *pg, double transform[3][3], msym_thresholds_t *thresholds){
    msym_error_t ret = MSYM_SUCCESS;

    if(MSYM_SUCCESS != (ret = generateSymmetryOperations(pg->type, pg->n, pg->order, &pg->sops))) goto err;
    
    if(isLinearPointGroup(pg)){
        pg->perm = NULL;
    } else {
        if(MSYM_SUCCESS != (ret = findSymmetryOperationPermutations(pg->order,pg->sops, thresholds, &pg->perm))) goto err;
    }
    
    memcpy(pg->transform, transform, sizeof(double[3][3]));
    
    double T[3][3];
    minv(pg->transform, T);
    
    for(msym_symmetry_operation_t *s = pg->sops;s < (pg->sops + pg->order);s++){
        if(pg->primary == NULL || (s->type == PROPER_ROTATION && s->order > pg->primary->order)) pg->primary = s;
        mvmul(s->v,T,s->v);
    }
    
    return ret;
    
err:
    free(pg->sops);
    pg->sops = NULL;
    // Need to free findSymmetryOperationPermutations if there is ever a possibility of an error after that call
    return ret;
}

msym_error_t pointGroupFromType(msym_point_group_type_t type, int n, msym_point_group_t *pg){
    msym_error_t ret = MSYM_SUCCESS;
    pg->type = type;
    switch (pg->type) {
        case MSYM_POINT_GROUP_TYPE_Cs:
        case MSYM_POINT_GROUP_TYPE_Ci:
            n = 1;
            break;
        case MSYM_POINT_GROUP_TYPE_T:
        case MSYM_POINT_GROUP_TYPE_Td:
        case MSYM_POINT_GROUP_TYPE_Th:
            n = 3;
            break;
        case MSYM_POINT_GROUP_TYPE_O:
        case MSYM_POINT_GROUP_TYPE_Oh:
            n = 4;
            break;
        case MSYM_POINT_GROUP_TYPE_I:
        case MSYM_POINT_GROUP_TYPE_Ih:
            n = 5;
            break;
        default:
            break;
    }
    pg->n = n;
    if(MSYM_SUCCESS != (ret = getPointGroupOrder(pg->type, pg->n, &pg->order))) goto err;
    if(MSYM_SUCCESS != (ret = getPointGroupName(pg->type, pg->n, sizeof(pg->name)/sizeof(char), pg->name))) goto err;

err:
    return ret;
}

msym_error_t pointGroupFromName(const char *name, msym_point_group_t *pg){
    msym_error_t ret = MSYM_SUCCESS;
    int n = 0, gi = 0, ri = 0;;
    char g = 0, r = 0;
    
    int map[7][6];
    
    
    struct _pg_map {
        int i;
        msym_point_group_type_t type;
    } pg_map[] = {
        {1,  MSYM_POINT_GROUP_TYPE_Cn},
        {2,  MSYM_POINT_GROUP_TYPE_Cnv},
        {3,  MSYM_POINT_GROUP_TYPE_Cnh},
        {4,  MSYM_POINT_GROUP_TYPE_Ci},
        {5,  MSYM_POINT_GROUP_TYPE_Cs},
        {6,  MSYM_POINT_GROUP_TYPE_Dn},
        {7,  MSYM_POINT_GROUP_TYPE_Dnh},
        {8,  MSYM_POINT_GROUP_TYPE_Dnd},
        {9,  MSYM_POINT_GROUP_TYPE_Sn},
        {10, MSYM_POINT_GROUP_TYPE_T},
        {11, MSYM_POINT_GROUP_TYPE_Th},
        {12, MSYM_POINT_GROUP_TYPE_Td},
        {13, MSYM_POINT_GROUP_TYPE_O},
        {14, MSYM_POINT_GROUP_TYPE_Oh},
        {15, MSYM_POINT_GROUP_TYPE_I},
        {16, MSYM_POINT_GROUP_TYPE_Ih},
        {17, MSYM_POINT_GROUP_TYPE_K},
        {18, MSYM_POINT_GROUP_TYPE_Kh}
        
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
                msymSetErrorDetails("Improper rotation order (%d) must be even and >= 4",n);
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
    
    return pointGroupFromType(pg_map[fi].type, n, pg);
err:
    return ret;
    
}


msym_error_t getPointGroupName(msym_point_group_type_t type, int n, size_t max, char name[max]){
    msym_error_t ret = MSYM_SUCCESS;
    switch(type) {
        case MSYM_POINT_GROUP_TYPE_Ci  : snprintf(name,max,"Ci"); break;
        case MSYM_POINT_GROUP_TYPE_Cs  : snprintf(name,max,"Cs"); break;
        case MSYM_POINT_GROUP_TYPE_Cn  : snprintf(name,max,"C%d",n); break;
        case MSYM_POINT_GROUP_TYPE_Cnh : snprintf(name,max,"C%dh",n); break;
        case MSYM_POINT_GROUP_TYPE_Cnv : snprintf(name,max,"C%dv",n); break;
        case MSYM_POINT_GROUP_TYPE_Dn  : snprintf(name,max,"D%d",n); break;
        case MSYM_POINT_GROUP_TYPE_Dnh : snprintf(name,max,"D%dh",n); break;
        case MSYM_POINT_GROUP_TYPE_Dnd : snprintf(name,max,"D%dd",n); break;
        case MSYM_POINT_GROUP_TYPE_Sn  : snprintf(name,max,"S%d",n); break;
        case MSYM_POINT_GROUP_TYPE_T   : snprintf(name,max,"T"); break;
        case MSYM_POINT_GROUP_TYPE_Td  : snprintf(name,max,"Td"); break;
        case MSYM_POINT_GROUP_TYPE_Th  : snprintf(name,max,"Th"); break;
        case MSYM_POINT_GROUP_TYPE_O   : snprintf(name,max,"O"); break;
        case MSYM_POINT_GROUP_TYPE_Oh  : snprintf(name,max,"Oh"); break;
        case MSYM_POINT_GROUP_TYPE_I   : snprintf(name,max,"I"); break;
        case MSYM_POINT_GROUP_TYPE_Ih  : snprintf(name,max,"Ih"); break;
        case MSYM_POINT_GROUP_TYPE_K   : snprintf(name,max,"K"); break;
        case MSYM_POINT_GROUP_TYPE_Kh  : snprintf(name,max,"Kh"); break;
        default :
            msymSetErrorDetails("Unknown point group when determining name");
            ret = MSYM_POINT_GROUP_ERROR;
            goto err;
    }
err:
    return ret;
}


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
    
    
    //TODO: this may still be needed, in broken symmetry
    //if(MSYM_SUCCESS != (ret = generateSymmetryOperationsImpliedRot(sopsl, sops, pg->order, thresholds, &sopsl))) goto err;
    
    if(MSYM_SUCCESS != (ret = transformAxes(type, n, primary, sopsl, sops, thresholds, pg->transform))) goto err;
    
    if(MSYM_SUCCESS != (ret = generateSymmetryOperations(type, n, pg->order, &pg->sops))) goto err;
    
    if(isLinearPointGroup(pg)){
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
        case (MSYM_POINT_GROUP_TYPE_Cs)  :
        case (MSYM_POINT_GROUP_TYPE_Ci)  : *order = 2; break;
        case (MSYM_POINT_GROUP_TYPE_Cn)  :
        case (MSYM_POINT_GROUP_TYPE_Sn)  : *order = n; break;
        case (MSYM_POINT_GROUP_TYPE_Cnh) :
        case (MSYM_POINT_GROUP_TYPE_Dn)  : *order = 2*n; break;
        case (MSYM_POINT_GROUP_TYPE_Cnv) : *order = (n == 0 ? 2 : 2*n); break;
        case (MSYM_POINT_GROUP_TYPE_Dnh) : *order = (n == 0 ? 4 : 4*n); break;
        case (MSYM_POINT_GROUP_TYPE_Dnd) : *order = 4*n; break;
        case (MSYM_POINT_GROUP_TYPE_T)   : *order = 12; break;
        case (MSYM_POINT_GROUP_TYPE_Td)  :
        case (MSYM_POINT_GROUP_TYPE_Th)  :
        case (MSYM_POINT_GROUP_TYPE_O)   : *order = 24; break;
        case (MSYM_POINT_GROUP_TYPE_Oh)  : *order = 48; break;
        case (MSYM_POINT_GROUP_TYPE_I)   : *order = 60; break;
        case (MSYM_POINT_GROUP_TYPE_Ih)  : *order = 120; break;
        case (MSYM_POINT_GROUP_TYPE_K)   :
        case (MSYM_POINT_GROUP_TYPE_Kh)  : *order = 0;
        default                :
            msymSetErrorDetails("Point group has no primary axis for reorientation");
            goto err;
    }    
err:
    return ret;
}

msym_error_t findPointGroup(int sopsl, msym_symmetry_operation_t *sops, msym_thresholds_t *thresholds, msym_point_group_t **pg){
    
    msym_error_t ret = MSYM_SUCCESS;
    
    msym_point_group_t tpg = {.type = MSYM_POINT_GROUP_TYPE_Kh, .n = 0, .primary = NULL};
    
    msym_symmetry_operation_t *tsops = NULL;
    
    int tsopsl = sopsl;
    
    if(MSYM_SUCCESS != (ret = determinePointGroup(sopsl, sops, thresholds, &tpg))) goto err;
    if(MSYM_SUCCESS != (ret = getPointGroupOrder(tpg.type, tpg.n, &tpg.order))) goto err;
    
    if(tpg.order < sopsl){
        int length = 2*sopsl > 121 ? 2*sopsl : 121;
        tsops = calloc(length, sizeof(msym_symmetry_operation_t));
        memcpy(tsops, sops, sizeof(msym_symmetry_operation_t[sopsl]));
        if(MSYM_SUCCESS != (ret = generateSymmetryOperationsImpliedRot(tsopsl, tsops, length, thresholds, &tsopsl))) goto err;
        tpg.primary = NULL;
        tpg.n = 0;
        tpg.type = MSYM_POINT_GROUP_TYPE_Kh;
        if(MSYM_SUCCESS != (ret = determinePointGroup(tsopsl, tsops, thresholds, &tpg))) goto err;
        if(MSYM_SUCCESS != (ret = getPointGroupOrder(tpg.type, tpg.n, &tpg.order))) goto err;
        
        if(tpg.order < tsopsl){
            char buf[4] = {0,0,0,0};
            getPointGroupName(tpg.type, tpg.n, sizeof(buf),buf);
            ret = MSYM_POINT_GROUP_ERROR;
            msymSetErrorDetails("More symmetry operations than order of point group (%s). Order: %d Number of operations: %d", buf, tpg.order,tsopsl);
            goto err;
        }
        
        if(MSYM_SUCCESS != (ret = generatePointGroup(tpg.type, tpg.n, tpg.primary, tsopsl, tsops, thresholds, pg))) goto err;
        
    } else {
        if(MSYM_SUCCESS != (ret = generatePointGroup(tpg.type, tpg.n, tpg.primary, sopsl, sops, thresholds, pg))) goto err;
    }
    
err:
    free(tsops);
    return ret;

}

msym_error_t determinePointGroup(int sopsl, msym_symmetry_operation_t *sops, msym_thresholds_t *thresholds, msym_point_group_t *pg){
    msym_error_t ret = MSYM_SUCCESS;
    int inversion = 0, sigma = 0, nC[6] = {0,0,0,0,0,0}, linear = 0;
    msym_symmetry_operation_t *primary = NULL;
    msym_symmetry_operation_t *s = NULL;
    
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
        
        pg->n = NULL == primary ? 1 : primary->order;
        pg->primary = primary;
        
        if(nC[3] >= 2) {
            if(nC[5] >= 2){
                if(inversion){
                    pg->type = MSYM_POINT_GROUP_TYPE_Ih;
                } else {
                    pg->type = MSYM_POINT_GROUP_TYPE_I;
                }
            } else if (nC[4] >= 2) {
                if(inversion){
                    pg->type = MSYM_POINT_GROUP_TYPE_Oh;
                } else {
                    pg->type = MSYM_POINT_GROUP_TYPE_O;
                }
                
            } else if (sigma){
                if(inversion){
                    pg->type = MSYM_POINT_GROUP_TYPE_Th;
                } else {
                    pg->type = MSYM_POINT_GROUP_TYPE_Td;
                }
            } else {
                pg->type = MSYM_POINT_GROUP_TYPE_T;
            }
            
        } else if (primary == NULL){
            if(sigma){
                pg->type = MSYM_POINT_GROUP_TYPE_Cs;
                pg->n = 1;
            } else if(inversion){
                pg->type = MSYM_POINT_GROUP_TYPE_Ci;
                pg->n = 1;
            } else {
                pg->type = MSYM_POINT_GROUP_TYPE_Cn;
                pg->n = 1;
                pg->sops = NULL;
                pg->order = 0;
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
            
            pg->n = primary->order;
            pg->primary = primary;
            
            if(nC2){ //actually nC2 == primary->order but less is acceptable here since we can generate the rest
                if(sigma_h){
                    pg->type = MSYM_POINT_GROUP_TYPE_Dnh;
                } else if (nsigma_v){ //actually nsigma_v == primary->order but less is acceptable here since we can generate the rest
                    pg->type = MSYM_POINT_GROUP_TYPE_Dnd;
                } else {
                    pg->type = MSYM_POINT_GROUP_TYPE_Dn;
                }
                
            } else if (sigma_h) {
                pg->type = MSYM_POINT_GROUP_TYPE_Cnh;
            } else if(nsigma_v) { //actually nsigma_v == primary->order but less is acceptable here since we can generate the rest
                pg->type = MSYM_POINT_GROUP_TYPE_Cnv;
            } else if(s != NULL){
                pg->n = s->order;
                pg->type = MSYM_POINT_GROUP_TYPE_Sn;
            } else {
                pg->type = MSYM_POINT_GROUP_TYPE_Cn;
            }
        }
    } else {
        for(int i = 0; i < sopsl;i++){
            inversion = inversion || sops[i].type == INVERSION;
            
            if(sops[i].type == PROPER_ROTATION && (sops[i].order == 0 || primary == NULL)){                primary = &(sops[i]);
            }
        }
        
        pg->n = primary->order;
        pg->primary = primary;
        
        if(inversion){
            pg->type = MSYM_POINT_GROUP_TYPE_Dnh;
            
        } else {
            pg->type = MSYM_POINT_GROUP_TYPE_Cnv;
        }
    }
    
    return ret;
    
//err:
//    return ret;
    
}


msym_error_t findSubgroup(msym_subgroup_t *subgroup, msym_thresholds_t *thresholds){
    msym_error_t ret = MSYM_SUCCESS;
    int inversion = 0, sigma = 0, nC[6] = {0,0,0,0,0,0}, linear = 0;
    msym_symmetry_operation_t *primary = NULL;
    msym_symmetry_operation_t *s = NULL, *sop = NULL;;
    
    
    for(int i = 0;i < subgroup->order;i++){
        if(subgroup->sops[i]->type == PROPER_ROTATION && subgroup->sops[i]->order == 0){
            linear = 1;
            break;
        } else if (subgroup->sops[i]->type == PROPER_ROTATION && subgroup->sops[i]->order > 2){
            break;
        }
    }
    if(!linear){
        for(int i = 0;i < subgroup->order;i++){
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
                    subgroup->type = MSYM_POINT_GROUP_TYPE_Ih;
                    subgroup->n = primary->order;
                } else {
                    subgroup->type = MSYM_POINT_GROUP_TYPE_I;
                    subgroup->n = primary->order;
                }
            } else if (nC[4] >= 2) {
                if(inversion){
                    subgroup->type = MSYM_POINT_GROUP_TYPE_Oh;
                    subgroup->n = primary->order;
                } else {
                    subgroup->type = MSYM_POINT_GROUP_TYPE_O;
                    subgroup->n = primary->order;
                }
                
            } else if (sigma){
                if(inversion){
                    subgroup->type = MSYM_POINT_GROUP_TYPE_Th;
                    subgroup->n = primary->order;
                } else {
                    subgroup->type = MSYM_POINT_GROUP_TYPE_Td;
                    subgroup->n = primary->order;
                }
            } else {
                subgroup->type = MSYM_POINT_GROUP_TYPE_T;
                subgroup->n = primary->order;
            }
            
        } else if (primary == NULL){
            if(sigma){
                subgroup->type = MSYM_POINT_GROUP_TYPE_Cs;
                subgroup->n = 1;
            } else if(inversion){
                subgroup->type = MSYM_POINT_GROUP_TYPE_Ci;
                subgroup->n = 1;
            } else {
                subgroup->type = MSYM_POINT_GROUP_TYPE_Cn;
                subgroup->n = 1;
            }
        } else {
            int nC2 = 0;
            int sigma_h = 0;
            int nsigma_v = 0;
            
            if(primary->order == 2 && s != NULL && !vparallel(primary->v, s->v, thresholds->angle)){
                for(int i = 0; i < subgroup->order;i++){
                    sop = subgroup->sops[i];
                    if(sop->power > 1) continue;
                    if(sop->type == PROPER_ROTATION && sop->order == 2 && vparallel(sop->v, s->v, thresholds->angle)){
                        primary = sop;
                        break;
                    }
                }
            }
            
            for(int i = 0; i < subgroup->order;i++){
                sop = subgroup->sops[i];
                if(sop->power > 1) continue;
                nC2 += sop->type == PROPER_ROTATION && sop->order == 2 && vperpendicular(primary->v,sop->v, thresholds->angle);
                sigma_h = sigma_h || (sop->type == REFLECTION && vparallel(primary->v,sop->v,thresholds->angle));
                nsigma_v += (sop->type == REFLECTION && vperpendicular(primary->v,sop->v,thresholds->angle));
            }
            if(nC2 == primary->order){
                if(sigma_h){
                    subgroup->type = MSYM_POINT_GROUP_TYPE_Dnh;
                    subgroup->n = primary->order;
                } else if (nsigma_v == primary->order){
                    subgroup->type = MSYM_POINT_GROUP_TYPE_Dnd;
                    subgroup->n = primary->order;
                } else {
                    subgroup->type = MSYM_POINT_GROUP_TYPE_Dn;
                    subgroup->n = primary->order;
                }
                
            } else if (sigma_h) {
                subgroup->type = MSYM_POINT_GROUP_TYPE_Cnh;
                subgroup->n = primary->order;
            } else if(nsigma_v == primary->order) {
                subgroup->type = MSYM_POINT_GROUP_TYPE_Cnv;
                subgroup->n = primary->order;
            } else if(s != NULL){
                subgroup->type = MSYM_POINT_GROUP_TYPE_Sn;
                subgroup->n = s->order;
            } else {
                subgroup->type = MSYM_POINT_GROUP_TYPE_Cn;
                subgroup->n = primary->order;
            }
        }
    } else {
        for(int i = 0; i < subgroup->order;i++){
            inversion = inversion || subgroup->sops[i]->type == INVERSION;
            
            if(subgroup->sops[i]->type == PROPER_ROTATION && (subgroup->sops[i]->order == 0 || primary == NULL)){
                primary = subgroup->sops[i];
            }
        }
        
        if(inversion){
            subgroup->type = MSYM_POINT_GROUP_TYPE_Dnh;
            subgroup->n = primary->order;
        } else {
            subgroup->type = MSYM_POINT_GROUP_TYPE_Cnv;
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
        case (MSYM_POINT_GROUP_TYPE_Ci)  :
        case (MSYM_POINT_GROUP_TYPE_K)   :
        case (MSYM_POINT_GROUP_TYPE_Kh)  :
            break;
        case (MSYM_POINT_GROUP_TYPE_Cs)  :
            for(primary = sops; primary < (primary + sopsl) && primary->type != REFLECTION; primary++){};
        case (MSYM_POINT_GROUP_TYPE_Cn)  :
        case (MSYM_POINT_GROUP_TYPE_Cnh) :
        case (MSYM_POINT_GROUP_TYPE_Sn) :
            if(MSYM_SUCCESS != (ret = reorientAxes(primary, sopsl, sops, thresholds))) goto err;
            if(MSYM_SUCCESS != (ret = transformPrimary(primary, sopsl, sops, thresholds, transform))) goto err;
            break;
        case (MSYM_POINT_GROUP_TYPE_Cnv) :
        case (MSYM_POINT_GROUP_TYPE_Dnh) :
            if(MSYM_SUCCESS != (ret = reorientAxes(primary, sopsl, sops, thresholds))) goto err;
            if(MSYM_SUCCESS != (ret = transformPrimary(primary, sopsl, sops, thresholds, transform))) goto err;
            if(n > 0){
                if(MSYM_SUCCESS != (ret = transformSecondary(type, primary, sopsl, sops, thresholds, transform))) goto err;
            }
            break;
        case (MSYM_POINT_GROUP_TYPE_Dn)  :
        case (MSYM_POINT_GROUP_TYPE_Dnd) :
        case (MSYM_POINT_GROUP_TYPE_O)   :
        case (MSYM_POINT_GROUP_TYPE_Oh)  :
            if(MSYM_SUCCESS != (ret = reorientAxes(primary, sopsl, sops, thresholds))) goto err;
            if(MSYM_SUCCESS != (ret = transformPrimary(primary, sopsl, sops, thresholds, transform))) goto err;
            if(MSYM_SUCCESS != (ret = transformSecondary(type, primary, sopsl, sops, thresholds, transform))) goto err;
            break;
        case (MSYM_POINT_GROUP_TYPE_T)   :
        case (MSYM_POINT_GROUP_TYPE_Td)  :
        case (MSYM_POINT_GROUP_TYPE_Th)  :
        case (MSYM_POINT_GROUP_TYPE_I)   :
        case (MSYM_POINT_GROUP_TYPE_Ih)  : {
            msym_symmetry_operation_t *dprimary = NULL;
            msym_symmetry_operation_t *sop;
            double z = -2.0;
            for(sop = sops; sop < (sops + sopsl); sop++){
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

int isLinearPointGroup(msym_point_group_t *pg){
    return 0 == pg->n && (MSYM_POINT_GROUP_TYPE_Dnh == pg->type || MSYM_POINT_GROUP_TYPE_Cnv == pg->type);
}

int isLinearSubgroup(msym_point_group_t *pg){
    int sub = 0;
    switch(pg->type){
        case MSYM_POINT_GROUP_TYPE_Cnv: sub = pg->n == 0 && pg->order > 2; break;
        case MSYM_POINT_GROUP_TYPE_Dnh: sub = pg->n == 0 && pg->order > 4; break;
        default: break;
    }
    return sub;
}

msym_error_t reduceLinearPointGroup(msym_point_group_t *pg, int n, msym_thresholds_t *thresholds){
    msym_error_t ret = MSYM_SUCCESS;
    int order = 0;
    msym_permutation_t *perm = NULL;
    msym_symmetry_operation_t *sops = NULL;
    msym_symmetry_operation_t *primary = NULL;
    if(!isLinearPointGroup(pg)){
        msymSetErrorDetails("Trying to reduce non linear point group");
        ret = MSYM_POINT_GROUP_ERROR;
        goto err;
    }
    
    if(n == 0) n = 2;
    
    if(MSYM_SUCCESS != (ret = getPointGroupOrder(pg->type, n, &order))) goto err;
    
    if(MSYM_SUCCESS != (ret = generateSymmetryOperations(pg->type, 0, order, &sops))) goto err;
    
    for(int i = 0;i < order;i++){
        if(PROPER_ROTATION == sops[i].type && n == sops[i].order && HORIZONTAL == sops[i].orientation && 1 == sops[i].power){
            primary = &sops[i];
            break;
        }
    }
    
    if(NULL == primary){
        msymSetErrorDetails("Cannot find primary axis when reducing linear group");
        ret = MSYM_POINT_GROUP_ERROR;
        goto err;
    }
    
    
    double T[3][3];
    minv(pg->transform, T);
    
    for(int i = 0; i < order;i++){
        mvmul(sops[i].v,T,sops[i].v);
    }
    
    
    
    if(MSYM_SUCCESS != (ret = findSymmetryOperationPermutations(order,sops, thresholds, &perm))) goto err;
    
    for(int i = 0;i < pg->order && pg->perm != NULL;i++){
        freePermutationData(&pg->perm[i]);
    }
    
    free(pg->sops);
    
    
    pg->sops = sops;
    pg->order = order;
    pg->perm = perm;
    pg->primary = primary;
    
    return ret;
err:
    free(sops);
    return ret;

}

/* Really should split the orientation and class from the symmetry operation
 * and move it to the group so we don't have to regenerate */
msym_error_t pointGroupFromSubgroup(const msym_subgroup_t *sg, msym_thresholds_t *thresholds, msym_point_group_t **opg){
    msym_error_t ret = MSYM_SUCCESS;
    *opg = calloc(1,sizeof(msym_point_group_t));
    msym_point_group_t *pg = *opg;
    pg->type = sg->type;
    pg->n = sg->n;
    pg->sops = malloc(sizeof(msym_symmetry_operation_t[sg->order]));
    memcpy(pg->name,sg->name,sizeof(pg->name));
    
    if(MSYM_SUCCESS != (ret = getPointGroupOrder(pg->type, pg->n, &pg->order))) goto err;
    if(pg->order != sg->order){
        msymSetErrorDetails("Point group order %d does not equal nuber of operations in subgroup %d for point gropu %s",pg->order,sg->order,pg->name);
        ret = MSYM_POINT_GROUP_ERROR;
        goto err;
    }
    
    for(int i = 0;i < sg->order;i++){
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
    
    if(isLinearPointGroup(pg)){
        pg->perm = NULL;
    } else {
        if(MSYM_SUCCESS != (ret = findSymmetryOperationPermutations(pg->order,pg->sops, thresholds, &pg->perm))) goto err;
    }
    
    double T[3][3];
    minv(pg->transform, T);
    
    for(int i = 0; i < pg->order;i++){
        // vector and orientation may not be equal
        if(NULL != sg->primary &&
           NULL == pg->primary &&
           pg->sops[i].type == sg->primary->type &&
           pg->sops[i].order == sg->primary->order &&
           pg->sops[i].power == sg->primary->power){
            
            pg->primary = &pg->sops[i];
        }
        mvmul(pg->sops[i].v,T,pg->sops[i].v);
    }
    
    return ret;
err:
    *opg = NULL;
    free(pg->sops);
    free(pg);
    return ret;
}


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
        case (MSYM_POINT_GROUP_TYPE_Cnv) :
            if(MSYM_SUCCESS != (ret = findSecondaryAxisSigma(primary, sopsl, sops, thresholds, axis))) goto err;
            break;
        case (MSYM_POINT_GROUP_TYPE_O)   :
        case (MSYM_POINT_GROUP_TYPE_Oh)  :
            if(MSYM_SUCCESS != (ret = findSecondaryAxisC4(primary, sopsl, sops, thresholds, axis))) goto err;
            break;
        case (MSYM_POINT_GROUP_TYPE_Dn)  :
        case (MSYM_POINT_GROUP_TYPE_Dnh) :
        case (MSYM_POINT_GROUP_TYPE_Dnd) :
        case (MSYM_POINT_GROUP_TYPE_T)   :
        case (MSYM_POINT_GROUP_TYPE_Td)  :
        case (MSYM_POINT_GROUP_TYPE_Th)  :
            if(MSYM_SUCCESS != (ret = findSecondaryAxisC2(primary, sopsl, sops, thresholds, axis))) goto err;
            break;
        case (MSYM_POINT_GROUP_TYPE_I)   :
        case (MSYM_POINT_GROUP_TYPE_Ih)  :
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
        
        [ 0] = {MSYM_POINT_GROUP_TYPE_Ci,  generateSymmetryOperationsCi},
        [ 1] = {MSYM_POINT_GROUP_TYPE_Cs,  generateSymmetryOperationsCs},
        [ 2] = {MSYM_POINT_GROUP_TYPE_Cn,  generateSymmetryOperationsCn},
        [ 3] = {MSYM_POINT_GROUP_TYPE_Cnh, generateSymmetryOperationsCnh},
        [ 4] = {MSYM_POINT_GROUP_TYPE_Cnv, generateSymmetryOperationsCnv},
        [ 5] = {MSYM_POINT_GROUP_TYPE_Dn,  generateSymmetryOperationsDn},
        [ 6] = {MSYM_POINT_GROUP_TYPE_Dnh, generateSymmetryOperationsDnh},
        [ 7] = {MSYM_POINT_GROUP_TYPE_Dnd, generateSymmetryOperationsDnd},
        [ 8] = {MSYM_POINT_GROUP_TYPE_Sn,  generateSymmetryOperationsSn},
        [ 9] = {MSYM_POINT_GROUP_TYPE_T,   generateSymmetryOperationsT},
        [10] = {MSYM_POINT_GROUP_TYPE_Td,  generateSymmetryOperationsTd},
        [11] = {MSYM_POINT_GROUP_TYPE_Th,  generateSymmetryOperationsTh},
        [12] = {MSYM_POINT_GROUP_TYPE_O,   generateSymmetryOperationsO},
        [13] = {MSYM_POINT_GROUP_TYPE_Oh,  generateSymmetryOperationsOh},
        [14] = {MSYM_POINT_GROUP_TYPE_I,   generateSymmetryOperationsI},
        [15] = {MSYM_POINT_GROUP_TYPE_Ih,  generateSymmetryOperationsIh},
        [16] = {MSYM_POINT_GROUP_TYPE_K,   generateSymmetryOperationsUnknown},
        [17] = {MSYM_POINT_GROUP_TYPE_Kh,  generateSymmetryOperationsUnknown}
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
        clean_debug_printf("i = %d m = %d index = %d ",i,m,index);
        printSymmetryOperation(&sops[index]);
    }
/*
    int ri = k + (((m >> 1)-1) << 1);
    sops[ri].cla = cla + (m >> 1) - 1;
    sops[ri].power = 1;
    sops[ri].p.orientation = HORIZONTAL;
    sops[ri].type = REFLECTION;
    vcopy(z,sops[ri].v);
    clean_debug_printf("replacing symmetry operation %d\n",ri);
    printSymmetryOperation(&sops[ri]);*/
    
    for(int i = 1;i < m >> 1;i++){
        int index = k + 1 + ((i-1) << 1);
        symopPow(&sn, m-i, &sops[index]);
        sops[index].cla = cla + i - 1;
        clean_debug_printf("i = %d m = %d index = %d ",i,m,index);
        printSymmetryOperation(&sops[index]);
        
    }
    k += m - 1;
    cla += m >> 1;
    
    clean_debug_printf("------ Sn %d operations %d classes------\n",k-*pk, cla-*pcla);
    *pk = k; *pcla = cla;
    
    return ret;
err:
    return ret;
}

msym_error_t generateSymmetryOperationsCs(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla;
    
    if(k > l){ret = MSYM_POINT_GROUP_ERROR; msymSetErrorDetails("Too many operations when generating Cs symmetry operations"); goto err;}
    
    
    msym_symmetry_operation_t sigma = {.type = REFLECTION, .orientation = HORIZONTAL, .order = 1, .power = 1, .cla = cla, .v = {0,0,1}};
    
    memcpy(&(sops[k]), &sigma, sizeof(msym_symmetry_operation_t));
    
    k++;
    cla++;
    
    clean_debug_printf("------ Cs %d operations %d classes------\n",k-*pk, cla-*pcla);
    *pk = k; *pcla = cla;
    
    return ret;
err:
    return ret;

}

msym_error_t generateSymmetryOperationsCi(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla;
    
    if(k > l){ret = MSYM_POINT_GROUP_ERROR; msymSetErrorDetails("Too many operations when generating Ci symmetry operations"); goto err;}
    
    
    msym_symmetry_operation_t inv = {.type = INVERSION, .orientation = NONE, .order = 1, .power = 1, .cla = cla, .v = {0,0,0}};
    
    memcpy(&(sops[k]), &inv, sizeof(msym_symmetry_operation_t));
    
    k++;
    cla++;
    
    clean_debug_printf("------ Ci %d operations %d classes------\n",k-*pk, cla-*pcla);
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
        clean_debug_printf("i = %d index = %d ",i,index);
        printSymmetryOperation(&sops[index]);
    }
    
    for(int i = 1;i < (n >> 1) + (n & 1);i++){
        int index = k + (i << 1) - 1;
        symopPow(&cn, n-i, &sops[index]);
        sops[index].cla = cla + ((index - 1) >> 1);
        clean_debug_printf("i = %d index = %d ",i,index);
        printSymmetryOperation(&sops[index]);
    }
    
    k += n - 1;
    cla += n >> 1;
    
    clean_debug_printf("------ Cn %d operations %d classes------\n",k-*pk, cla-*pcla);
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
    for(s = n;s % 2 == 0;s = s >> 1){
        cn.order = s;
        for(int i = 1;i <= s >> 1;i += 2){
            int index = k + ((i >> 1) << 1);
            symopPow(&cn, i, &sops[index]);
            sops[index].cla = cla + (i >> 1);
            printSymmetryOperation(&sops[index]);
        }
        for(int i = 1;i < s >> 1;i += 2){
            int index = k + 1 + ((i >> 1) << 1);
            symopPow(&cn, s-i, &sops[index]);
            sops[index].cla = cla + (i >> 1);
            printSymmetryOperation(&sops[index]);
        }
        
        k += (s >> 1);
        cla += (s >> 2) + ((s >> 1) & 1);
        
        sn.order = s;
        for(int i = 1;i <= s >> 1;i += 2){
            int index = k + ((i >> 1) << 1);
            symopPow(&sn, i, &sops[index]);
            sops[index].cla = cla + (i >> 1);
            printSymmetryOperation(&sops[index]);
        }
        
        for(int i = 1;i < s >> 1;i += 2){
            int index = k + 1 + ((i >> 1) << 1);
            symopPow(&sn, s-i, &sops[index]);
            sops[index].cla = cla + (i >> 1);
            printSymmetryOperation(&sops[index]);
        }
        
        k += (s >> 1);
        cla += (s >> 2) + ((s >> 1) & 1);
        
    }
    
    if(MSYM_SUCCESS != (ret = generateSymmetryOperationsSn(s,l,sops,&k,&cla))) goto err;
    
    clean_debug_printf("------ Cnh %d operations %d classes------\n",k-*pk, cla-*pcla);
    *pk = k; *pcla = cla;
    
    return ret;
err:
    return ret;
    
}

msym_error_t generateSymmetryOperationsCnv(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla;
    
    if(n == 0 && l == 2){ // normal c0v
        if(k + 1 > l){
            ret = MSYM_POINT_GROUP_ERROR;
            msymSetErrorDetails("Too many operations when generating C%dv symmetry operations",n);
            goto err;
        }
        msym_symmetry_operation_t c0 = {.type = PROPER_ROTATION, .order = n, .power = 1, .orientation = HORIZONTAL, .cla = cla, .v = {0,0,1}};
        memcpy(&sops[k], &c0, sizeof(msym_symmetry_operation_t));
        k += 1;
        cla += 1;
        
    } else if(n == 0){ // c0v represented by a subclass cnv where n is even
        int n0 = l/2;
        double z[3] = {0,0,1};
        if(n0 & 1){
            ret = MSYM_POINT_GROUP_ERROR;
            msymSetErrorDetails("Cannot generate an odd representative (C%dv) of C0v",n0);
            goto err;
        }
        
        if(MSYM_SUCCESS != (ret = generateSymmetryOperationsCn(n0,l,sops,&k,&cla))) goto err;
        
        if(k + n0 > l){
            ret = MSYM_POINT_GROUP_ERROR;
            msymSetErrorDetails("Too many operations when generating D0h (D%dh) symmetry operations",n0);
            goto err;
        }
        
        // Can't use generateReflectionPlanes they'll generate vertical and dihedral
        msym_symmetry_operation_t sigma = {.type = REFLECTION, .power = 1, .order = 1, .orientation = VERTICAL, .v = {0,1,0}};
        for(int i = 0;i < n0;i++){
            int index = k+i;
            memcpy(&(sops[index]), &sigma, sizeof(msym_symmetry_operation_t));
            vrotate(i*M_PI/n0, sigma.v, z, sops[index].v);
            sops[index].cla = cla;
        }
        
        k += n0;
        cla += 1;
        
        clean_debug_printf("\n------ C0v %d operations %d classes------\n",k-*pk, cla-*pcla);
    } else {
        if(k + (n << 1) - 1 > l){
            ret = MSYM_POINT_GROUP_ERROR;
            msymSetErrorDetails("Too many operations when generating C%dv symmetry operations",n);
            goto err;
        }
        if(MSYM_SUCCESS != (ret = generateSymmetryOperationsCn(n,l,sops,&k,&cla))) goto err;
        if(MSYM_SUCCESS != (ret = generateReflectionPlanes(n,l,sops,&k,&cla))) goto err;
        clean_debug_printf("------ Cnv %d operations %d classes------\n",k-*pk, cla-*pcla);
        
    }
    
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
    
    clean_debug_printf("\n------ Dn %d operations %d classes------\n",k-*pk, cla-*pcla);
    *pk = k; *pcla = cla;
    
    return ret;
err:
    return ret;
}

msym_error_t generateSymmetryOperationsDnh(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla;
    
    if(n == 0 && l == 4){ // standard d0h
        if(k + 3 > l){
            ret = MSYM_POINT_GROUP_ERROR;
            msymSetErrorDetails("Too many operations when generating D%dh symmetry operations",n);
            goto err;
        }
        if(MSYM_SUCCESS != (ret = generateSymmetryOperationsCnv(n,2,sops,&k,&cla))) goto err;
        msym_symmetry_operation_t sigma = {.type = REFLECTION, .order = 1, .power = 1, .orientation = HORIZONTAL, .cla = cla, .v = {0,0,1}};
        msym_symmetry_operation_t inversion = {.type = INVERSION, .order = 1, .power = 1, .orientation = NONE, .cla = cla+1, .v = {0,0,1}};
        memcpy(&sops[k], &sigma, sizeof(msym_symmetry_operation_t));
        memcpy(&sops[k+1], &inversion, sizeof(msym_symmetry_operation_t));
        k += 2;
        cla += 2;
    }
    else if(n == 0){ // d0h represented by a subclass cnv where n is even
        int n0 = l/4;
        double z[3] = {0,0,1};
        if(n0 & 1){
            ret = MSYM_POINT_GROUP_ERROR;
            msymSetErrorDetails("Cannot generate an odd representative (D%dh) of D0h",n0);
            goto err;
        }
        
        if(MSYM_SUCCESS != (ret = generateSymmetryOperationsCnh(n0,l,sops,&k,&cla))) goto err;
        
        if(k + (n0 << 1) > l){
            ret = MSYM_POINT_GROUP_ERROR;
            msymSetErrorDetails("Too many operations when generating D0h (D%dh) symmetry operations",n0);
            goto err;
        }
        
        // Can't use generateReflectionPlanes/generateC2Axes they'll generate vertical and dihedral
        msym_symmetry_operation_t c2 = {.type = PROPER_ROTATION, .power = 1, .order = 2, .orientation = VERTICAL, .v = {1,0,0}};
        msym_symmetry_operation_t sigma = {.type = REFLECTION, .power = 1, .order = 1, .orientation = VERTICAL, .v = {0,1,0}};
        for(int i = 0;i < n0;i++){
            int index = k+i;
            memcpy(&(sops[index]), &sigma, sizeof(msym_symmetry_operation_t));
            vrotate(i*M_PI/n0, sigma.v, z, sops[index].v);
            sops[index].cla = cla;
            clean_debug_printf("generated sigma at %d\n",index);
            index += n0;
            memcpy(&(sops[index]), &c2, sizeof(msym_symmetry_operation_t));
            vrotate(i*M_PI/n0, c2.v, z, sops[index].v);
            sops[index].cla = cla+1;
        }
        
        k += n0 << 1;
        cla += 2;

        clean_debug_printf("\n------ D0h %d operations %d classes------\n",k-*pk, cla-*pcla);
    } else {
        
        if(k + (n << 2) - 1 > l){
            ret = MSYM_POINT_GROUP_ERROR;
            msymSetErrorDetails("Too many operations when generating D%dh symmetry operations",n);
            goto err;
        }
        if(MSYM_SUCCESS != (ret = generateSymmetryOperationsCnh(n,l,sops,&k,&cla))) goto err;
        if(MSYM_SUCCESS != (ret = generateReflectionPlanes(n,l,sops,&k,&cla))) goto err;
        if(MSYM_SUCCESS != (ret = generateC2Axes(n,l,sops,&k,&cla))) goto err;
        clean_debug_printf("\n------ Dnh %d operations %d classes------\n",k-*pk, cla-*pcla);
        
    }
    
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
    
    
    clean_debug_printf("\n------ Dnd %d operations %d classes------\n",k-*pk, cla-*pcla);
    *pk = k; *pcla = cla;
    
    return ret;
err:
    return ret;
}

msym_error_t generateSymmetryOperationsT(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla;
    
    msym_symmetry_operation_t c2[1] = {
        [0] = {.type = PROPER_ROTATION, .order = 2, .power = 1, .cla = cla, .orientation = HORIZONTAL}
    };
    
    msym_symmetry_operation_t c3[2] = {
        [0] = {.type = PROPER_ROTATION, .order = 3, .power = 1, .cla = cla+1, .orientation = NONE},
        [1] = {.type = PROPER_ROTATION, .order = 3, .power = 2, .cla = cla+1, .orientation = NONE}
    };
    
    if(MSYM_SUCCESS != (ret = generateSymmetryOperationsTetrahedral(l,sops, 1, c2, 0, NULL, 2, c3, &k))) goto err;
    
    cla += 2;
    
    clean_debug_printf("------ T %d operations %d classes------\n",k-*pk, cla-*pcla);
    *pk = k; *pcla = cla;
    
    return ret;
err:
    return ret;
}

msym_error_t generateSymmetryOperationsTd(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla;
    
    msym_symmetry_operation_t c2[3] = {
        [0] = {.type = PROPER_ROTATION, .order = 2, .power = 1, .cla = cla, .orientation = HORIZONTAL},
        [1] = {.type = IMPROPER_ROTATION, .order = 4, .power = 1, .cla = cla+1, .orientation = HORIZONTAL},
        [2] = {.type = IMPROPER_ROTATION, .order = 4, .power = 3, .cla = cla+1, .orientation = HORIZONTAL}
    };
    
    msym_symmetry_operation_t cs[4] = {
        [0] = {.type = REFLECTION, .order = 1, .power = 1, .cla = cla+2, .orientation = DIHEDRAL},
    };
    
    msym_symmetry_operation_t c3[2] = {
        [0] = {.type = PROPER_ROTATION, .order = 3, .power = 1, .cla = cla+3, .orientation = NONE},
        [1] = {.type = PROPER_ROTATION, .order = 3, .power = 2, .cla = cla+3, .orientation = NONE}
    };
    
    if(MSYM_SUCCESS != (ret = generateSymmetryOperationsTetrahedral(l,sops, 3, c2, 1, cs, 2, c3, &k))) goto err;
    
    cla += 4;
    
    clean_debug_printf("------ Td %d operations %d classes------\n",k-*pk, cla-*pcla);
    *pk = k; *pcla = cla;
    
    return ret;
err:
    return ret;
}

msym_error_t generateSymmetryOperationsTh(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla;
    
    msym_symmetry_operation_t c2[2] = {
        [0] = {.type = PROPER_ROTATION, .order = 2, .power = 1, .cla = cla, .orientation = HORIZONTAL},
        [1] = {.type = REFLECTION, .order = 1, .power = 1, .cla = cla+1, .orientation = HORIZONTAL}
    };
    
    msym_symmetry_operation_t c3[4] = {
        [0] = {.type = PROPER_ROTATION, .order = 3, .power = 1, .cla = cla+2, .orientation = NONE},
        [1] = {.type = PROPER_ROTATION, .order = 3, .power = 2, .cla = cla+2, .orientation = NONE},
        [2] = {.type = IMPROPER_ROTATION, .order = 6, .power = 1, .cla = cla+3, .orientation = NONE},
        [3] = {.type = IMPROPER_ROTATION, .order = 6, .power = 5, .cla = cla+3, .orientation = NONE}
    };
    
    if(MSYM_SUCCESS != (ret = generateSymmetryOperationsTetrahedral(l,sops, 2, c2, 0, NULL, 4, c3, &k))) goto err;
    
    if(k - 1 > l){ret = MSYM_POINT_GROUP_ERROR; msymSetErrorDetails("Too many operations when generating Th operations %d >= %d",k, l); goto err;}
    
    sops[k].type = INVERSION;
    sops[k].power = 1;
    sops[k].order = 1;
    sops[k].orientation = NONE;
    sops[k].cla = cla+4;
    k++;
    
    cla += 5;
    
    clean_debug_printf("------ Th %d operations %d classes------\n",k-*pk, cla-*pcla);
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
        [0] = {.type = PROPER_ROTATION, .order = 2, .power = 1, .cla = cla+2, .orientation = HORIZONTAL},
        [1] = {.type = PROPER_ROTATION, .order = 4, .power = 1, .cla = cla+3, .orientation = HORIZONTAL},
        [2] = {.type = PROPER_ROTATION, .order = 4, .power = 3, .cla = cla+3, .orientation = HORIZONTAL}
    };
    
    if(MSYM_SUCCESS != (ret = generateSymmetryOperationsOctahedral(l,sops, 1, c2, 2, c3, 3, c4, &k))) goto err;
    
    cla += 4;
    
    clean_debug_printf("------ O %d operations %d classes------\n",k-*pk, cla-*pcla);
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
        [0] = {.type = PROPER_ROTATION, .order = 2, .power = 1, .cla = cla+4, .orientation = HORIZONTAL},
        [1] = {.type = PROPER_ROTATION, .order = 4, .power = 1, .cla = cla+5, .orientation = HORIZONTAL},
        [2] = {.type = PROPER_ROTATION, .order = 4, .power = 3, .cla = cla+5, .orientation = HORIZONTAL},
        [3] = {.type = IMPROPER_ROTATION, .order = 4, .power = 1, .cla = cla+6, .orientation = HORIZONTAL},
        [4] = {.type = IMPROPER_ROTATION, .order = 4, .power = 3, .cla = cla+6, .orientation = HORIZONTAL},
        [5] = {.type = REFLECTION, .order = 1, .power = 1, .cla = cla+7, .orientation = HORIZONTAL}
    };
    
    if(MSYM_SUCCESS != (ret = generateSymmetryOperationsOctahedral(l,sops, 2, c2, 4, c3, 6, c4, &k))) goto err;
    
    if(k - 1 > l){ret = MSYM_POINT_GROUP_ERROR; msymSetErrorDetails("Too many operations when generating Oh operations %d >= %d",k, l); goto err;}
    
    sops[k].type = INVERSION;
    sops[k].power = 1;
    sops[k].order = 1;
    sops[k].orientation = NONE;
    sops[k].cla = cla+8;
    k++;
    
    cla += 9;
    
    clean_debug_printf("------ O %d operations %d classes------\n",k-*pk, cla-*pcla);
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
    
    clean_debug_printf("------ I %d operations %d classes------\n",k-*pk, cla-*pcla);
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
    
    if(k - 1 > l){ret = MSYM_POINT_GROUP_ERROR; msymSetErrorDetails("Too many operations when generating Ih operations %d >= %d",k + 12, l); goto err;}
    
    sops[k].type = INVERSION;
    sops[k].power = 1;
    sops[k].order = 1;
    sops[k].orientation = NONE;
    sops[k].cla = cla+8;
    k++;
    
    cla += 9;
    
    clean_debug_printf("------ Ih %d operations %d classes------\n",k-*pk, cla-*pcla);
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
    
    clean_debug_printf("------ Tetra %d operations 0 classes------\n",k-*pk);
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
    
    clean_debug_printf("------ Octa %d operations 0 classes------\n",k-*pk);
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
    
    clean_debug_printf("------ Icosa %d operations 0 classes------\n",k-*pk);
    *pk = k;
    
    return ret;
err:
    return ret;
    
}

msym_error_t generateReflectionPlanes(int n, int l, msym_symmetry_operation_t sops[l], int *pk, int *pcla){
    msym_error_t ret = MSYM_SUCCESS;
    int k = *pk, cla = *pcla;
    double z[3] = {0.0,0.0,1.0}, y[3] = {0.0,1.0,0.0};
    msym_symmetry_operation_t sigma = {.type = REFLECTION, .power = 1, .order = 1};
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
    clean_debug_printf("------ R %d operations %d classes------\n",k-*pk, cla-*pcla);
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
    clean_debug_printf("------ C2 %d operations %d classes------\n",k-*pk, cla-*pcla);
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
    
    if(isLinearSubgroup(pg)){
        switch (pg->type) {
            case MSYM_POINT_GROUP_TYPE_Cnv : n = pg->order / 4; break;
            case MSYM_POINT_GROUP_TYPE_Dnh : n = pg->order / 2; break;
            default: break;
        }
    }
    
    switch (pg->type) {
        case MSYM_POINT_GROUP_TYPE_Kh  : size = -1; break;
        case MSYM_POINT_GROUP_TYPE_K   : size = -1; break;
        case MSYM_POINT_GROUP_TYPE_Ih  : size = 162; break;
        case MSYM_POINT_GROUP_TYPE_I   : size = 57; break;
        case MSYM_POINT_GROUP_TYPE_Oh  : size = 96; break;
        case MSYM_POINT_GROUP_TYPE_O   : size = 28; break;
        case MSYM_POINT_GROUP_TYPE_Th  : size = 24; break;
        case MSYM_POINT_GROUP_TYPE_Td  : size = 28; break;
        case MSYM_POINT_GROUP_TYPE_T   : size = 9; break;
        case MSYM_POINT_GROUP_TYPE_Ci  : size = 0; break;
        case MSYM_POINT_GROUP_TYPE_Cs  : size = 0; break;
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
                case MSYM_POINT_GROUP_TYPE_Cnh : {
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
                case MSYM_POINT_GROUP_TYPE_Dn  :
                case MSYM_POINT_GROUP_TYPE_Cnv : size = n + ndiv + sdiv; break;
                case MSYM_POINT_GROUP_TYPE_Cn  : size = ndiv - 1; break;
                case MSYM_POINT_GROUP_TYPE_Dnh : {
                    if(n % 2 == 0) size = 4*n + 2*ndiv + 3*sdiv + 4 + neven + seven;
                    else size = 3*(n+sdiv+1) + 2*ndiv;
                    break;
                }
                case MSYM_POINT_GROUP_TYPE_Dnd : {
                    if(n % 2 == 0) size = 2*n + 3 + ndiv + 2*sdiv + nodd + sodd;
                    else size = 3*(n+sdiv+1) + 2*ndiv;
                    break;
                }
                case MSYM_POINT_GROUP_TYPE_Sn : size = ndiv - 1; break;
                default : break;
            }
        }
    }
    
    return size;
}

