//
//  character_table.c
//  libmsym
//
//  Created by Marcus Johansson on 28/11/14.
//  Copyright (c) 2014 Marcus Johansson. 
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "character_table.h"
#include "point_group.h"

#include "debug.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288419716939937510582
#endif

typedef struct _msym_representation {
    enum {IRREDUCIBLE = 1, REDUCIBLE = 2} type;
    int d;
    struct {
        int p, v, h, i, l;
    } eig;
    char name[8];
} msym_representation_t;

msym_error_t verifyCharacterTable(msym_character_table_t *ct);

msym_error_t setRepresentationName(msym_representation_t *rep);
msym_error_t representationCharacter(int n, msym_symmetry_operation_t *sop, msym_representation_t *rep, double *c);

msym_error_t getRepresentationsCi(int n, int rl, msym_representation_t rep[rl]);
msym_error_t getRepresentationsCs(int n, int rl, msym_representation_t rep[rl]);
msym_error_t getRepresentationsCn(int n, int rl, msym_representation_t rep[rl]);
msym_error_t getRepresentationsCnh(int n, int rl, msym_representation_t rep[rl]);
msym_error_t getRepresentationsCnv(int n, int rl, msym_representation_t rep[rl]);
msym_error_t getRepresentationsSn(int n, int rl, msym_representation_t rep[rl]);
msym_error_t getRepresentationsDn(int n, int rl, msym_representation_t rep[rl]);
msym_error_t getRepresentationsDnh(int n, int rl, msym_representation_t rep[rl]);
msym_error_t getRepresentationsDnd(int n, int rl, msym_representation_t rep[rl]);

msym_error_t getCharacterTableT(int sopsl, msym_symmetry_operation_t sops[sopsl], msym_character_table_t *ct);
msym_error_t getCharacterTableTd(int sopsl, msym_symmetry_operation_t sops[sopsl], msym_character_table_t *ct);
msym_error_t getCharacterTableTh(int sopsl, msym_symmetry_operation_t sops[sopsl], msym_character_table_t *ct);
msym_error_t getCharacterTableO(int sopsl, msym_symmetry_operation_t sops[sopsl], msym_character_table_t *ct);
msym_error_t getCharacterTableOh(int sopsl, msym_symmetry_operation_t sops[sopsl], msym_character_table_t *ct);
msym_error_t getCharacterTableI(int sopsl, msym_symmetry_operation_t sops[sopsl], msym_character_table_t *ct);
msym_error_t getCharacterTableIh(int sopsl, msym_symmetry_operation_t sops[sopsl], msym_character_table_t *ct);
msym_error_t getCharacterTableUnknown(int sopsl, msym_symmetry_operation_t sops[sopsl], msym_character_table_t *ct);


msym_error_t getPredefinedCharacterTable(int sopsl, msym_symmetry_operation_t sops[sopsl], int l, const msym_symmetry_operation_t tsops[l], const char *tname[l], const int tdim[l], const int tred[l], const double (*table)[l], msym_character_table_t *ct);

msym_error_t getRepresentationName(msym_point_group_type_t type, int n, msym_representation_t *rep, int l, char name[l]);


void decomposeRepresentation(msym_character_table_t *ct, double rspan[ct->d], double dspan[ct->d]){
    int order = 0;
    double (*ctable)[ct->d] = ct->table;
    memset(dspan,0, sizeof(double[ct->d]));
    
    for(int k = 0;k < ct->d;k++){
        order += ct->classc[k];
        for(int j = 0; j < ct->d;j++) dspan[k] += ct->classc[j]*rspan[j]*ctable[k][j];
    }
    for(int k = 0;k < ct->d;k++) dspan[k] /= order;
}


void directProduct(int l, double irrep1[l], double irrep2[l], double pspan[l]){
    for(int i = 0;i < l;i++) pspan[i] = irrep1[i]*irrep2[i];
}

msym_error_t generateCharacterTable(msym_point_group_type_t type, int n, int sopsl, msym_symmetry_operation_t sops[sopsl], msym_character_table_t **oct){
    msym_error_t ret = MSYM_SUCCESS;
    int d = sops[sopsl-1].cla + 1; //max cla
    int linear = n == 0 && ((MSYM_POINT_GROUP_TYPE_Dnh == type && sopsl == 4) || (MSYM_POINT_GROUP_TYPE_Cnv == type && sopsl == 2));
    
    msym_character_table_t *ct = calloc(1, sizeof(msym_character_table_t) + sizeof(int[d]) + sizeof(msym_symmetry_species_t[d]) + sizeof(msym_symmetry_operation_t*[d]) + sizeof(double[d][d]));
    
    
    ct->table = (double (*)[d])(ct + 1);
    ct->s = (msym_symmetry_species_t*)((double (*)[d])ct->table + d);
    ct->sops = (msym_symmetry_operation_t **)(ct->s + d);
    ct->classc = (int *)(ct->sops + d);
    
    ct->d = d;

    msym_representation_t *rep = calloc(ct->d, sizeof(msym_representation_t));;
    double (*table)[ct->d] = (double (*)[ct->d]) ct->table;
    
    
    const struct _fmap {
        msym_point_group_type_t type;
        enum {REP, TAB} c;
        union {
            msym_error_t (*fr)(int, int, msym_representation_t *);
            msym_error_t (*ft)(int, msym_symmetry_operation_t *, msym_character_table_t *);
        } f;
    } fmap[18] = {
        [ 0] = {MSYM_POINT_GROUP_TYPE_Ci,  REP, .f.fr = getRepresentationsCi},
        [ 1] = {MSYM_POINT_GROUP_TYPE_Cs,  REP, .f.fr = getRepresentationsCs},
        [ 2] = {MSYM_POINT_GROUP_TYPE_Cn,  REP, .f.fr = getRepresentationsCn},
        [ 3] = {MSYM_POINT_GROUP_TYPE_Cnh, REP, .f.fr = getRepresentationsCnh},
        [ 4] = {MSYM_POINT_GROUP_TYPE_Cnv, REP, .f.fr = getRepresentationsCnv},
        [ 5] = {MSYM_POINT_GROUP_TYPE_Dn,  REP, .f.fr = getRepresentationsDn},
        [ 6] = {MSYM_POINT_GROUP_TYPE_Dnh, REP, .f.fr = getRepresentationsDnh},
        [ 7] = {MSYM_POINT_GROUP_TYPE_Dnd, REP, .f.fr = getRepresentationsDnd},
        [ 8] = {MSYM_POINT_GROUP_TYPE_Sn,  REP, .f.fr = getRepresentationsSn},
        [ 9] = {MSYM_POINT_GROUP_TYPE_T,   TAB, .f.ft = getCharacterTableT},
        [10] = {MSYM_POINT_GROUP_TYPE_Td,  TAB, .f.ft = getCharacterTableTd},
        [11] = {MSYM_POINT_GROUP_TYPE_Th,  TAB, .f.ft = getCharacterTableTh},
        [12] = {MSYM_POINT_GROUP_TYPE_O,   TAB, .f.ft = getCharacterTableO},
        [13] = {MSYM_POINT_GROUP_TYPE_Oh,  TAB, .f.ft = getCharacterTableOh},
        [14] = {MSYM_POINT_GROUP_TYPE_I,   TAB, .f.ft = getCharacterTableI},
        [15] = {MSYM_POINT_GROUP_TYPE_Ih,  TAB, .f.ft = getCharacterTableIh},
        [16] = {MSYM_POINT_GROUP_TYPE_K,   TAB, .f.ft = getCharacterTableUnknown},
        [17] = {MSYM_POINT_GROUP_TYPE_Kh,  TAB, .f.ft = getCharacterTableUnknown}
        
    };
    
    int fi, fil = sizeof(fmap)/sizeof(fmap[0]);
    for(fi = 0; fi < fil;fi++){
        if(fmap[fi].type == type) {
            if(MSYM_SUCCESS != (ret = fmap[fi].c == REP ? fmap[fi].f.fr(n,ct->d,rep) : fmap[fi].f.ft(sopsl,sops,ct))) goto err;
            break;
        }
    }
    
    if(fi == fil){
        msymSetErrorDetails("Unknown point group when generating character table");
        ret = MSYM_POINT_GROUP_ERROR;
        goto err;
    }
    
    for(int i = 0; i < sopsl;i++){
        ct->classc[sops[i].cla]++;
    }
    
    if(fmap[fi].c == REP){
        for(int i = 0;i < ct->d;i++){
            if(MSYM_SUCCESS != (ret = getRepresentationName(type, n, &rep[i], sizeof(ct->s[i].name), ct->s[i].name))) goto err;
            ct->s[i].d = rep[i].d;
            ct->s[i].r = rep[i].type;
            int nc = -1;
            for(int j = 0;j < sopsl;j++){
                if(nc < sops[j].cla){
                    nc = sops[j].cla;
                    if(MSYM_SUCCESS != (ret = representationCharacter(n,&sops[j],&rep[i],&table[i][nc]))) goto err;
                }
            }
        }
    }
    
    
    for (int i = 0; i < ct->d; i++) {
        for(int j = 0;j < sopsl;j++){
            if(sops[j].cla == i){
                ct->sops[i] = &sops[j];
                break;
            }
        }
    }

    debug_printCharacterTable(ct);
    
    if(!linear && MSYM_SUCCESS != (ret = verifyCharacterTable(ct))) goto err;
    
    *oct = ct;
    
    free(rep);
    return ret;
err:
    free(rep);
    free(ct);
    return ret;
}

#define CHARACTER_TABLE_VERIFICATION_THRESHOLD 1e-10
msym_error_t verifyCharacterTable(msym_character_table_t *ct){
    msym_error_t ret = MSYM_SUCCESS;
    double (*table)[ct->d] = (double (*)[ct->d]) ct->table;
    for(int i = 0;i < ct->d && ret == MSYM_SUCCESS;i++){
        for(int j = i+1;j < ct->d;j++){
            double r = 0.0;
            for(int k = 0;k < ct->d;k++){
                r += ct->classc[k]*table[i][k]*table[j][k];
            }
            if(r > CHARACTER_TABLE_VERIFICATION_THRESHOLD){
                msymSetErrorDetails("Character table verification failed irrep %s(%d) and %s(%d) are not orthogonal, product %e > %e",ct->s[i].name,i,ct->s[j].name,j,r,CHARACTER_TABLE_VERIFICATION_THRESHOLD);
                ret = MSYM_INVALID_CHARACTER_TABLE;
            }
        }
    }
    return ret;
}

msym_error_t getRepresentationsCi(int n, int rl, msym_representation_t rep[rl]){
    msym_error_t ret = MSYM_SUCCESS;
    int r = 0;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
    r++;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.v = rep[r].eig.h = 1;
    rep[r].eig.i = -1;
    r++;
    
    return ret;
//err:
//    return ret;
}

msym_error_t getRepresentationsCs(int n, int rl, msym_representation_t rep[rl]){
    msym_error_t ret = MSYM_SUCCESS;
    int r = 0;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
    r++;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.v = rep[r].eig.i = 1;
    rep[r].eig.h = -1;
    r++;
    
    return ret;
//err:
//    return ret;
}

msym_error_t getRepresentationsCn(int n, int rl, msym_representation_t rep[rl]){
    msym_error_t ret = MSYM_SUCCESS;
    int r = 0;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
    r++;
    if(!(n & 1)){
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
        rep[r].eig.p = -1;
        r++;
    }
    for(int i = 1;r < rl;i++, r++){
        rep[r].type = REDUCIBLE;
        rep[r].d = 2;
        rep[r].eig.l = i;
        rep[r].eig.p = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
    }
    
    return ret;
//err:
//    return ret;
}

msym_error_t getRepresentationsCnh(int n, int rl, msym_representation_t rep[rl]){
    msym_error_t ret = MSYM_SUCCESS;
    int r = 0;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
    r++;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.v = 1;
    rep[r].eig.h = rep[r].eig.i = -1;
    r++;
    if(!(n & 1)){
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.v = rep[r].eig.i = 1;
        rep[r].eig.p = -1;
        rep[r].eig.h = 1 - (n & 2);
        r++;
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.v = 1;
        rep[r].eig.p = rep[r].eig.i = -1;
        rep[r].eig.h = -1 + (n & 2);
        r++;
    }
    for(int i = 1;r < rl;i++){
        rep[r].type = REDUCIBLE;
        rep[r].d = 2;
        rep[r].eig.l = i;
        rep[r].eig.p = rep[r].eig.v = rep[r].eig.h = 1;
        rep[r].eig.i = 1 - ((i & 1) << 1);
        r++;
        rep[r].type = REDUCIBLE;
        rep[r].d = 2;
        rep[r].eig.l = i;
        rep[r].eig.p = rep[r].eig.v = 1;
        rep[r].eig.h = -1;
        rep[r].eig.i = -1 + ((i & 1) << 1);
        r++;
    }
    
    return ret;
//err:
//    return ret;
}

msym_error_t getRepresentationsCnv(int n, int rl, msym_representation_t rep[rl]){
    msym_error_t ret = MSYM_SUCCESS;
    
    int r = 0;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
    r++;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.h = rep[r].eig.i = 1;
    rep[r].eig.v = -1;
    r++;
    if(n && !(n & 1)){
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.v = rep[r].eig.i = rep[r].eig.h = 1;
        rep[r].eig.p = -1;
        r++;
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.h = rep[r].eig.i = 1 ;
        rep[r].eig.p = rep[r].eig.v = -1;
        r++;
    }
    for(int i = 1;r < rl;i++, r++){
        rep[r].type = (n == 0 && (r == rl - 1)) ? REDUCIBLE : IRREDUCIBLE;
        rep[r].d = 2;
        rep[r].eig.l = i;
        rep[r].eig.p = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
    }
    
    return ret;
//err:
//    return ret;
}

msym_error_t getRepresentationsSn(int n, int rl, msym_representation_t rep[rl]){
    msym_error_t ret = MSYM_SUCCESS;
    int r = 0;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
    r++;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.v = 1;
    rep[r].eig.h = rep[r].eig.i = -1;
    r++;
    
    if(!((n >> 1) & 1)){
        for(int i = 1;r < rl;i++, r++){
            rep[r].type = REDUCIBLE;
            rep[r].d = 2;
            rep[r].eig.l = i;
            rep[r].eig.p = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
        }
    } else {
        for(int i = 1;r < rl;i++){
            rep[r].type = REDUCIBLE;
            rep[r].d = 2;
            rep[r].eig.l = i;
            rep[r].eig.p = rep[r].eig.v = rep[r].eig.i = 1;
            rep[r].eig.h = 1 - ((i % 2) << 1);
            r++;
            rep[r].type = REDUCIBLE;
            rep[r].d = 2;
            rep[r].eig.l = i;
            rep[r].eig.p = rep[r].eig.v = 1;
            rep[r].eig.h = -1 + ((i % 2) << 1);
            rep[r].eig.i = -1;
            r++;
        }
    }
    
    return ret;
//err:
//    return ret;
}

msym_error_t getRepresentationsDn(int n, int rl, msym_representation_t rep[rl]){
    msym_error_t ret = MSYM_SUCCESS;
    int r = 0;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
    r++;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.h = rep[r].eig.i = 1;
    rep[r].eig.v = -1;
    r++;
    if(~n & 1){
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.v = rep[r].eig.i = rep[r].eig.h = 1;
        rep[r].eig.p = -1;
        r++;
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.h = rep[r].eig.i = 1 ;
        rep[r].eig.p = rep[r].eig.v = -1;
        r++;
    }
    for(int i = 1;r < rl;i++, r++){
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 2;
        rep[r].eig.l = i;
        rep[r].eig.p = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
    }
    
    return ret;
//err:
//    return ret;
}

msym_error_t getRepresentationsDnh(int n, int rl, msym_representation_t rep[rl]){
    msym_error_t ret = MSYM_SUCCESS;
    
    int r = 0;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
    r++;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.h = rep[r].eig.i = 1;
    rep[r].eig.v = -1;
    r++;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.v = 1;
    rep[r].eig.h = rep[r].eig.i = -1;
    r++;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = 1;
    rep[r].eig.h = rep[r].eig.i = rep[r].eig.v = -1;
    r++;
    if(n && !(n & 1)){
        
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.v = rep[r].eig.i = 1;
        rep[r].eig.p = -1;
        rep[r].eig.h = 1 - (n & 2);
        r++;
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.v = 1;
        rep[r].eig.p = rep[r].eig.i = -1;
        rep[r].eig.h = -1 + (n & 2);
        r++;
        
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.i = 1;
        rep[r].eig.p = rep[r].eig.v = -1;
        rep[r].eig.h = 1 - (n & 2);
        r++;
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = 1;
        rep[r].eig.p = rep[r].eig.v = rep[r].eig.i = -1;
        rep[r].eig.h = -1 + (n & 2);
        r++;
    }
    for(int i = 1;r < rl;i++){
        rep[r].type = (n == 0 && (r == rl - 2)) ? REDUCIBLE : IRREDUCIBLE;
        rep[r].d = 2;
        rep[r].eig.l = i;
        rep[r].eig.p = rep[r].eig.v = rep[r].eig.h = 1;
        rep[r].eig.i = 1 - ((i & 1) << 1);
        r++;
        rep[r].type = (n == 0 && (r == rl - 1)) ? REDUCIBLE : IRREDUCIBLE;
        rep[r].d = 2;
        rep[r].eig.l = i;
        rep[r].eig.p = rep[r].eig.v = 1;
        rep[r].eig.h = -1;
        rep[r].eig.i = -1 + ((i & 1) << 1);
        r++;
    }
    
    return ret;
//err:
//    return ret;
}

msym_error_t getRepresentationsDnd(int n, int rl, msym_representation_t rep[rl]){
    msym_error_t ret = MSYM_SUCCESS;
    int r = 0;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
    r++;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.h = rep[r].eig.i = 1;
    rep[r].eig.v = -1;
    r++;
    if(!(n & 1)){
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.v = rep[r].eig.i = rep[r].eig.p = 1;
        rep[r].eig.h = -1;
        r++;
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.i = rep[r].eig.p = 1 ;
        rep[r].eig.h = rep[r].eig.v = -1;
        r++;
        for(int i = 1;r < rl;i++, r++){
            rep[r].type = IRREDUCIBLE;
            rep[r].d = 2;
            rep[r].eig.l = i;
            rep[r].eig.p = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
        }
    } else {
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.v = rep[r].eig.p = 1;
        rep[r].eig.h = rep[r].eig.i = -1;
        r++;
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.p = 1 ;
        rep[r].eig.h = rep[r].eig.i = rep[r].eig.v = -1;
        r++;
        for(int i = 1;r < rl;i++){
            rep[r].type = IRREDUCIBLE;
            rep[r].d = 2;
            rep[r].eig.l = i;
            rep[r].eig.p = rep[r].eig.v = rep[r].eig.i = 1;
            rep[r].eig.h = 1 - ((i % 2) << 1);
            r++;
            rep[r].type = IRREDUCIBLE;
            rep[r].d = 2;
            rep[r].eig.l = i;
            rep[r].eig.p = rep[r].eig.v = 1;
            rep[r].eig.h = -1 + ((i % 2) << 1);
            rep[r].eig.i = -1;
            r++;
        }
        
    }
    
    
    return ret;
//err:
//    return ret;
}


msym_error_t representationCharacter(int n, msym_symmetry_operation_t *sop, msym_representation_t *rep, double *c){
    msym_error_t ret = MSYM_SUCCESS;
    double x = 0;
    if(sop->orientation == HORIZONTAL){
        switch(rep->d)
        {
            case 1: {
                switch (sop->type) {
                    case IDENTITY           : x = 1; break;
                    case REFLECTION         : x = rep->eig.h; break;
                    case INVERSION          : x = rep->eig.i; break;
                    case PROPER_ROTATION    : x = !n || ((n/sop->order) & 1) ? rep->eig.p : 1 ; break;
                    case IMPROPER_ROTATION  : x = rep->eig.h*(!n || ((n/sop->order) & 1) ? rep->eig.p : 1); break;
                    default :
                        ret = MSYM_INVALID_CHARACTER_TABLE;
                        msymSetErrorDetails("Invalid symmetry operation when building character table");
                        goto err;
                        
                }
                break;
            }
            case 2 : {
                switch (sop->type) {
                    case IDENTITY           : x = 2; break;
                    case REFLECTION         : x = 2*rep->eig.h; break;
                    case INVERSION          : x = 2*rep->eig.i; break;
                    case PROPER_ROTATION    : x = 2*cos(2*rep->eig.l*sop->power*(M_PI/sop->order)); break;
                    case IMPROPER_ROTATION  : x = rep->eig.h*2*cos(2*rep->eig.l*sop->power*(M_PI/sop->order)); break;
                    default :
                        ret = MSYM_INVALID_CHARACTER_TABLE;
                        msymSetErrorDetails("Invalid symmetry operation when building character table");
                        goto err;
                        
                }
                break;
            }
            default :
                ret = MSYM_INVALID_CHARACTER_TABLE;
                msymSetErrorDetails("Invalid dimension (%d) of irreducible representation for point group",rep->d);
                goto err;
        }
    } else {
        switch(rep->d)
        {
            case 1: {
                switch (sop->type) {
                    case IDENTITY           : x = 1; break;
                    case INVERSION          : x = rep->eig.i; break;
                    case REFLECTION         : x = sop->orientation == VERTICAL ? rep->eig.v*rep->eig.h : rep->eig.p*rep->eig.v*rep->eig.h; break;
                    case PROPER_ROTATION    : x = sop->orientation == VERTICAL ? rep->eig.v : rep->eig.p*rep->eig.v; break;
                    case IMPROPER_ROTATION  :
                    default :
                        ret = MSYM_INVALID_CHARACTER_TABLE;
                        msymSetErrorDetails("Invalid symmetry operation when building character table");
                        goto err;
                }
                break;
            }
            case 2 : {
                switch (sop->type) {
                    case IDENTITY           : x = 2; break;
                    case REFLECTION         : x = 0; break;
                    case INVERSION          : x = 2*rep->eig.i; break;
                    case PROPER_ROTATION    : x = 0; break;
                    case IMPROPER_ROTATION  :
                    default :
                        ret = MSYM_INVALID_CHARACTER_TABLE;
                        msymSetErrorDetails("Invalid symmetry operation when building character table");
                        goto err;
                        
                }
                break;
            }
            default :
                ret = MSYM_INVALID_CHARACTER_TABLE;
                msymSetErrorDetails("Invalid dimension (%d) of irreducible representation for point group",rep->d);
                goto err;
        }
        
    }
    
    *c = x;
err:
    return ret;
}

msym_error_t getCharacterTableUnknown(int sopsl, msym_symmetry_operation_t sops[sopsl], msym_character_table_t *ct){
    msym_error_t ret = MSYM_SUCCESS;
    msymSetErrorDetails("Character table NYI");
    ret = MSYM_INVALID_CHARACTER_TABLE;
    
//err:
    return ret;
}

msym_error_t getCharacterTableT(int sopsl, msym_symmetry_operation_t sops[sopsl], msym_character_table_t *ct){
    msym_error_t ret = MSYM_SUCCESS;
    
    const msym_symmetry_operation_t tsops[3] = {
        [0] = {.type = IDENTITY, .order = 1, .power = 1, .orientation = NONE},
        [1] = {.type = PROPER_ROTATION, .order = 3, .power = 1, .orientation = NONE},
        [2] = {.type = PROPER_ROTATION, .order = 2, .power = 1, .orientation = HORIZONTAL}
    };
    
    const char *tname[3] = {"A","E","T"};
    const int tdim[3] = {1,2,3};
    const int tred[3] = {1,2,1};
    
    const double table[][3] = {
        [0] = {1,  1,  1},
        [1] = {2, -1,  2}, // Reducible to {1 e e* 1},{1 e* e 1} where e = e^(i2pi/3)
        [2] = {3,  0, -1}
    };
    
    if(MSYM_SUCCESS != (ret = getPredefinedCharacterTable(sopsl, sops, 3, tsops, tname, tdim, tred, table, ct))) goto err;
    
err:
    return ret;
    
}

msym_error_t getCharacterTableTd(int sopsl, msym_symmetry_operation_t sops[sopsl], msym_character_table_t *ct){
    msym_error_t ret = MSYM_SUCCESS;
    
    const msym_symmetry_operation_t tsops[5] = {
        [0] = {.type = IDENTITY, .order = 1, .power = 1, .orientation = NONE},
        [1] = {.type = PROPER_ROTATION, .order = 2, .power = 1, .orientation = HORIZONTAL},
        [2] = {.type = PROPER_ROTATION, .order = 3, .power = 1, .orientation = NONE},
        [3] = {.type = IMPROPER_ROTATION, .order = 4, .power = 1, .orientation = HORIZONTAL},
        [4] = {.type = REFLECTION, .order = 1, .power = 1, .orientation = DIHEDRAL},
    };
    
    const char *tname[5] = {"A1","A2","E","T1","T2"};
    const int tdim[5] = {1,1,2,3,3};
    const int tred[5] = {1,1,1,1,1};
    
    const double table[][5] = {
        [0] = {1,  1,  1,  1,  1},
        [1] = {1,  1,  1, -1, -1},
        [2] = {2,  2, -1,  0,  0},
        [3] = {3, -1,  0,  1, -1},
        [4] = {3, -1,  0, -1,  1}
    };
    
    if(MSYM_SUCCESS != (ret = getPredefinedCharacterTable(sopsl, sops, 5, tsops, tname, tdim, tred, table, ct))) goto err;
    
err:
    return ret;
}

msym_error_t getCharacterTableTh(int sopsl, msym_symmetry_operation_t sops[sopsl], msym_character_table_t *ct){
    msym_error_t ret = MSYM_SUCCESS;
    
    const msym_symmetry_operation_t tsops[6] = {
        
        [0] = {.type = IDENTITY, .order = 1, .power = 1, .orientation = NONE},
        [1] = {.type = PROPER_ROTATION, .order = 3, .power = 1, .orientation = NONE},
        [2] = {.type = PROPER_ROTATION, .order = 2, .power = 1, .orientation = HORIZONTAL},
        [3] = {.type = INVERSION, .order = 1, .power = 1, .orientation = NONE},
        [4] = {.type = IMPROPER_ROTATION, .order = 6, .power = 1, .orientation = NONE},
        [5] = {.type = REFLECTION, .order = 1, .power = 1, .orientation = HORIZONTAL}
    };
    
    
    const char *tname[6] = {"Ag","Au","Eg","Eu","Tg","Tu"};
    const int tdim[6] = {1,1,2,2,3,3};
    const int tred[6] = {1,1,2,2,1,1};
    
    const double table[][6] = {
        [0] = {1,  1,  1,  1,  1,  1},
        [1] = {1,  1,  1, -1, -1, -1},
        [2] = {2, -1,  2,  2, -1,  2},
        [3] = {2, -1,  2, -2,  1, -2},
        [4] = {3,  0, -1,  3,  0, -1},
        [5] = {3,  0, -1, -3,  0,  1}
    };
    
    if(MSYM_SUCCESS != (ret = getPredefinedCharacterTable(sopsl, sops, 6, tsops, tname, tdim, tred, table, ct))) goto err;
    
err:
    return ret;
}


msym_error_t getCharacterTableO(int sopsl, msym_symmetry_operation_t sops[sopsl], msym_character_table_t *ct){
    msym_error_t ret = MSYM_SUCCESS;
    const msym_symmetry_operation_t tsops[5] = {
        [0] = {.type = IDENTITY, .order = 1, .power = 1, .orientation = NONE},
        [1] = {.type = PROPER_ROTATION, .order = 3, .power = 1, .orientation = NONE},
        [2] = {.type = PROPER_ROTATION, .order = 2, .power = 1, .orientation = VERTICAL},
        [3] = {.type = PROPER_ROTATION, .order = 4, .power = 1, .orientation = HORIZONTAL},
        [4] = {.type = PROPER_ROTATION, .order = 2, .power = 1, .orientation = HORIZONTAL}
    };
    
    const double table[][5] = {
        [0] = {1,  1,  1,  1,  1},
        [1] = {1,  1, -1, -1,  1},
        [2] = {2, -1,  0,  0,  2},
        [3] = {3,  0, -1,  1, -1},
        [4] = {3,  0,  1, -1, -1}
    };
    
    const char *tname[5] = {"A1","A2","E","T1","T2"};
    const int tdim[5] = {1,1,2,3,3};
    const int tred[5] = {1,1,1,1,1};
    
    
    if(MSYM_SUCCESS != (ret = getPredefinedCharacterTable(sopsl, sops, 5, tsops, tname, tdim, tred, table, ct))) goto err;
    
err:
    return ret;
}

msym_error_t getCharacterTableOh(int sopsl, msym_symmetry_operation_t sops[sopsl], msym_character_table_t *ct){
    msym_error_t ret = MSYM_SUCCESS;
    const msym_symmetry_operation_t tsops[10] = {
        [0] = {.type = IDENTITY, .order = 1, .power = 1, .orientation = NONE},
        [1] = {.type = PROPER_ROTATION, .order = 4, .power = 1, .orientation = HORIZONTAL},
        [2] = {.type = PROPER_ROTATION, .order = 2, .power = 1, .orientation = HORIZONTAL},
        [3] = {.type = PROPER_ROTATION, .order = 3, .power = 1, .orientation = NONE},
        [4] = {.type = PROPER_ROTATION, .order = 2, .power = 1, .orientation = VERTICAL},
        [5] = {.type = INVERSION, .order = 1, .power = 1, .orientation = NONE},
        [6] = {.type = IMPROPER_ROTATION, .order = 4, .power = 1, .orientation = HORIZONTAL},
        [7] = {.type = IMPROPER_ROTATION, .order = 6, .power = 1, .orientation = NONE},
        [8] = {.type = REFLECTION, .order = 1, .power = 1, .orientation = HORIZONTAL},
        [9] = {.type = REFLECTION, .order = 1, .power = 1, .orientation = DIHEDRAL}
    };
    
    const double table[][10] = {
        [0] = {1,  1,  1,  1,  1,  1,  1,  1,  1,  1},
        [1] = {1,  1,  1,  1,  1, -1, -1, -1, -1, -1},
        [2] = {1, -1,  1,  1, -1,  1, -1,  1,  1, -1},
        [3] = {1, -1,  1,  1, -1, -1,  1, -1, -1,  1},
        [4] = {2,  0,  2, -1,  0,  2,  0, -1,  2,  0},
        [5] = {2,  0,  2, -1,  0, -2,  0,  1, -2,  0},
        [6] = {3,  1, -1,  0, -1,  3,  1,  0, -1, -1},
        [7] = {3,  1, -1,  0, -1, -3, -1,  0,  1,  1},
        [8] = {3, -1, -1,  0,  1,  3, -1,  0, -1,  1},
        [9] = {3, -1, -1,  0,  1, -3,  1,  0,  1, -1}
    };
    
    
    
    const char *tname[10] = {"A1g", "A1u", "A2g", "A2u", "Eg", "Eu", "T1g", "T1u", "T2g", "T2u"};
    const int tdim[10] = {1,1,1,1,2,2,3,3,3,3};
    const int tred[10] = {1,1,1,1,1,1,1,1,1,1};

    
    if(MSYM_SUCCESS != (ret = getPredefinedCharacterTable(sopsl, sops, 10, tsops, tname, tdim, tred, table, ct))) goto err;
    
err:
    return ret;
}

#define C4PI (-1.61803398874989484820458683436563811772030917980576286213544) //(2*cos(4*M_PI/5))
#define C2PI (0.618033988749894848204586834365638117720309179805762862135448) //(2*cos(2*M_PI/5))

msym_error_t getCharacterTableI(int sopsl, msym_symmetry_operation_t sops[sopsl], msym_character_table_t *ct){
    msym_error_t ret = MSYM_SUCCESS;
    
    const msym_symmetry_operation_t tsops[5] = {
        [0] = {.type = IDENTITY, .order = 1, .power = 1, .orientation = NONE},
        [1] = {.type = PROPER_ROTATION, .order = 2, .power = 1, .orientation = NONE},
        [2] = {.type = PROPER_ROTATION, .order = 3, .power = 1, .orientation = NONE},
        [3] = {.type = PROPER_ROTATION, .order = 5, .power = 1, .orientation = NONE},
        [4] = {.type = PROPER_ROTATION, .order = 5, .power = 2, .orientation = NONE},
    };
    
    const char *tname[5] = {"A","T1","T2","G","H"};
    const int tdim[5] = {1,3,3,4,5};
    const int tred[5] = {1,1,1,1,1};
    
    
    //         E     C2     C3    C5     C52
    const double table[][5] = {
        [0] = {1,     1,     1,    1,     1},
        [1] = {3,    -1,     0,   -C4PI, -C2PI},
        [2] = {3,    -1,     0,   -C2PI, -C4PI},
        [3] = {4,     0,     1,   -1,    -1},
        [4] = {5,     1,    -1,    0,     0}
    };

    if(MSYM_SUCCESS != (ret = getPredefinedCharacterTable(sopsl, sops, 5, tsops, tname, tdim, tred, table, ct))) goto err;
    
err:
    return ret;
}




msym_error_t getCharacterTableIh(int sopsl, msym_symmetry_operation_t sops[sopsl], msym_character_table_t *ct){
    msym_error_t ret = MSYM_SUCCESS;
    
    const msym_symmetry_operation_t tsops[10] = {
        [0] = {.type = IDENTITY, .order = 1, .power = 1, .orientation = NONE},
        [1] = {.type = PROPER_ROTATION, .order = 2, .power = 1, .orientation = NONE},
        [2] = {.type = REFLECTION, .order = 1, .power = 1, .orientation = NONE},
        [3] = {.type = IMPROPER_ROTATION, .order = 6, .power = 1, .orientation = NONE},
        [4] = {.type = PROPER_ROTATION, .order = 5, .power = 1, .orientation = NONE},
        [5] = {.type = IMPROPER_ROTATION, .order = 10, .power = 1, .orientation = NONE},
        [6] = {.type = PROPER_ROTATION, .order = 5, .power = 2, .orientation = NONE},
        [7] = {.type = INVERSION, .order = 1, .power = 1, .orientation = NONE},
        [8] = {.type = PROPER_ROTATION, .order = 3, .power = 1, .orientation = NONE},
        [9] = {.type = IMPROPER_ROTATION, .order = 10, .power = 3, .orientation = NONE},
    };
    
    const char *tname[10] = {"Ag","Au","T1g","T1u","T2g","T2u","Gg","Gu","Hg","Hu"};
    const int tdim[10] = {1,1,3,3,3,3,4,4,5,5};
    const int tred[10] = {1,1,1,1,1,1,1,1,1,1};
    
    //          E    C2   R    S6   C5   S10  C52  i    C3   S103
    static const double table[][10] =
    {
        [0] = {1,     1,      1,      1,      1,      1,      1,      1,      1,      1     },
        [1] = {1,     1,     -1,     -1,      1,     -1,      1,     -1,      1,     -1     },
        [2] = {3,    -1,     -1,      0,     -C4PI,  -C2PI,  -C2PI,   3,      0,     -C4PI  },
        [3] = {3,    -1,      1,      0,     -C4PI,   C2PI,  -C2PI,  -3,      0,      C4PI  },
        [4] = {3,    -1,     -1,      0,     -C2PI,  -C4PI,  -C4PI,   3,      0,     -C2PI  },
        [5] = {3,    -1,      1,      0,     -C2PI,   C4PI,  -C4PI,  -3,      0,      C2PI  },
        [6] = {4,     0,      0,      1,     -1,     -1,     -1,      4,      1,     -1     },
        [7] = {4,     0,      0,     -1,     -1,      1,     -1,     -4,      1,      1     },
        [8] = {5,     1,      1,     -1,      0,      0,      0,      5,     -1,      0     },
        [9] = {5,     1,     -1,      1,      0,      0,      0,     -5,     -1,      0     }
    };
    
    if(MSYM_SUCCESS != (ret = getPredefinedCharacterTable(sopsl, sops, 10, tsops, tname, tdim, tred, table, ct))) goto err;
    
err:
    return ret;
}

msym_error_t getPredefinedCharacterTable(int sopsl, msym_symmetry_operation_t sops[sopsl], int l, const msym_symmetry_operation_t tsops[l], const char *tname[l], const int tdim[l], const int tred[l], const double (*table)[l], msym_character_table_t *ct){
    msym_error_t ret = MSYM_SUCCESS;
    
    double (*ctable)[l] = ct->table;
    
    if(ct->d != l){
        msymSetErrorDetails("Unexpected size of character table %d != %d",l,ct->d);
        ret = MSYM_INVALID_CHARACTER_TABLE;
        goto err;
    }
    
    for(int i = 0; i < l;i++){
        ct->s[i].d = tdim[i];
        ct->s[i].r = tred[i];
        snprintf(ct->s[i].name, sizeof(ct->s[i].name), "%s",tname[i]);
        msym_symmetry_operation_t *sop = NULL;
        for(sop = sops;sop < (sops + sopsl);sop++){
            if(tsops[i].type == sop->type && tsops[i].order == sop->order && tsops[i].power == sop->power && tsops[i].orientation == sop->orientation){
                int cla = sop->cla;
                if(cla >= l){
                    msymSetErrorDetails("Conjugacy class exceeds character table size %d >= %d",sop->cla,l);
                    ret = MSYM_INVALID_CHARACTER_TABLE;
                    goto err;
                }
                
                for(int j = 0;j < l;j++){
                    ctable[j][cla] = table[j][i];
                }
                break;
            }
        }
        if(sop >= (sops + sopsl)){
            msymSetErrorDetails("Could not find representative symmetry operation when generating character table");
            ret = MSYM_INVALID_CHARACTER_TABLE;
            goto err;
        }
    }
    
err:
    return ret;
}


msym_error_t getRepresentationName(msym_point_group_type_t type, int n, msym_representation_t *rep, int l, char name[l]){
    msym_error_t ret = MSYM_SUCCESS;
    if(rep->d < 1 || rep->d > 5 || abs(rep->eig.p) > 1 || abs(rep->eig.v) > 1 || abs(rep->eig.h) > 1 || abs(rep->eig.i) > 1) {
        ret = MSYM_INVALID_CHARACTER_TABLE;
        msymSetErrorDetails("Invalid character table represenation");
        goto err;
    }
    
    int eindex[4] = {rep->eig.p, rep->eig.h,rep->eig.v,rep->eig.i};
    switch (type) {
        case MSYM_POINT_GROUP_TYPE_Ci :
            eindex[1] = eindex[2] = 0;
            break;
        case MSYM_POINT_GROUP_TYPE_Cs :
            eindex[2] = eindex[3] = 0;
            break;
        case MSYM_POINT_GROUP_TYPE_Cn :
            eindex[1] = eindex[2] = eindex[3] = 0;
            break;
        case MSYM_POINT_GROUP_TYPE_Cnv :
            eindex[1] = eindex[3] = 0;
            break;
        case MSYM_POINT_GROUP_TYPE_Cnh :
            if(n & 1){eindex[3] = 0;}
            else {eindex[1] = 0;}
            eindex[2] = 0;
            break;
        case MSYM_POINT_GROUP_TYPE_Sn :
            if(!((n >> 1) & 1)){eindex[3] = 0; eindex[0] = rep->eig.h;}
            eindex[1] = eindex[2] = 0;
            break;
            break;
        case MSYM_POINT_GROUP_TYPE_Dn  :
            if(n == 2 && eindex[0]+eindex[2] == 0){
                eindex[0] = eindex[0]*eindex[2];
                eindex[2] = ((eindex[2] + 1) >> 1) + 1;
            }
            eindex[1] = 0;
            eindex[3] = 0;
            break;
        case MSYM_POINT_GROUP_TYPE_Dnd :
            if(!(n & 1)){eindex[3] = 0; eindex[0] = rep->eig.h;}
            eindex[1] = 0;
            break;
        case MSYM_POINT_GROUP_TYPE_Dnh :
            if(n == 2 && eindex[0]+eindex[2] == 0){
                eindex[0] = eindex[0]*eindex[2];
                eindex[2] = ((eindex[2] + 1) >> 1) + 1;
            }
            if(n & 1){eindex[3] = 0;}
            else {eindex[1] = 0;}
            break;
        default:
            break;
    }
    
    char types[] = {'A','B','E','T','G','H'}, *si[] = {"u","","g"}, *sv[] = {"2", "", "1", "3"}, *sh[] = {"''", "", "'"};
    char rtype = rep->d == 1 ? types[(1 - eindex[0]) >> 1] : types[rep->d];
    if(rep->d == 1){
        snprintf(name,sizeof(char[l]),"%c%s%s%s",rtype,sv[eindex[2]+1],si[eindex[3]+1],sh[eindex[1]+1]);
    }
    else if (rep->eig.l >  0){
        snprintf(name,sizeof(char[l]),"%c%d%s%s",rtype,rep->eig.l,si[eindex[3]+1],sh[eindex[1]+1]);
    } else {
        snprintf(name,sizeof(char[l]),"%c%s%s",rtype,si[eindex[3]+1],sh[eindex[1]+1]);
    }
err:
    return ret;
}


