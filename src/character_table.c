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
#include "linalg.h"

typedef struct _msym_representation {
    enum {IRREDUCIBLE, REDUCIBLE} type;
    int d;
    struct {
        int p, v, h, i, l;
    } eig;
    char name[8];
} msym_representation_t;

msym_error_t verifyCharacterTable(msym_character_table_t *ct);

msym_error_t setRepresentationName(msym_representation_t *rep);
msym_error_t representationCharacter(int n, msym_symmetry_operation_t *sop, msym_representation_t *rep, double *c);

msym_error_t getRepresentationsCn(int n, int rl, msym_representation_t rep[rl]);
msym_error_t getRepresentationsCnh(int n, int rl, msym_representation_t rep[rl]);
msym_error_t getRepresentationsCnv(int n, int rl, msym_representation_t rep[rl]);
msym_error_t getRepresentationsDn(int n, int rl, msym_representation_t rep[rl]);
msym_error_t getRepresentationsDnh(int n, int rl, msym_representation_t rep[rl]);
msym_error_t getRepresentationsDnd(int n, int rl, msym_representation_t rep[rl]);
msym_error_t getRepresentationsUnknown(int n, int rl, msym_representation_t rep[rl]);


void decomposeRepresentation(CharacterTable *ct, double rspan[ct->l], double dspan[ct->l]){
    int order = 0;
    memset(dspan,0, sizeof(double[ct->l]));
    for(int k = 0;k < ct->l;k++){
        order += ct->classc[k];
        for(int j = 0; j < ct->l;j++) dspan[k] += ct->classc[j]*rspan[j]*ct->irrep[k].v[j];
    }
    for(int k = 0;k < ct->l;k++) dspan[k] /= order;
}

void decomposeRepresentation2(msym_character_table_t *ct, double rspan[ct->d], double dspan[ct->d]){
    int order = 0;
    double (*ctable)[ct->d] = ct->table;
    memset(dspan,0, sizeof(double[ct->d]));
    
    for(int k = 0;k < ct->d;k++){
        order += ct->classc[k];
        for(int j = 0; j < ct->d;j++) dspan[k] += ct->classc[j]*rspan[j]*ctable[k][j];
    }
    for(int k = 0;k < ct->d;k++) dspan[k] /= order;
}

void directProduct(int l, IrreducibleRepresentation *irrep1, IrreducibleRepresentation *irrep2, double pspan[l]){
    for(int i = 0;i < l;i++) pspan[i] = irrep1->v[i]*irrep2->v[i];
}

void directProduct2(int l, double irrep1[l], double irrep2[l], double pspan[l]){
    for(int i = 0;i < l;i++) pspan[i] = irrep1[i]*irrep2[i];
}

msym_error_t characterTableUnknown(int n, CharacterTable *ct){
    printf("WARNING UNKOWN\n");
    msymSetErrorDetails("Character table unknown");
    //return MSYM_INVALID_CHARACTER_TABLE;
    return MSYM_SUCCESS;
}

msym_error_t characterTableTd(int n, CharacterTable *ct){
    ct->l = sizeof(TdIrrep)/sizeof(TdIrrep[0]);
    ct->irrep = malloc(ct->l*sizeof(IrreducibleRepresentation));
    for(int i = 0; i < ct->l;i++){
        enum IrreducibleRepresentationEnum irrep = TdIrrep[i];
        ct->irrep[i].name = IrreducibleRepresentationName[irrep];
        ct->irrep[i].v = TdTable[irrep];
        ct->irrep[i].d = Degeneracy[irrep];
        ct->irrep[i].l = sizeof(TdTable[irrep])/sizeof(TdTable[irrep][0]);
    }
    return MSYM_SUCCESS;
}

msym_error_t characterTableIh(int n, CharacterTable *ct){
    ct->l = sizeof(IhIrrep)/sizeof(IhIrrep[0]);
    ct->irrep = malloc(ct->l*sizeof(IrreducibleRepresentation));
    for(int i = 0; i < ct->l;i++){
        enum IrreducibleRepresentationEnum irrep = IhIrrep[i];
        ct->irrep[i].name = IrreducibleRepresentationName[irrep];
        ct->irrep[i].v = IhTable[irrep];
        ct->irrep[i].d = Degeneracy[irrep];
        ct->irrep[i].l = sizeof(IhTable[irrep])/sizeof(IhTable[irrep][0]);
    }
    return MSYM_SUCCESS;
}

msym_error_t characterTableCnv(int n, CharacterTable *ct){
    msym_error_t ret = MSYM_SUCCESS;
    switch(n) {
        case 3 : {
            ct->l = sizeof(C3vIrrep)/sizeof(C3vIrrep[0]);
            ct->irrep = malloc(ct->l*sizeof(IrreducibleRepresentation));
            for(int i = 0; i < ct->l;i++){
                enum IrreducibleRepresentationEnum irrep = C3vIrrep[i];
                ct->irrep[i].name = IrreducibleRepresentationName[irrep];
                ct->irrep[i].v = C3vTable[irrep];
                ct->irrep[i].d = Degeneracy[irrep];
                ct->irrep[i].l = sizeof(C3vTable[irrep])/sizeof(C3vTable[irrep][0]);
            }
            break;
        }
        case 4 : {
            ct->l = sizeof(C4vIrrep)/sizeof(C4vIrrep[0]);
            ct->irrep = malloc(ct->l*sizeof(IrreducibleRepresentation));
            for(int i = 0; i < ct->l;i++){
                enum IrreducibleRepresentationEnum irrep = C4vIrrep[i];
                ct->irrep[i].name = IrreducibleRepresentationName[irrep];
                ct->irrep[i].v = C4vTable[irrep];
                ct->irrep[i].d = Degeneracy[irrep];
                ct->irrep[i].l = sizeof(C4vTable[irrep])/sizeof(C4vTable[irrep][0]);
            }
            break;
        }

        default:
            msymSetErrorDetails("Cannot find C%dv character table",n);
            ret = MSYM_INVALID_CHARACTER_TABLE;
    }
    
    return ret;
}

msym_error_t characterTableCnh(int n, CharacterTable *ct){
    msym_error_t ret = MSYM_SUCCESS;
    switch(n) {
        case 2 : {
            ct->l = sizeof(C2hIrrep)/sizeof(C2hIrrep[0]);
            ct->irrep = malloc(ct->l*sizeof(IrreducibleRepresentation));
            for(int i = 0; i < ct->l;i++){
                enum IrreducibleRepresentationEnum irrep = C2hIrrep[i];
                ct->irrep[i].name = IrreducibleRepresentationName[irrep];
                ct->irrep[i].v = C2hTable[irrep];
                ct->irrep[i].d = Degeneracy[irrep];
                ct->irrep[i].l = sizeof(C2hTable[irrep])/sizeof(C2hTable[irrep][0]);
            }
            break;
        }
            
        default:
            msymSetErrorDetails("Cannot find C%dh character table",n);
            ret = MSYM_INVALID_CHARACTER_TABLE;
    }
    
    return ret;
}

msym_error_t characterTableDnh(int n, CharacterTable *ct){
    msym_error_t ret = MSYM_SUCCESS;
    switch(n) {
        case 2 : {
            ct->l = sizeof(D2hIrrep)/sizeof(D2hIrrep[0]);
            ct->irrep = malloc(ct->l*sizeof(IrreducibleRepresentation));
            for(int i = 0; i < ct->l;i++){
                enum IrreducibleRepresentationEnum irrep = D2hIrrep[i];
                ct->irrep[i].name = IrreducibleRepresentationName[irrep];
                ct->irrep[i].v = D2hTable[irrep];
                ct->irrep[i].d = Degeneracy[irrep];
                ct->irrep[i].l = sizeof(D2hTable[irrep])/sizeof(D2hTable[irrep][0]);
            }
            break;
        }
        case 4 : {
            ct->l = sizeof(D4hIrrep)/sizeof(D4hIrrep[0]);
            ct->irrep = malloc(ct->l*sizeof(IrreducibleRepresentation));
            for(int i = 0; i < ct->l;i++){
                enum IrreducibleRepresentationEnum irrep = D4hIrrep[i];
                ct->irrep[i].name = IrreducibleRepresentationName[irrep];
                ct->irrep[i].v = D4hTable[irrep];
                ct->irrep[i].d = Degeneracy[irrep];
                ct->irrep[i].l = sizeof(D4hTable[irrep])/sizeof(D4hTable[irrep][0]);
            }
            break;
        }
        case 6 : {
            ct->l = sizeof(D6hIrrep)/sizeof(D6hIrrep[0]);
            ct->irrep = malloc(ct->l*sizeof(IrreducibleRepresentation));
            for(int i = 0; i < ct->l;i++){
                enum IrreducibleRepresentationEnum irrep = D6hIrrep[i];
                ct->irrep[i].name = IrreducibleRepresentationName[irrep];
                ct->irrep[i].v = D6hTable[irrep];
                ct->irrep[i].d = Degeneracy[irrep];
                ct->irrep[i].l = sizeof(D6hTable[irrep])/sizeof(D6hTable[irrep][0]);
            }
            break;
        }
        default:
            msymSetErrorDetails("Cannot find D%dh character table",n);
            ret = MSYM_INVALID_CHARACTER_TABLE;
    }
    return ret;
}

double getCharacterCnv(int n, int k){
    return 0;
}

msym_error_t generateCharacterTable(msym_point_group_type_t type, int n, int sopsl, msym_symmetry_operation_t sops[sopsl], msym_character_table_t **oct){
    msym_error_t ret = MSYM_SUCCESS;
    msym_character_table_t *ct = calloc(1, sizeof(msym_character_table_t));
    ct->d = sops[sopsl-1].cla + 1;
    ct->table = calloc(ct->d, sizeof(double[ct->d]));
    ct->classc = calloc(ct->d, sizeof(int));
    ct->s = calloc(ct->d, sizeof(msym_symmetry_species_t));
    msym_representation_t *rep = calloc(ct->d, sizeof(msym_representation_t));
    double (*table)[ct->d] = (double (*)[ct->d]) ct->table;
    
    const struct _fmap {
        msym_point_group_type_t type;
        msym_error_t (*f)(int, int, msym_representation_t *);
    } fmap[18] = {
        
        [ 0] = {POINT_GROUP_Ci,  getRepresentationsUnknown},
        [ 1] = {POINT_GROUP_Cs,  getRepresentationsUnknown},
        [ 2] = {POINT_GROUP_Cn,  getRepresentationsCn},
        [ 3] = {POINT_GROUP_Cnh, getRepresentationsCnh},
        [ 4] = {POINT_GROUP_Cnv, getRepresentationsCnv},
        [ 5] = {POINT_GROUP_Dn,  getRepresentationsDn},
        [ 6] = {POINT_GROUP_Dnh, getRepresentationsDnh},
        [ 7] = {POINT_GROUP_Dnd, getRepresentationsDnd},
        [ 8] = {POINT_GROUP_Sn,  getRepresentationsUnknown},
        [ 9] = {POINT_GROUP_T,   getRepresentationsUnknown},
        [10] = {POINT_GROUP_Td,  getRepresentationsUnknown},
        [11] = {POINT_GROUP_Th,  getRepresentationsUnknown},
        [12] = {POINT_GROUP_O,   getRepresentationsUnknown},
        [13] = {POINT_GROUP_Oh,  getRepresentationsUnknown},
        [14] = {POINT_GROUP_I,   getRepresentationsUnknown},
        [15] = {POINT_GROUP_Ih,  getRepresentationsUnknown},
        [16] = {POINT_GROUP_K,   getRepresentationsUnknown},
        [17] = {POINT_GROUP_Kh,  getRepresentationsUnknown}
    };
    
    int fi, fil = sizeof(fmap)/sizeof(fmap[0]);
    for(fi = 0; fi < fil;fi++){
        if(fmap[fi].type == type) {
            if(MSYM_SUCCESS != (ret = fmap[fi].f(n,ct->d,rep))) goto err;
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
    
    for(int i = 0;i < ct->d;i++){
        snprintf(ct->s[i].name, sizeof(ct->s[i].name), "%s",rep[i].name);
        ct->s[i].d = rep[i].d;
        int nc = -1;
        for(int j = 0;j < sopsl;j++){
            if(nc < sops[j].cla){
                nc = sops[j].cla;
                if(MSYM_SUCCESS != (ret = representationCharacter(n,&sops[j],&rep[i],&table[i][nc]))) goto err;
            }
        }
    }
    
    
    
    int nc = -1;
    for(int j = 0;j < sopsl;j++){
        //if(j == 0) printf("\t");
        if(nc < sops[j].cla){
            nc = sops[j].cla;
            char buf[12];
            symmetryOperationName(&sops[j], 12, buf);
            printf("%d%s",ct->classc[sops[j].cla],buf);
            /*if(sops[j].order == 2 && sops[j].type == PROPER_ROTATION && sops[j].v[2] != 1.0){
                if(sops[j].power == 1) printf("'");
                else printf("''");
            }
            if(sops[j].type == REFLECTION){
                if(sops[j].v[2] == 1.0) printf("h");
                else if(sops[j].power == 1) printf("v");
                else printf("d");
            }*/
            printf("\t\t");
        }
    }
    printf("\n");
    for(int i = 0;i < ct->d;i++){
        
        printf("%s\t",rep[i].name);
        for(int j = 0;j < ct->d;j++){
            printf("% .3lf\t\t",table[i][j]);
        }
        printf("\n");
    }
    
    if(MSYM_SUCCESS != (ret = verifyCharacterTable(ct))) goto err;
    
    *oct = ct;
    
    free(rep);
    return ret;
err:
    free(ct->table);
    free(ct->s);
    free(ct->classc);
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

msym_error_t getRepresentationsUnknown(int n, int rl, msym_representation_t rep[rl]){
    msymSetErrorDetails("Character table representation NYI");
    return MSYM_INVALID_CHARACTER_TABLE;
}

msym_error_t getRepresentationsCn(int n, int rl, msym_representation_t rep[rl]){
    msym_error_t ret = MSYM_SUCCESS;
    int r = 0;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
    if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
    r++;
    if(~n & 1){
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
        rep[r].eig.p = -1;
        if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;;
        r++;
    }
    for(int i = 1;r < rl;i++, r++){
        rep[r].type = REDUCIBLE;
        rep[r].d = 2;
        rep[r].eig.l = i;
        rep[r].eig.p = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
        if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;;
    }
    
    return ret;
err:
    return ret;
}

msym_error_t getRepresentationsCnh(int n, int rl, msym_representation_t rep[rl]){
    msym_error_t ret = MSYM_SUCCESS;
    int r = 0;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
    if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
    r++;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.v = 1;
    rep[r].eig.h = rep[r].eig.i = -1;
    if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
    r++;
    if(~n & 1){
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.v = rep[r].eig.i = rep[r].eig.h = 1;
        rep[r].eig.p = -1;
        if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
        r++;
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.v = 1;
        rep[r].eig.p = rep[r].eig.i = rep[r].eig.h = -1;
        if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
        r++;
    }
    for(int i = 1;r < rl;i++, r++){
        rep[r].type = REDUCIBLE;
        rep[r].d = 2;
        rep[r].eig.l = i;
        rep[r].eig.p = rep[r].eig.v = rep[r].eig.h = 1;
        rep[r].eig.i = -1;
        if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;;
        r++;
        rep[r].type = REDUCIBLE;
        rep[r].d = 2;
        rep[r].eig.l = i;
        rep[r].eig.p = rep[r].eig.v = rep[r].eig.i = 1;
        rep[r].eig.h = -1;
        if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;;
    }
    
    return ret;
err:
    return ret;
}

msym_error_t getRepresentationsCnv(int n, int rl, msym_representation_t rep[rl]){
    msym_error_t ret = MSYM_SUCCESS;
    int r = 0;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
    if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
    r++;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.h = rep[r].eig.i = 1;
    rep[r].eig.v = -1;
    if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
    r++;
    if(~n & 1){
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.v = rep[r].eig.i = rep[r].eig.h = 1;
        rep[r].eig.p = -1;
        if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
        r++;
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.h = rep[r].eig.i = 1 ;
        rep[r].eig.p = rep[r].eig.v = -1;
        if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
        r++;
    }
    for(int i = 1;r < rl;i++, r++){
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 2;
        rep[r].eig.l = i;
        rep[r].eig.p = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
        if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;;
    }
    
    return ret;
err:
    return ret;
}

msym_error_t getRepresentationsDn(int n, int rl, msym_representation_t rep[rl]){
    msym_error_t ret = MSYM_SUCCESS;
    int r = 0;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
    if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
    r++;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.h = rep[r].eig.i = 1;
    rep[r].eig.v = -1;
    if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
    r++;
    if(~n & 1){
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.v = rep[r].eig.i = rep[r].eig.h = 1;
        rep[r].eig.p = -1;
        if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
        r++;
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.h = rep[r].eig.i = 1 ;
        rep[r].eig.p = rep[r].eig.v = -1;
        if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
        r++;
    }
    for(int i = 1;r < rl;i++, r++){
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 2;
        rep[r].eig.l = i;
        rep[r].eig.p = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
        if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;;
    }
    
    return ret;
err:
    return ret;
}

msym_error_t getRepresentationsDnh(int n, int rl, msym_representation_t rep[rl]){
    msym_error_t ret = MSYM_SUCCESS;
    int r = 0;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
    if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
    r++;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.h = rep[r].eig.i = 1;
    rep[r].eig.v = -1;
    if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
    r++;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.v = 1;
    rep[r].eig.h = rep[r].eig.i = -1;
    if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
    r++;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = 1;
    rep[r].eig.h = rep[r].eig.i = rep[r].eig.v = -1;
    if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
    r++;
    if(~n & 1){
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.v = rep[r].eig.i = 1;
        rep[r].eig.p = rep[r].eig.h = -1;
        if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
        r++;
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.i = 1;
        rep[r].eig.p = rep[r].eig.h = rep[r].eig.v = -1;
        if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
        r++;
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.h = rep[r].eig.v = 1;
        rep[r].eig.p = rep[r].eig.i = -1;
        if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
        r++;
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.h = 1;
        rep[r].eig.p = rep[r].eig.i = rep[r].eig.v = -1;
        if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
        r++;
    }
    for(int i = 1;r < rl;i++, r++){
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 2;
        rep[r].eig.l = i;
        rep[r].eig.p = rep[r].eig.v = rep[r].eig.h = 1;
        rep[r].eig.i = 1 - ((i & 1) << 1);
        if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;;
        r++;
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 2;
        rep[r].eig.l = i;
        rep[r].eig.p = rep[r].eig.v = 1;
        rep[r].eig.h = -1;
        rep[r].eig.i = -1 + ((i & 1) << 1);
        if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;;
    }
    
    return ret;
err:
    return ret;
}

msym_error_t getRepresentationsDnd(int n, int rl, msym_representation_t rep[rl]){
    msym_error_t ret = MSYM_SUCCESS;
    int r = 0;
    printf("\n\n\nprimary is a sn axis here e messed up!\n\n\n\n");
    exit(1);
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
    if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
    r++;
    rep[r].type = IRREDUCIBLE;
    rep[r].d = 1;
    rep[r].eig.p = rep[r].eig.l = rep[r].eig.h = rep[r].eig.i = 1;
    rep[r].eig.v = -1;
    if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
    r++;
    if(~n & 1){
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.v = rep[r].eig.i = 1;
        rep[r].eig.h = rep[r].eig.p = -1;
        if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
        r++;
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.i = 1 ;
        rep[r].eig.h = rep[r].eig.p = rep[r].eig.v = -1;
        if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
        r++;
        for(int i = 1;r < rl;i++, r++){
            rep[r].type = IRREDUCIBLE;
            rep[r].d = 2;
            rep[r].eig.l = i;
            rep[r].eig.p = rep[r].eig.v = rep[r].eig.h = rep[r].eig.i = 1;
            if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;;
        }
    } else {
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.v = rep[r].eig.p = 1;
        rep[r].eig.h = rep[r].eig.i = -1;
        if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
        r++;
        rep[r].type = IRREDUCIBLE;
        rep[r].d = 1;
        rep[r].eig.l = rep[r].eig.p = 1 ;
        rep[r].eig.h = rep[r].eig.i = rep[r].eig.v = -1;
        if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
        r++;
        for(int i = 1;r < rl;i++, r++){
            rep[r].type = IRREDUCIBLE;
            rep[r].d = 2;
            rep[r].eig.l = i;
            rep[r].eig.p = rep[r].eig.v = rep[r].eig.i = 1;
            rep[r].eig.h = -1;
            if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;
            r++;
            rep[r].type = IRREDUCIBLE;
            rep[r].d = 2;
            rep[r].eig.l = i;
            rep[r].eig.p = rep[r].eig.v = rep[r].eig.h = 1;
            rep[r].eig.i = -1;
            if(MSYM_SUCCESS != (ret = setRepresentationName(&rep[r]))) goto err;;
        }
        
    }
    
    
    return ret;
err:
    return ret;
}


msym_error_t representationCharacter(int n, msym_symmetry_operation_t *sop, msym_representation_t *rep, double *c){
    msym_error_t ret = MSYM_SUCCESS;
    double x = 0;
    
    //printSymmetryOperation(sop);
    if(sop->orientation == HORIZONTAL){
    //if(sop->v[2] == 1.0){
        switch(rep->d)
        {
            case 1: {
                switch (sop->type) {
                    case IDENTITY           : x = 1; break;
                    case REFLECTION         : x = rep->eig.h; break;
                    case INVERSION          : x = rep->eig.i; break;
                    case PROPER_ROTATION    : x = ((n/sop->order) & 1) ? rep->eig.p : 1 ; break;
                    case IMPROPER_ROTATION  : x = rep->eig.h*(((n/sop->order) & 1) ? rep->eig.p : 1); break;
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
                    case REFLECTION         : x = rep->eig.p == 1 || sop->orientation == VERTICAL ? rep->eig.h*rep->eig.v : -rep->eig.h*rep->eig.v; break;
                    case PROPER_ROTATION    : x = rep->eig.p == 1 || sop->orientation == VERTICAL ? rep->eig.v : -rep->eig.v; break;
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

msym_error_t symmetryOperationCharacterDnh(msym_symmetry_operation_t *sop, msym_representation_t *irrep, int n, int *character){
    msym_error_t ret = MSYM_SUCCESS;
    int x = 0;
    switch (sop->type) {
        case IDENTITY: x = irrep->d; break;
        case REFLECTION : {
            if(sop->v[2] == 1.0){
                x = irrep->eig.h;
            } else if(sop->power == 1){
                x = irrep->eig.v*irrep->eig.h;
            } else if(sop->power == -1){
                x = irrep->eig.v*irrep->eig.h*irrep->eig.p;
            }
            break;
        }
        case PROPER_ROTATION : {
            if(sop->v[2] == 1.0){
                x = 0;
                printf("asddiiv\n");
                exit(1);
            }
            if(irrep->d == 1 && sop->power == 1){
                x = (irrep->eig.p + ((~sop->order & 1) << 1));
            } else if(irrep->d == 1 && sop->power == -1){
                x = (irrep->eig.p + ((~sop->order & 1) << 1));
            } else {
                x = 2*cos(irrep->eig.l*2*M_PI/sop->order);
            }
            
        }
        case IMPROPER_ROTATION : {
            if(irrep->d == 1){
                x = (irrep->eig.p + ((~sop->order & 1) << 1))*irrep->eig.h;
            } else {
                x = 2*cos(irrep->eig.l*2*M_PI/sop->order)*irrep->eig.h;
            }
        }
            
        default:
            break;
    }
    
    *character = x;
    
    return ret;
}


msym_error_t setRepresentationName(msym_representation_t *rep){
    msym_error_t ret = MSYM_SUCCESS;
    if(rep->d < 1 || rep->d > 5 || abs(rep->eig.p) > 1 || abs(rep->eig.v) > 1 || abs(rep->eig.h) > 1 || abs(rep->eig.i) > 1) {
        ret = MSYM_INVALID_CHARACTER_TABLE;
        msymSetErrorDetails("Invalid irreducible represenation");
        goto err;
    }
    
    char types[] = {'A','B','E','T','G','H'}, *si[] = {"u","","g"}, *sv[] = {"2", "", "1"}, *sh[] = {"''", "", "'"};
    char type = rep->d == 1 ? types[(1 - rep->eig.p) >> 1] : types[rep->d];

    if(rep->d == 1){
        snprintf(rep->name,sizeof(rep->name),"%c%s%s%s",type,sv[rep->eig.v+1],si[rep->eig.i+1],sh[rep->eig.h+1]);
    }
    else if (rep->eig.l >  0){
        snprintf(rep->name,sizeof(rep->name),"%s%c%d%s%s",rep->type == IRREDUCIBLE ? "" : "*",type,rep->eig.l,si[rep->eig.i+1],sh[rep->eig.h+1]);
    } else {
        snprintf(rep->name,sizeof(rep->name),"%s%c%s%s",rep->type == IRREDUCIBLE ? "" : "*",type,si[rep->eig.i+1],sh[rep->eig.h+1]);
    }
err:
    return ret;
}

// \u03C0 = pi
// \u03C3 = sigma
void printCharacterTable(CharacterTable *ct){
    printf("Character Table:\n");
    for(int i = 0; i < ct->l;i++){
        printf("\t %d%s",ct->classc[i],ct->name[i]);
    }
    printf("\n");
    for(int i = 0; i < ct->l;i++){
        printf("%s:\t",ct->irrep[i].name);
        for(int j = 0; j < ct->irrep[i].l; j++){
            char *pre = signbit(ct->irrep[i].v[j]) == 1 ? "" : " ";
            printf("%s%.3lf\t",pre,ct->irrep[i].v[j]);
        }
        printf("\n");
    
    }
}

void convertNewCharacterTable(msym_point_group_t *pg){
    msym_character_table_t *ct = malloc(sizeof(msym_character_table_t));
    ct->d = pg->ct->l;
    double (*t)[ct->d] = malloc(sizeof(double[ct->d][ct->d]));
    msym_symmetry_species_t *ssp = malloc(sizeof(msym_symmetry_species_t[ct->d]));
    ct->s = ssp;
    ct->table = t;
    ct->classc = pg->ct->classc;
    for(int i = 0;i < pg->ct->l;i++){
        memcpy(t[i], pg->ct->irrep[i].v, sizeof(double[ct->d]));
        ssp[i].d = pg->ct->irrep[i].d;
        sprintf(ssp[i].name, "%s",pg->ct->irrep[i].name);
    }
    pg->ct2 = ct;
}
