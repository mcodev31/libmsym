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

#include "character_table.h"
#include "point_group.h"

msym_error_t characterTableUnknown(int n, CharacterTable *ct){
    msymSetErrorDetails("Character table unknown");
    return MSYM_INVALID_CHARACTER_TABLE;
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


void testirrepname(){
    char buf[16];
    for(int d = 1;d < 6;d++)
        for(int p = 0;p < 3;p++)
            for(int v = 0;v < 3;v++)
                for(int h = 0;h < 3;h++)
                    for(int l = 0;l < 3+d;l++)
                        for(int i = 0;i < 3;i++){
                            msym_irreducible_representation_t irrep;
                            irrep.d = d;
                            irrep.eig.p = p-1;
                            irrep.eig.v = v-1;
                            irrep.eig.i = i-1;
                            irrep.eig.h = h-1;
                            irrep.eig.l = l;
                            irreducibleRepresenationName(&irrep,16,buf);
                            printf("d = %d, p = %d, i = %d, h = %d v = %d -> ",
                                   irrep.d,irrep.eig.p,irrep.eig.i,irrep.eig.h,irrep.eig.v);
                            printf("%s\n",buf);
                        }
    
}

#define CHARACTER_THRESHOLD 0.001
msym_error_t symmetryOperationCharacterDnh(msym_symmetry_operation_t *sop, msym_irreducible_representation_t *irrep, int n, int *character){
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


msym_error_t irreducibleRepresenationName(msym_irreducible_representation_t *irrep, int l, char name[l]){
    msym_error_t ret = MSYM_SUCCESS;
    if(irrep->d < 1 || irrep->d > 5 || abs(irrep->eig.p) > 1 || abs(irrep->eig.v) > 1 || abs(irrep->eig.h) > 1 || abs(irrep->eig.i) > 1) {
        ret = MSYM_INVALID_CHARACTER_TABLE;
        msymSetErrorDetails("Invalid irreducible represenation");
        goto err;
    }
    
    char types[] = {'A','B','E','T','G','H'}, *si[] = {"u","","g"}, *sv[] = {"2", "", "1"}, *sh[] = {"''", "", "'"};
    char type = irrep->d == 1 ? types[(1 - irrep->eig.p) >> 1] : types[irrep->d];

    if(irrep->d == 1){
        snprintf(name,l,"%c%s%s%s",type,sv[irrep->eig.v+1],si[irrep->eig.i+1],sh[irrep->eig.h+1]);
    }
    else if (irrep->eig.l >  0){
        snprintf(name,l,"%c%d%s%s",type,irrep->eig.l,si[irrep->eig.i+1],sh[irrep->eig.h+1]);
    } else {
        snprintf(name,l,"%c%s%s",type,si[irrep->eig.i+1],sh[irrep->eig.h+1]);
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
