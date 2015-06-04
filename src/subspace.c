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
#include <stdarg.h>

#include "msym.h"
#include "subspace.h"

void tabprintf(char *format, int indent, ...);
void printSubspaceTree(CharacterTable *ct, msym_subspace_t *ss,int indent);
void tabPrintTransform(int r, int c, double M[r][c],int indent);

//Should filter so we don't get subspaces with just one subspace as well
int filterSubspace(msym_subspace_t *ss){
    int ret = 0;
    if(ss->subspacel == 0){
        ret = ss->d > 0 && ss->basisl > 0;
    } else {
        for(int i = 0;i < ss->subspacel;i++){
            if(!filterSubspace(&ss->subspace[i])){
                ss->subspacel--;
                
                if(ss->subspacel == 0){
                    free(ss->subspace);
                    ss->subspace = NULL;
                    break;
                } else {
                    memcpy(&ss->subspace[i], &ss->subspace[ss->subspacel], sizeof(msym_subspace_t));
                    ss->subspace = realloc(ss->subspace, sizeof(msym_subspace_t[ss->subspacel]));
                    i--;
                }
                
            }
        }
        ret = ss->subspacel > 0;
    }
    return ret;
}


void freeSubspace(msym_subspace_t *ss){
    free(ss->basis.o);
    free(ss->space);
    for(int i = 0; i < ss->subspacel;i++){
        freeSubspace(&ss->subspace[i]);
    }
    
    free(ss->subspace);
}

void printSubspace(CharacterTable *ct, msym_subspace_t *ss){
    printSubspaceTree(ct,ss,0);
}



void printSubspaceTree(CharacterTable *ct, msym_subspace_t *ss,int indent){
    if(ct == NULL){
        tabprintf("Subspace irrep: %d\n", indent,ss->irrep);
    } else {
        tabprintf("Subspace irrep: %s\n", indent,ct->irrep[ss->irrep].name);
    }
    if(ss->subspacel == 0){
        if(ss->d > 0 && ss->basisl > 0){
            tabprintf("", indent);
            for(int i = 0;i < ss->basisl;i++) {
                char *namev[3] = {"x","y","z"};
                if(ss->type == ATOMIC_ORBITAL){
                    printf("  %s\t",ss->basis.o[i]->name);
                } else if (ss->type == MASS_WEIGHTED_COORDINATES){
                    char *name = namev[(int)round(ss->basis.q[i].v[1]+2*ss->basis.q[i].v[2])];
                    printf("  %s%s\t",ss->basis.q[i].element->name,name);
                }
            }
            printf("\n");
            double (*space)[ss->basisl] = (double (*)[ss->basisl]) ss->space;
            tabPrintTransform(ss->d,ss->basisl,space,indent);
        } else {
            tabprintf("No subspaces spaned\n", indent);
        }
    } else {
        for(int i = 0; i < ss->subspacel;i++){
            printSubspaceTree(ct,&ss->subspace[i],indent+1);
        }
    }
}

void tabprintf(char *format, int indent, ...){
    for(int i = 0; i < indent;i++) printf("\t");
    va_list args;
    va_start (args, indent);
    vprintf (format, args);
    va_end (args);
}

void tabPrintTransform(int r, int c, double M[r][c],int indent) {
    if(r == 0 || c == 0) {tabprintf("[]\n",indent);return;}
    //printf("\n");
    tabprintf("[",indent);
    for(int i = 0;i < r;i++){
        for(int j = 0;j<c;j++){
            char *pre = signbit(M[i][j]) ? "" : " ";
            char *post1 = "\t";
            char *post2 = (j == (c - 1)) ? (i == (r - 1)) ? "" : ";" : " ";
            
            printf("%s%.8lf%s%s",pre,M[i][j],post1,post2);
        }
        printf("%s",(i == (r - 1)) ? "]\n" : "\n ");
        tabprintf(" ", indent);
    }
    printf("\n");
    
}

void printTransform(int r, int c, double M[r][c]) {
    printf("\n[");
    for(int i = 0;i < r;i++){
        for(int j = 0;j<c;j++){
            char *pre = signbit(M[i][j]) ? "" : " ";
            char *post1 = "";
            char *post2 = (j == (c - 1)) ? (i == (r - 1)) ? "" : ";" : " ";
            
            printf("%s%.8lf%s%s",pre,M[i][j],post1,post2);
        }
        printf("%s",(i == (r - 1)) ? "]\n" : "\n ");
    }
    
}