//
//  debug.c
//  libmsym
//
//  Created by Marcus Johansson on 25/10/15.
//  Copyright (c) 2015 Marcus Johansson.
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#include "debug.h"

#ifdef LIBMSYM_DEBUG

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>



void tabPrintTransform(int r, int c, double M[r][c],int indent);
void tabprintf(char *format, int indent, ...);

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

void printSubspace(msym_character_table_t *ct, int l, msym_subrepresentation_space_t srs[l]){
    for(int k = 0;k < l;k++){
        printf("Subspace %d %s\n",k,ct->s[srs[k].s].name);
        for(int i = 0;i < srs[k].salcl;i++){
            for(int j = 0;j < srs[k].salc[i].fl;j++){
                msym_basis_function_t *bf = srs[k].salc[i].f[j];
                if(bf == NULL){
                    printf("error bf\n");
                    exit(1);
                }
                printf("\t  %s%s\t\t",bf->element->name,bf->name);
            }
            printf("\n");
            
            double (*space)[srs[k].salc[i].fl] = (double (*)[srs[k].salc[i].fl]) srs[k].salc[i].pf;
            if(space == NULL){
                printf("error space\n");
                exit(1);
            }
            tabPrintTransform(srs[k].salc[i].d,srs[k].salc[i].fl,space,1);
        }
    }
}

void printPermutation(msym_permutation_t *perm){
    int l = perm->p_length;
    printf("(");
    for(int j = 0; j < l; j++){
        printf(j == l -1 ? "%d" : "%d\t",j);
    }
    printf(")\n(");
    for(int j = 0; j < l; j++){
        printf(j == l -1 ? "%d" : "%d\t",perm->p[j]);
    }
    printf(")\n");
    
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

void printCharacterTable(msym_character_table_t *ct){
    msym_symmetry_operation_t **sops = ct->sops;
    int sopsl = ct->d;
    double (*table)[ct->d] = (double (*)[ct->d]) ct->table;
    printf("\t");
    for(int j = 0;j < sopsl;j++){
        char buf[12];
        symmetryOperationName(sops[j], 12, buf);
        printf("%d%s",ct->classc[sops[j]->cla],buf);
        printf("\t\t");
    }
    
    printf("\n");
    
    for(int i = 0;i < ct->d;i++){
        
        printf("%s\t",ct->s[i].name);
        for(int j = 0;j < ct->d;j++){
            printf("% .3lf\t\t",table[i][j]);
        }
        printf("\n");
    }
    
}
#endif
