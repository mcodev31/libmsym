//
//  basis_function.c
//  libmsym
//
//  Created by Marcus Johansson on 07/11/14.
//  Copyright (c) 2014 Marcus Johansson. 
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#include <stdlib.h>
#include <string.h>

#include "basis_function.h"

//These are real shperical harmonics cannot handle complex
msym_error_t basisFunctionFromQuantumNumbers(int n, int l, int m, msym_basis_function_t *bf){
    
    if(l > n || abs(m) > l) goto err;
    
    bf->f.rsh.n = n;
    bf->f.rsh.l = l;
    bf->f.rsh.m = m;
    
    memset(bf->name,0,sizeof(bf->name));
    
    switch(l) {
        case 0 :
            snprintf(bf->name, sizeof(bf->name), "%ds",n);
            break;
        case 1 : {
            char *d = "?";
            switch(m) {
                case -1 :
                    d = "y";
                    break;
                case 1 :
                    d = "x";
                    break;
                case 0 :
                    d = "z";
                    break;
            }
            snprintf(bf->name, sizeof(bf->name), "%dp%s",n,d);
            break;
        }
        case 2 : {
            char *d = (m < 0 ? "-" : "+");
            snprintf(bf->name, sizeof(bf->name), "%dd%d%s",n,abs(m),d);
            break;
        }
        default : {
            char t = l <= 20 ? ('f' - 3 + l + (l >= 7) + (l >= 12) + (l >= 14)) : '?';
            char *d = (m == 0 ? "" : m < 0 ? "-" : "+");
            snprintf(bf->name, sizeof(bf->name), "%d%c%d%s",n,t,abs(m),d);
        }
    }
    return MSYM_SUCCESS;
err:
    msymSetErrorDetails("Invalid orbital quantum numbers n:%d l:%d m:%d",n,l,m);
    return MSYM_INVALID_BASIS_FUNCTIONS;
}

msym_error_t basisFunctionFromName(char *name, msym_basis_function_t *bf){
    int n, l, m;
    char cl, cm1 = '\0', cm2 = '\0';
    
    sscanf(name,"%d%c%c%c",&n,&cl,&cm1,&cm2);
    
    switch(cl) {
        case 's' : l = m = 0; break;
        case 'p' : {
            l = 1;
            switch(cm1) {
                case 'x' : m = 1; break;
                case 'y' : m = -1; break;
                case 'z' : m = 0; break;
                default : goto err;
            }
            break;
        }
        default :
            if(cl < 'd' || cl == 'e' || cl == 'j' || cl > 'z') goto err;
            
            l = cl - 'd' + 2 - (cl > 'e') - (cl > 'j') - (cl > 'p') - (cl > 's');
            
            m = (cm1 - (int)'0')*(cm2 == '-' ? -1 : 1);
    }
    
    return basisFunctionFromQuantumNumbers(n,l,m,bf);
    
err:
    msymSetErrorDetails("Invalid orbital name %s",name);
    return MSYM_INVALID_BASIS_FUNCTIONS;
    
}
