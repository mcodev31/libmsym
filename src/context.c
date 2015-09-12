//
//  context.c
//  libmsym
//
//  Created by Marcus Johansson on 30/01/15.
//  Copyright (c) 2015 Marcus Johansson. 
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "msym.h"
#include "context.h"
#include "elements.h"
#include "equivalence_set.h"
#include "geometry.h"
#include "linalg.h"
#include "subspace.h"

msym_thresholds_t default_thresholds = {
    .zero = DEFAULT_ZERO_THRESHOLD,
    .geometry = DEFAULT_GEOMETRY_THRESHOLD,
    .angle = DEFAULT_ANGLE_THRESHOLD,
    .equivalence = DEFAULT_EQUIVALENCE_THRESHOLD,
    .eigfact = DEFAULT_EIGFACT_THRESHOLD,
    .permutation = DEFAULT_PERMUTATION_THRESHOLD,
    .orthogonalization = DEFAULT_ORTHOGONALIZATION_THRESHOLD
};


struct _msym_context {
    msym_thresholds_t *thresholds;
    msym_element_t *elements;
    msym_element_t **pelements;
    msym_orbital_t *orbitals;
    msym_orbital_t **porbitals;
    msym_basis_function_t *basis;
    msym_equivalence_set_t *es;
    msym_equivalence_set_t **eesmap;
    msym_permutation_t **es_perm;
    msym_subspace_t *oss;
    msym_subspace_t *dss;
    msym_subspace_2_t *salc_ss;
    int *oss_span;
    int *dss_span;
    int *salc_span;
    int el;
    int ol;
    int basisl;
    int esl;
    int ossl;
    int dssl;
    int salc_ssl;
    int es_perml;
    int sgl;
    msym_point_group_t *pg;
    msym_subgroup_t *sg;
    double cm[3];
    msym_geometry_t geometry;
    double eigval[3];
    double eigvec[3][3];
    struct _external_data {
        msym_element_t *elements;
        msym_orbital_t *orbitals;
        msym_orbital_t **porbitals;
        msym_subspace_2_t *salc_ss;
        msym_basis_function_t *basis;
        msym_symmetry_operation_t *sops;
        msym_subgroup_t *sg;
        msym_equivalence_set_t *es;
    } ext;
};

/***********************
 * Public API
 ***********************/

msym_context msymCreateContext(){
    msym_context ctx = malloc(sizeof(struct _msym_context));
    msym_thresholds_t *threshols = malloc(sizeof(msym_thresholds_t));
    
    //Don't generally handle allocation errors...
    if(ctx == NULL) {msymSetErrorDetails("Context memory allocation failed"); goto err;}
    if(threshols == NULL) {msymSetErrorDetails("Thresholds memory allocation failed"); goto err;}    
    
    memset(ctx, 0, sizeof(struct _msym_context));
    
    ctx->geometry = GEOMETRY_UNKNOWN;
    
    ctx->thresholds = threshols;
    msymSetThresholds(ctx, &default_thresholds);
    
    return ctx;
    
    err :
    free(ctx);
    free(threshols);
    return NULL;
}

msym_error_t msymSetThresholds(msym_context ctx, msym_thresholds_t *thresholds){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(thresholds->angle < 1.0 && !signbit(thresholds->angle) &&
       thresholds->equivalence < 1.0 && !signbit(thresholds->equivalence) &&
       thresholds->geometry < 1.0 && !signbit(thresholds->geometry) &&
       !signbit(thresholds->eigfact) &&
       !signbit(thresholds->orthogonalization) &&
       !signbit(thresholds->zero) &&
       !signbit(thresholds->permutation)){
        memcpy(ctx->thresholds, thresholds, sizeof(msym_thresholds_t));
    } else {
        ret = MSYM_INVALID_THRESHOLD;
    }
    
err:
    return ret;
}

msym_error_t msymGetThresholds(msym_context ctx, msym_thresholds_t **thresholds){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->thresholds == NULL) ret = MSYM_INVALID_THRESHOLD;
    *thresholds = ctx->thresholds;
err:
    return ret;
}


msym_error_t msymSetElements(msym_context ctx, int length, msym_element_t elements[length]){
    msym_error_t ret = MSYM_SUCCESS;
    msym_thresholds_t *thresholds = NULL;
    struct {msym_orbital_t *s; msym_orbital_t *e;} aorb = {.s= NULL, .e = NULL} ;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    ctxDestroyElements(ctx);
    
    if(MSYM_SUCCESS != (ret = msymGetThresholds(ctx, &thresholds))) goto err;
    
    ctx->elements = malloc(sizeof(msym_element_t[length]));
    ctx->pelements = malloc(sizeof(msym_element_t *[length]));
    
    for(int i = 0; i < length;i++){
        ctx->pelements[i] = &ctx->elements[i];
        memcpy(&ctx->elements[i], &elements[i], sizeof(msym_element_t));
        ctx->elements[i].ao = NULL;
        ctx->ol += elements[i].aol;
        for(msym_orbital_t **o = elements[i].ao;o < elements[i].ao + elements[i].aol;o++){
            if(aorb.s == NULL || *o < aorb.s) aorb.s = *o;
            if(aorb.e == NULL || *o >= aorb.e) aorb.e = *o+1;
        }
    }
    
    if(ctx->ol > 0){
        if(aorb.e - aorb.s != ctx->ol) { //Make sure we're pointing into a sequential array
            ret = MSYM_INVALID_ORBITALS;
            msymSetErrorDetails("Orbital data not in continuous memory block");
            goto err;
        }
        ctx->orbitals = malloc(ctx->ol*sizeof(msym_orbital_t));
        ctx->porbitals = malloc(ctx->ol*sizeof(msym_orbital_t*));
        msym_orbital_t **porb = ctx->porbitals;
        
        for(int i = 0; i < length;i++){
            ctx->elements[i].ao = porb;
            for(int j = 0;j < ctx->elements[i].aol;j++){
                ctx->elements[i].ao[j] = elements[i].ao[j] - aorb.s + ctx->orbitals;
                if(elements[i].ao[j]->n <= 0){
                    if(MSYM_SUCCESS != (ret = orbitalFromName(elements[i].ao[j]->name,ctx->elements[i].ao[j]))) goto err;
                } else {
                    if(MSYM_SUCCESS != (ret = orbitalFromQuantumNumbers(elements[i].ao[j]->n, elements[i].ao[j]->l, elements[i].ao[j]->m, ctx->elements[i].ao[j]))) goto err;
                }
            }
            porb += ctx->elements[i].aol;
        }
    }
    
    ctx->el = length;
    
    for(int i = 0;i < length;i++){
        if(MSYM_SUCCESS != (ret = complementElementData(&ctx->elements[i]))) goto err;
    }
    
    
    if(MSYM_SUCCESS != (ret = findCenterOfMass(ctx->el,ctx->pelements,ctx->cm))) goto err;
    
    for(msym_element_t *a = ctx->elements; a < (ctx->elements+length); a++){
        vsub(a->v,ctx->cm,a->v);
    }
    
    double zero[3] = {0,0,0};
    
    if(MSYM_SUCCESS != (ret = findGeometry(length, ctx->pelements, zero, thresholds, &ctx->geometry, ctx->eigval, ctx->eigvec))) goto err;
    
    return ret;
err:
    free(ctx->orbitals);
    free(ctx->porbitals);
    free(ctx->elements);
    free(ctx->pelements);
    ctx->elements = NULL;
    ctx->pelements = NULL;
    ctx->orbitals = NULL;
    ctx->porbitals = NULL;
    ctx->el = 0;
    ctx->ossl = 0;
    ctx->ol = 0;

    return ret;
}

/* We don't really need to copy this every time */
msym_error_t msymGetElements(msym_context ctx, int *length, msym_element_t **elements){
    msym_error_t ret = MSYM_SUCCESS;
    msym_element_t *relements = NULL;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->elements == NULL) {ret = MSYM_INVALID_ELEMENTS;goto err;}
    if(ctx->ext.elements == NULL) ctx->ext.elements = malloc(sizeof(msym_element_t[ctx->el]));
    if(ctx->orbitals != NULL) {
        if(ctx->ext.orbitals == NULL) ctx->ext.orbitals = malloc(sizeof(msym_orbital_t[ctx->ol]));
        memcpy(ctx->ext.orbitals,ctx->orbitals,sizeof(msym_orbital_t[ctx->ol]));
    }
    if(ctx->porbitals != NULL){
        if(ctx->ext.porbitals == NULL) ctx->ext.orbitals = calloc(ctx->ol,sizeof(msym_orbital_t*));
    }
    
    memcpy(ctx->ext.elements,ctx->elements,sizeof(msym_element_t[ctx->el]));
    msym_orbital_t **porb = ctx->ext.porbitals;
    for(msym_element_t *a = ctx->ext.elements; a < (ctx->ext.elements+ctx->el); a++){
        vadd(a->v,ctx->cm,a->v);
        for(int i = 0;i < a->aol && ctx->ext.orbitals != NULL && porb != NULL;i++){
            porb[i] = a->ao[i] - ctx->orbitals + ctx->ext.orbitals;
        }
        if(porb != NULL){
            a->ao = porb;
            porb += a->aol;
        }
    }
    
    *elements = ctx->ext.elements;
    *length = ctx->el;
    return ret;
err:
    free(relements);
    *elements = NULL;
    *length = 0;
    return ret;
}

msym_error_t msymGetEquivalenceSets(msym_context ctx, int *length, msym_equivalence_set_t **es){
    msym_error_t ret = MSYM_SUCCESS;
    msym_element_t *elements = NULL;
    int el = 0;
    
    *es = NULL;
    
    if(ctx->es == NULL){ret = MSYM_INVALID_EQUIVALENCE_SET;goto err;}
    if(ctx->ext.es == NULL){
        if(MSYM_SUCCESS != (ret = msymGetElements(ctx, &el, &elements))) goto err; //a bit lazy
        if(MSYM_SUCCESS != (ret = copyEquivalenceSets(ctx->esl,ctx->es,&ctx->ext.es))) goto err;
        for(int i = 0; i < ctx->esl;i++){
            for(int j = 0;j < ctx->es[i].length;j++){
                ctx->ext.es[i].elements[j] = ctx->ext.es[i].elements[j] - ctx->elements + elements;
            }
        }
    }
    
    *es = ctx->ext.es;
    *length = ctx->esl;
    
err:
    return ret;
}

msym_error_t msymGetBasisFunctions(msym_context ctx, int *length, msym_basis_function_t **basis){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->basis == NULL) {ret = MSYM_INVALID_ORBITALS;goto err;}
    if(ctx->ext.basis == NULL) ctx->ext.basis = malloc(sizeof(msym_basis_function_t[ctx->basisl]));
    memcpy(ctx->ext.basis,ctx->basis,sizeof(msym_basis_function_t[ctx->basisl]));
    for(int i = 0;i < ctx->basisl;i++) ctx->ext.basis[i].element = ctx->basis[i].element - ctx->elements + ctx->ext.elements;
    
    *length = ctx->basisl;
    *basis = ctx->ext.basis;
err:
    return ret;
}

msym_error_t msymSetBasisFunctions(msym_context ctx, int length, msym_basis_function_t *basis){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->elements == NULL) {ret = MSYM_INVALID_ELEMENTS;goto err;}
    ctxDestroyBasisFunctions(ctx);
    ctx->basis = malloc(sizeof(msym_basis_function_t[length]));
    memcpy(ctx->basis, basis, sizeof(msym_basis_function_t[length]));
    for(int i = 0;i < length;i++){
        msym_basis_function_t *bf = &ctx->basis[i];
        if(bf->element >= ctx->ext.elements && ctx->basis[i].element < ctx->ext.elements + ctx->el){
            bf->element = bf->element - ctx->ext.elements + ctx->elements;
        } else {
            ret = MSYM_INVALID_ORBITALS;
            goto err;
        }
        
        if(bf->type != REAL_SPHERICAL_HARMONIC){
            ret = MSYM_INVALID_ORBITALS;
            goto err;
        }
        
        msym_spherical_harmonic_t *sh = &bf->f.sh;
        
        if(sh->n <= 0){
            if(MSYM_SUCCESS != (ret = basisFunctionFromName(bf->name,bf))) goto err;
        } else {
            if(MSYM_SUCCESS != (ret = basisFunctionFromQuantumNumbers(sh->n, sh->l, sh->m, bf))) goto err;
        }
    }
    
    ctx->basisl = length;
    
    return ret;

err:
    free(ctx->basis);
    ctx->basis = NULL;
    return ret;
}

msym_error_t msymGetPointGroupName(msym_context ctx, int l, char buf[l]){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->pg == NULL) {ret = MSYM_INVALID_POINT_GROUP;goto err;}
    snprintf(buf, l, "%s",ctx->pg->name);
err:
    return ret;
}

msym_error_t msymGetSubgroups(msym_context ctx, int *sgl, msym_subgroup_t **sg){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->pg == NULL) {ret = MSYM_INVALID_POINT_GROUP;goto err;}
    if(ctx->pg->perm == NULL) {ret = MSYM_INVALID_PERMUTATION;goto err;}
    
    if(ctx->ext.sops == NULL){
        msym_symmetry_operation_t *extsops = NULL;
        int extsopsl = 0;
        if(MSYM_SUCCESS != (ret = msymGetSymmetryOperations(ctx, &extsopsl, &extsops))) goto err;
    }
    if(ctx->sg == NULL){
        int sgmax = numberOfSubgroups(ctx->pg);
        if(MSYM_SUCCESS != (ret = findPermutationSubgroups(ctx->pg->order, ctx->pg->perm, sgmax, ctx->pg->sops, &ctx->sgl, &ctx->sg))) goto err;

        for(int i = 0;i < ctx->sgl;i++){
            if(MSYM_SUCCESS != (ret = findSubgroup(&ctx->sg[i], ctx->thresholds))) goto err;
        }
    }
    
    if(ctx->ext.sg == NULL){
        ctx->ext.sg = malloc(sizeof(msym_subgroup_t[ctx->sgl]));
        memcpy(ctx->ext.sg, ctx->sg, sizeof(msym_subgroup_t[ctx->sgl]));
        for(int i = 0;i < ctx->sgl;i++){
            ctx->ext.sg[i].sops = malloc(sizeof(msym_symmetry_operation_t *[ctx->sg[i].sopsl]));
            for(int j = 0;j < ctx->sg[i].sopsl;j++){
                ctx->ext.sg[i].sops[j] = ctx->sg[i].sops[j] - ctx->pg->sops + ctx->ext.sops;
                ctx->ext.sg[i].subgroup[0] = ctx->sg[i].subgroup[0] == NULL ? NULL : ctx->sg[i].subgroup[0] - ctx->sg + ctx->ext.sg;
                ctx->ext.sg[i].subgroup[1] = ctx->sg[i].subgroup[1] == NULL ? NULL : ctx->sg[i].subgroup[1] - ctx->sg + ctx->ext.sg;
            }
        }
    }
    

    *sgl = ctx->sgl;
    *sg = ctx->ext.sg;
    return ret;
    
err:
    return ret;
    
    
    
}

msym_error_t msymGetSALCSubspaces(msym_context ctx, int *l, msym_subspace_2_t **ss){
    msym_error_t ret = MSYM_SUCCESS;
    if(NULL == ctx) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(NULL == ctx->salc_ss){
        if(MSYM_SUCCESS != (ret = msymGenerateSALCSubspaces(ctx))) goto err;
        if(NULL == ctx->salc_ss){ret = MSYM_INVALID_ORBITALS;goto err;}
    }
    
    msym_basis_function_t *bfs;
    int bfsl = 0;
    if(MSYM_SUCCESS != (ret = msymGetBasisFunctions(ctx, &bfsl, &bfs))) goto err;
    msym_subspace_2_t *salc_ss = ctx->ext.salc_ss = calloc(ctx->salc_ssl, sizeof(msym_subspace_2_t));
    memcpy(ctx->ext.salc_ss, ctx->salc_ss, sizeof(msym_subspace_2_t[ctx->salc_ssl]));
    for(int i = 0; i < ctx->salc_ssl;i++){
        salc_ss[i].salcl = ctx->salc_ss[i].salcl;
        salc_ss[i].s = ctx->salc_ss[i].s;
        salc_ss[i].salc = malloc(sizeof(msym_salc_t[salc_ss[i].salcl]));
        for(int j = 0; j < salc_ss[i].salcl;j++){
            salc_ss[i].salc[j].d = ctx->salc_ss[i].salc[j].d;
            salc_ss[i].salc[j].fl = ctx->salc_ss[i].salc[j].fl;
            size_t pfsize = sizeof(double[salc_ss[i].salc[j].d][salc_ss[i].salc[j].fl]);
            salc_ss[i].salc[j].pf = malloc(pfsize);
            salc_ss[i].salc[j].f = malloc(sizeof(msym_basis_function_t *[salc_ss[i].salc[j].fl]));
            memcpy(salc_ss[i].salc[j].pf, ctx->salc_ss[i].salc[j].pf, pfsize);
            for(int k = 0;k < salc_ss[i].salc[j].fl;k++){
                salc_ss[i].salc[j].f[k] = ctx->salc_ss[i].salc[j].f[k] - ctx->basis + ctx->ext.basis;
            }
        }
    }
    
    *ss = ctx->ext.salc_ss;
    *l = ctx->salc_ssl;
    
    
    return ret;
err:
    freeSALCSubspaces(ctx->salc_ssl, ctx->ext.salc_ss);
    ctx->ext.salc_ss = NULL;
    return ret;
}

msym_error_t msymGetPointGroup(msym_context ctx, msym_point_group_t **pg){
    msym_error_t ret = MSYM_SUCCESS;
    if(NULL == ctx) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(NULL == ctx->pg){ret = MSYM_INVALID_POINT_GROUP;goto err;}
    printf("WARNING returning internal pointer\n");
    *pg = ctx->pg;
    
err:
    return ret;
    
}

msym_error_t msymGetCharacterTable(msym_context ctx, msym_character_table_t **ct){
    msym_error_t ret = MSYM_SUCCESS;
    if(NULL == ctx) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(NULL == ctx->pg){ret = MSYM_INVALID_POINT_GROUP;goto err;}
    if(NULL == ctx->pg->ct2){
        msym_point_group_t *pg = ctx->pg;
        if(MSYM_SUCCESS != (ret = generateCharacterTable(pg->type, pg->n, pg->order, pg->sops, &pg->ct2))) goto err;
    }
    printf("WARNING returning internal pointer\n");
    *ct = ctx->pg->ct2;
    
err:
    return ret;
    
}

msym_error_t msymGetCenterOfMass(msym_context ctx, double v[3]){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->elements == NULL) {ret = MSYM_INVALID_ELEMENTS;goto err;}
    vcopy(ctx->cm, v);
err:
    return ret;
}

msym_error_t msymGetGeometry(msym_context ctx, msym_geometry_t *geometry){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->elements == NULL) {ret = MSYM_INVALID_ELEMENTS;goto err;}
    if(ctx->geometry == GEOMETRY_UNKNOWN) {ret = MSYM_INVALID_GEOMETRY;goto err;}
    *geometry = ctx->geometry;
err:
    return ret;
}

msym_error_t msymGetPrincipalMoments(msym_context ctx, double eigval[3]){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->elements == NULL) {ret = MSYM_INVALID_ELEMENTS;goto err;}
    vcopy(ctx->eigval, eigval);
err:
    return ret;
}
msym_error_t msymGetPrincipalAxes(msym_context ctx, double eigvec[3][3]){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->elements == NULL) {ret = MSYM_INVALID_ELEMENTS;goto err;}
    mcopy(ctx->eigvec, eigvec);
err:
    return ret;
}

msym_error_t msymGetRadius(msym_context ctx, double *radius){
    msym_error_t ret = MSYM_SUCCESS;
    double r = 0.0;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->elements == NULL) {ret = MSYM_INVALID_ELEMENTS;goto err;}
    for(int i = 0;i < ctx->el;i++){
        double abs = vabs(ctx->elements[i].v);
        r = r > abs ? r : abs;
    }
    *radius = r;

err:
    return ret;
}

msym_error_t msymGetSymmetryOperations(msym_context ctx, int *sopsl, msym_symmetry_operation_t **sops){
    msym_error_t ret = MSYM_SUCCESS;
    msym_symmetry_operation_t *rsops = NULL;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->pg == NULL || ctx->pg->sops == NULL) {ret = MSYM_INVALID_POINT_GROUP;goto err;}
    if(ctx->ext.sops == NULL) ctx->ext.sops = malloc(sizeof(msym_symmetry_operation_t[ctx->pg->order]));
    memcpy(ctx->ext.sops,ctx->pg->sops,sizeof(msym_symmetry_operation_t[ctx->pg->order]));
    *sops = ctx->ext.sops;
    *sopsl = ctx->pg->order;
    return ret;
err:
    free(rsops);
    *sops = NULL;
    *sopsl = 0;
    return ret;
}





//Fix this only for testing
msym_error_t msymSetOrbitalsMB(msym_context ctx){
    int cnt = 0;
    for(int i = 0;i < ctx->el;i++){
        cnt++;
        if(ctx->elements[i].n > 2){
            cnt += 4;
        }
    }
    ctx->orbitals = malloc(cnt*sizeof(msym_orbital_t));
    msym_orbital_t **porb = malloc(cnt*sizeof(msym_orbital_t*));
    
    ctx->ol = cnt;
    int cnt2 = 0;
    for(int i = 0;i < ctx->el;i++){
        //double v[3];
        ctx->elements[i].ao = porb + cnt2;
        ctx->elements[i].aol = 1;
        //vcopy(ctx->elements[i].v, v);
        
        ctx->elements[i].ao[0] = (msym_orbital_t*)&ctx->orbitals[cnt2];
        orbitalFromQuantumNumbers(1, 0, 0, (msym_orbital_t*)&ctx->orbitals[cnt2++]);
        
        if(ctx->elements[i].n > 2){
            ctx->elements[i].aol += 4;
            
            ctx->elements[i].ao[1] = (msym_orbital_t*)&ctx->orbitals[cnt2];
            orbitalFromQuantumNumbers(2, 0,  0, (msym_orbital_t*)&ctx->orbitals[cnt2++]);
            
            ctx->elements[i].ao[2] = (msym_orbital_t*)&ctx->orbitals[cnt2];
            orbitalFromQuantumNumbers(2, 1, -1, (msym_orbital_t*)&ctx->orbitals[cnt2++]);
            
            ctx->elements[i].ao[3] = (msym_orbital_t*)&ctx->orbitals[cnt2];
            orbitalFromQuantumNumbers(2, 1,  0, (msym_orbital_t*)&ctx->orbitals[cnt2++]);
            
            ctx->elements[i].ao[4] = (msym_orbital_t*)&ctx->orbitals[cnt2];
            orbitalFromQuantumNumbers(2, 1,  1, (msym_orbital_t*)&ctx->orbitals[cnt2++]);
        }
    }
    printf("generated %d orbitals on %d elements\n",ctx->ol, ctx->el);
    ctx->porbitals = porb;
    return MSYM_SUCCESS;
}

msym_error_t msymReleaseContext(msym_context ctx){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    free(ctx->thresholds);
    ctxDestroyElements(ctx);
    ctxDestroyPointGroup(ctx);
    free(ctx);
err:
    return ret;
}




/***********************
 * Private API
 ***********************/


msym_error_t ctxGetElements(msym_context ctx, int *l, msym_element_t **elements){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->elements == NULL) {ret = MSYM_INVALID_ELEMENTS; goto err;}
    *elements = (msym_element_t *) ctx->elements;
    *l = ctx->el;
err:
    return ret;
}

msym_error_t ctxGetInternalElement(msym_context ctx, msym_element_t *ext, msym_element_t **element){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->ext.elements == NULL) {ret = MSYM_INVALID_ELEMENTS;goto err;}
    if(ext < ctx->ext.elements || ext >= ctx->ext.elements+ctx->el){
        msymSetErrorDetails("Element pointer (%p) outside memory block (%p -> %p)", ext, ctx->ext.elements, ctx->ext.elements + ctx->el);
        ret = MSYM_INVALID_ELEMENTS;
        goto err;
    }
    *element = ext - ctx->ext.elements + ctx->elements;
err:
    return ret;
}

msym_error_t ctxGetSubgroups(msym_context ctx, int *sgl, msym_subgroup_t **sg){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->sg == NULL) {ret = MSYM_INVALID_SUBGROUPS;goto err;}
    *sg = ctx->sg;
    *sgl = ctx->sgl;
err:
    return ret;
}

msym_error_t ctxSetSubgroups(msym_context ctx, int sgl, msym_subgroup_t *sg){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    ctxDestroySubgroups(ctx);
    ctx->sg = sg;
    ctx->sgl = sgl;
err:
    return ret;
}

msym_error_t ctxGetInternalSubgroup(msym_context ctx, msym_subgroup_t *ext, msym_subgroup_t **sg){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->ext.sg == NULL) {ret = MSYM_INVALID_SUBGROUPS;goto err;}
    if(ext < ctx->ext.sg || ext >= ctx->ext.sg+ctx->sgl){
        msymSetErrorDetails("Subgroup pointer (%p) outside memory block (%p -> %p)", ext, ctx->ext.sg, ctx->ext.sg + ctx->sgl);
        ret = MSYM_INVALID_POINT_GROUP;
        goto err;
    }
    *sg = ext - ctx->ext.sg + ctx->sg;
err:
    return ret;
}

msym_error_t ctxGetElementPtrs(msym_context ctx, int *l, msym_element_t ***pelements){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->pelements == NULL) {ret = MSYM_INVALID_ELEMENTS; goto err;}
    *pelements = (msym_element_t **) ctx->pelements;
    *l = ctx->el;
err:
    return ret;
}

msym_error_t ctxGetOrbitals(msym_context ctx, int *l, msym_orbital_t **orbitals){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    if(ctx->orbitals == NULL) {ret = MSYM_INVALID_ORBITALS; goto err;}
    *orbitals = (msym_orbital_t *) ctx->orbitals;
    *l = ctx->ol;
err:
    return ret;
}

msym_error_t ctxGetBasisFunctions(msym_context ctx, int *l, msym_basis_function_t **basis){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    if(ctx->basis == NULL) {ret = MSYM_INVALID_ORBITALS; goto err;}
    *basis = (msym_basis_function_t *) ctx->basis;
    *l = ctx->basisl;
err:
    return ret;
}

msym_error_t ctxSetCenterOfMass(msym_context ctx, double cm[3]){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    vcopy(cm,ctx->cm);
err:
    return ret;
}

msym_error_t ctxSetPointGroup(msym_context ctx, msym_point_group_t *pg){
    msym_error_t ret = MSYM_SUCCESS;
    if(MSYM_SUCCESS != (ret = ctxDestroyPointGroup(ctx))) goto err;
    ctx->pg = pg;
err:
    return ret;
}

msym_error_t ctxGetPointGroup(msym_context ctx, msym_point_group_t **pg){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    if(ctx->pg == NULL) {ret = MSYM_INVALID_POINT_GROUP; goto err;}
    *pg = ctx->pg;
err:
    return ret;
}

msym_error_t ctxSetEquivalenceSets(msym_context ctx, int esl, msym_equivalence_set_t *es){
    msym_error_t ret = MSYM_SUCCESS;
    if(MSYM_SUCCESS != (ret = ctxDestroyEquivalcenceSets(ctx))) goto err;
    if(!ctx->eesmap) ctx->eesmap = calloc(ctx->el, sizeof(msym_equivalence_set_t *));
    for(int i = 0;i < esl;i++){
        for(int j = 0;j < es[i].length;j++){
            ctx->eesmap[es[i].elements[j] - ctx->elements] = &es[i];
        }
    }
    for(int i = 0; i < ctx->el;i++){
        if(!ctx->eesmap[i]){
            ret = MSYM_INVALID_EQUIVALENCE_SET;
            msymSetErrorDetails("Element %d does not map to any equivalence set",i);
            goto err;
        }
    }
    ctx->es = es;
    ctx->esl = esl;
err:
    return ret;
}

msym_error_t ctxGetEquivalenceSets(msym_context ctx, int *esl, msym_equivalence_set_t **es){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    if(ctx->es == NULL) {ret = MSYM_INVALID_EQUIVALENCE_SET;goto err;}
    *es = ctx->es;
    *esl = ctx->esl;
err:
    return ret;
}

msym_error_t ctxGetElementEquivalenceSetMap(msym_context ctx, msym_equivalence_set_t ***eesmap){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    if(ctx->eesmap == NULL) {ret = MSYM_INVALID_EQUIVALENCE_SET;goto err;}
    *eesmap = ctx->eesmap;
err:
    return ret;
}

msym_error_t ctxSetEquivalenceSetPermutations(msym_context ctx, int r, int c, msym_permutation_t **perm){
    msym_error_t ret = MSYM_SUCCESS;
    if(MSYM_SUCCESS != (ret = ctxDestroyEquivalcenceSetPermutations(ctx))) goto err;
    if(r != ctx->esl || ctx->pg == NULL || c != ctx->pg->order){ret = MSYM_INVALID_PERMUTATION; goto err;}
    ctx->es_perm = perm;
    ctx->es_perml = c;
err:
    return ret;
}

msym_error_t ctxGetEquivalenceSetPermutations(msym_context ctx, int *r, int *c, msym_permutation_t ***perm){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    if(ctx->es == NULL || ctx->es_perml == 0 || ctx->es_perm == NULL){ret = MSYM_INVALID_PERMUTATION;goto err;}
    *r = ctx->esl;
    *c = ctx->es_perml;
    *perm = ctx->es_perm;
err:
    return ret;
}

msym_error_t ctxGetOrbitalSubspaces(msym_context ctx, int *ssl, msym_subspace_t **ss, int **span){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    if(ctx->oss == NULL) {ret = MSYM_INVALID_SUBSPACE; goto err;}
    *ssl = ctx->ossl;
    *ss = ctx->oss;
    *span = ctx->oss_span;
err:
    return ret;
}

msym_error_t ctxSetOrbitalSubspaces(msym_context ctx, int ssl, msym_subspace_t *ss, int *span){
    msym_error_t ret = MSYM_SUCCESS;
    if(MSYM_SUCCESS != (ret = ctxDestroyOrbitalSubspaces(ctx))) goto err;
    ctx->ossl = ssl;
    ctx->oss = ss;
    ctx->oss_span = span;
err:
    return ret;
}

msym_error_t ctxGetSALCSubspaces(msym_context ctx, int *ssl, msym_subspace_2_t **ss, int **span){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    if(ctx->salc_ss == NULL) {ret = MSYM_INVALID_SUBSPACE; goto err;}
    *ssl = ctx->salc_ssl;
    *ss = ctx->salc_ss;
    *span = ctx->salc_span;
err:
    return ret;
}

msym_error_t ctxSetSALCSubspaces(msym_context ctx, int ssl, msym_subspace_2_t *ss, int *span){
    msym_error_t ret = MSYM_SUCCESS;
    if(MSYM_SUCCESS != (ret = ctxDestroySALCSubspaces(ctx))) goto err;
    ctx->salc_ssl = ssl;
    ctx->salc_ss = ss;
    ctx->salc_span = span;
err:
    return ret;
}

msym_error_t ctxGetDisplacementSubspaces(msym_context ctx, int *ssl, msym_subspace_t **ss, int **span){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    if(ctx->dss == NULL) {ret = MSYM_INVALID_SUBSPACE; goto err;}
    *ssl = ctx->dssl;
    *ss = ctx->dss;
    *span = ctx->dss_span;
err:
    return ret;
}

msym_error_t ctxSetDisplacementSubspaces(msym_context ctx, int ssl, msym_subspace_t *ss, int *span){
    msym_error_t ret = MSYM_SUCCESS;
    if(MSYM_SUCCESS != (ret = ctxDestroyDisplacementSubspaces(ctx))) goto err;
    ctx->dssl = ssl;
    ctx->dss = ss;
    ctx->dss_span = span;
err:
    return ret;
}

msym_error_t ctxGetGeometry(msym_context ctx, msym_geometry_t *g, double eigval[3], double eigvec[3][3]){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    if(ctx->geometry == GEOMETRY_UNKNOWN) {ret = MSYM_INVALID_GEOMETRY; goto err;}
    *g = ctx->geometry;
    mcopy(ctx->eigvec, eigvec);
    vcopy(ctx->eigval, eigval);
err:
    return ret;
}


msym_error_t ctxDestroyElements(msym_context ctx){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    ctxDestroyEquivalcenceSets(ctx);
    ctxDestroyOrbitalSubspaces(ctx);
    ctxDestroyDisplacementSubspaces(ctx);
    ctxDestroySALCSubspaces(ctx);
    ctxDestroyBasisFunctions(ctx);
    free(ctx->elements);
    free(ctx->pelements);
    free(ctx->orbitals);
    free(ctx->porbitals);
    free(ctx->eesmap);
    free(ctx->ext.elements);
    free(ctx->ext.orbitals);
    free(ctx->ext.porbitals);
    
    ctx->elements = NULL;
    ctx->pelements = NULL;
    ctx->orbitals = NULL;
    ctx->porbitals = NULL;
    ctx->eesmap = NULL;
    ctx->ext.elements = NULL;
    ctx->ext.orbitals = NULL;
    ctx->ext.porbitals = NULL;
    ctx->el = 0;
    ctx->ol = 0;
    ctx->geometry = GEOMETRY_UNKNOWN;
    memset(ctx->eigvec,0,sizeof(ctx->eigvec));
    memset(ctx->eigval,0,sizeof(ctx->eigval));
    memset(ctx->cm,0,sizeof(ctx->cm));
err:
    return ret;
}

msym_error_t ctxDestroyEquivalcenceSets(msym_context ctx){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    ctxDestroyEquivalcenceSetPermutations(ctx);
    if(ctx->eesmap) memset(ctx->eesmap, ctx->el, sizeof(msym_equivalence_set_t *));
    free(ctx->es);
    free(ctx->ext.es);
    ctx->ext.es = NULL;
    ctx->es = NULL;
    ctx->esl = 0;
err:
    return ret;
}

msym_error_t ctxDestroyEquivalcenceSetPermutations(msym_context ctx){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    for(int i = 0; i < ctx->esl;i++){
        for(int j = 0;j< ctx->es_perml;j++){
            freePermutationData(&ctx->es_perm[i][j]);
        }
    }
    free(ctx->es_perm);
    ctx->es_perm = NULL;
    ctx->es_perml = 0;
err:
    return ret;
}

msym_error_t ctxDestroyPointGroup(msym_context ctx){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    if(ctx->pg == NULL) goto err;
    ctxDestroyEquivalcenceSets(ctx);
    ctxDestroySubgroups(ctx);
    for(int i = 0;i < ctx->pg->order && ctx->pg->perm != NULL;i++){
        freePermutationData(&ctx->pg->perm[i]);
    }
    for(int i = 0;i < ctx->sgl && ctx->sg != NULL;i++){
        free(ctx->sg[i].sops);
    }
    for(int i = 0;i < ctx->sgl && ctx->ext.sg != NULL;i++){
        free(ctx->ext.sg[i].sops);
    }
    free(ctx->pg->perm);
    free(ctx->pg->ct);
    free(ctx->pg->sops);
    free(ctx->pg);
    free(ctx->ext.sops);
    
    
    ctx->pg = NULL;
    ctx->sg = NULL;
    ctx->ext.sops = NULL;
    ctx->ext.sg = NULL;
err:
    return ret;
}

msym_error_t ctxDestroySubgroups(msym_context ctx){
    msym_error_t ret = MSYM_SUCCESS;
    free(ctx->sg);
    free(ctx->ext.sg);
    ctx->sg = NULL;
    ctx->ext.sg = NULL;
    ctx->sgl = 0;
err:
    return ret;
}

msym_error_t ctxDestroyBasisFunctions(msym_context ctx){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    ctxDestroySALCSubspaces(ctx);
    free(ctx->basis);
    free(ctx->ext.basis);
    ctx->basis = NULL;
    ctx->ext.basis = NULL;
    ctx->basisl = 0;
err:
    return ret;
}

msym_error_t ctxDestroyOrbitalSubspaces(msym_context ctx){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    for(int i = 0;i < ctx->ossl && ctx->oss != NULL;i++){
        freeSubspace(&ctx->oss[i]);
    }
    free(ctx->oss);
    free(ctx->oss_span);
    ctx->oss_span = NULL;
    ctx->oss = NULL;
    ctx->ossl = 0;
err:
    return ret;
}

msym_error_t ctxDestroyDisplacementSubspaces(msym_context ctx){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    for(int i = 0;i < ctx->dssl && ctx->dss != NULL;i++){
        freeSubspace(&ctx->dss[i]);
    }
    free(ctx->dss);
    free(ctx->dss_span);
    ctx->dss_span = NULL;
    ctx->dss = NULL;
    ctx->dssl = 0;
err:
    return ret;
}

msym_error_t ctxDestroySALCSubspaces(msym_context ctx){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    freeSALCSubspaces(ctx->salc_ssl, ctx->salc_ss);
    freeSALCSubspaces(ctx->salc_ssl, ctx->ext.salc_ss);
    free(ctx->salc_span);
    ctx->salc_ss = NULL;
    ctx->ext.salc_ss = NULL;
    ctx->salc_span = NULL;
    ctx->salc_ssl = 0;
err:
    return ret;
}
