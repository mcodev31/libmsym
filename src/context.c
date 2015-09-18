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

const msym_thresholds_t default_thresholds = {
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
    msym_basis_function_t *basis;
    msym_equivalence_set_t *es;
    msym_permutation_t **es_perm;
    msym_subspace_t *salc_ss;
    int *oss_span;
    int *dss_span;
    int *salc_span;
    unsigned long int flags;
    int elementsl;
    int basisl;
    int esl;
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
        msym_equivalence_set_t **eesmap;
        msym_element_t *set_elements_ptr;
        msym_element_t *elements;
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
    
    ctx->geometry = MSYM_GEOMETRY_UNKNOWN;
    
    ctx->thresholds = threshols;
    msymSetThresholds(ctx, &default_thresholds);
    
    return ctx;
    
    err :
    free(ctx);
    free(threshols);
    return NULL;
}

const msym_thresholds_t *msymGetDefaultThresholds(){
    return &default_thresholds;
}

msym_error_t msymSetThresholds(msym_context ctx, const msym_thresholds_t *thresholds){
    msym_error_t ret = MSYM_SUCCESS;
    if(NULL == ctx) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(NULL != thresholds &&
       thresholds->angle < 1.0 && !signbit(thresholds->angle) &&
       thresholds->equivalence < 1.0 && !signbit(thresholds->equivalence) &&
       thresholds->geometry < 1.0 && !signbit(thresholds->geometry) &&
       !signbit(thresholds->eigfact) &&
       !signbit(thresholds->orthogonalization) &&
       !signbit(thresholds->zero) &&
       !signbit(thresholds->permutation)){
        if(ctx->thresholds != thresholds) memcpy(ctx->thresholds, thresholds, sizeof(msym_thresholds_t));
    } else {
        ret = MSYM_INVALID_THRESHOLD;
        goto err;
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
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}

    ctxDestroyElements(ctx);
    
    if(MSYM_SUCCESS != (ret = ctxGetThresholds(ctx, &thresholds))) goto err;
    
    ctx->elements = malloc(sizeof(msym_element_t[length]));
    ctx->pelements = malloc(sizeof(msym_element_t *[length]));
    memcpy(ctx->elements, elements, sizeof(msym_element_t[length]));
    ctx->elementsl = length;
    
    for(int i = 0;i < length;i++){
        ctx->pelements[i] = &(ctx->elements[i]);
        if(MSYM_SUCCESS != (ret = complementElementData(ctx->pelements[i]))) goto err;
        printElement(ctx->pelements[i]);
    }
    
    
    
    if(MSYM_SUCCESS != (ret = findCenterOfMass(ctx->elementsl,ctx->pelements,ctx->cm))) goto err;
    
    for(msym_element_t *a = ctx->elements; a < (ctx->elements+length); a++){
        vsub(a->v,ctx->cm,a->v);
    }
    
    double zero[3] = {0,0,0};
    
    if(MSYM_SUCCESS != (ret = findGeometry(length, ctx->pelements, zero, thresholds, &ctx->geometry, ctx->eigval, ctx->eigvec))) goto err;
    
    ctx->ext.elements = malloc(sizeof(msym_element_t[length]));
    memcpy(ctx->ext.elements, ctx->elements, sizeof(msym_element_t[length]));
    
    ctx->ext.set_elements_ptr = elements;
    
    return ret;
err:
    free(ctx->elements);
    free(ctx->pelements);
    free(ctx->ext.elements);
    ctx->ext.elements = NULL;
    ctx->elements = NULL;
    ctx->pelements = NULL;
    ctx->elementsl = 0;

    return ret;
}

msym_error_t msymGetElements(msym_context ctx, int *length, msym_element_t **elements){
    msym_error_t ret = MSYM_SUCCESS;
    msym_element_t *relements = NULL;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->elements == NULL || ctx->ext.elements == NULL) {ret = MSYM_INVALID_ELEMENTS;goto err;}
    
    *elements = ctx->ext.elements;
    *length = ctx->elementsl;
    return ret;
err:
    free(relements);
    *elements = NULL;
    *length = 0;
    return ret;
}

msym_error_t msymGetEquivalenceSets(msym_context ctx, int *length, const msym_equivalence_set_t **es){
    msym_error_t ret = MSYM_SUCCESS;
    
    if(ctx->ext.es == NULL){ret = MSYM_INVALID_EQUIVALENCE_SET;goto err;}
    
    *es = ctx->ext.es;
    *length = ctx->esl;
    
err:
    return ret;
}

msym_error_t msymGetBasisFunctions(msym_context ctx, int *length, msym_basis_function_t **basis){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->basis == NULL) {ret = MSYM_INVALID_ORBITALS;goto err;}
    
    *length = ctx->basisl;
    *basis = ctx->basis;
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
        if(bf->element >= ctx->ext.set_elements_ptr && bf->element < ctx->ext.set_elements_ptr + ctx->elementsl){
            bf->element = bf->element - ctx->ext.set_elements_ptr + ctx->ext.elements;
        }
        else if(!(bf->element >= ctx->ext.elements && ctx->basis[i].element < ctx->ext.elements + ctx->elementsl)){
            ret = MSYM_INVALID_ORBITALS;
            msymSetErrorDetails("Basis function element not set");
            goto err;
        }
        
        if(bf->type != MSYM_BASIS_TYPE_REAL_SPHERICAL_HARMONIC){
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

msym_error_t msymGetSubgroups(msym_context ctx, int *sgl, const msym_subgroup_t **sg){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->pg == NULL) {ret = MSYM_INVALID_POINT_GROUP;goto err;}
    if(ctx->pg->perm == NULL) {ret = MSYM_INVALID_PERMUTATION;goto err;}
    
    if(ctx->sg == NULL){
        int sgmax = numberOfSubgroups(ctx->pg);
        if(MSYM_SUCCESS != (ret = findPermutationSubgroups(ctx->pg->order, ctx->pg->perm, sgmax, ctx->pg->sops, &ctx->sgl, &ctx->sg))) goto err;

        for(int i = 0;i < ctx->sgl;i++){
            if(MSYM_SUCCESS != (ret = findSubgroup(&ctx->sg[i], ctx->thresholds))) goto err;
        }
    }

    *sgl = ctx->sgl;
    *sg = ctx->sg;
    return ret;
    
err:
    return ret;
    
    
    
}

msym_error_t msymGetSALCSubspaces(msym_context ctx, int *l, const msym_subspace_t **ss){
    msym_error_t ret = MSYM_SUCCESS;
    if(NULL == ctx) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(NULL == ctx->salc_ss){
        if(MSYM_SUCCESS != (ret = msymGenerateSALCSubspaces(ctx))) goto err;
        if(NULL == ctx->salc_ss){ret = MSYM_INVALID_ORBITALS;goto err;}
    }
    
    *ss = ctx->salc_ss;
    *l = ctx->salc_ssl;
    
    return ret;
err:
    return ret;
}

msym_error_t msymGetPointGroup(msym_context ctx, msym_point_group_t **pg){
    msym_error_t ret = MSYM_SUCCESS;
    if(NULL == ctx) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(NULL == ctx->pg){ret = MSYM_INVALID_POINT_GROUP;goto err;}
    *pg = ctx->pg;    
err:
    return ret;
    
}

msym_error_t msymGetCharacterTable(msym_context ctx, const msym_character_table_t **ct){
    msym_error_t ret = MSYM_SUCCESS;
    if(NULL == ctx) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(NULL == ctx->pg){ret = MSYM_INVALID_POINT_GROUP;goto err;}
    if(NULL == ctx->pg->ct){
        msym_point_group_t *pg = ctx->pg;
        if(MSYM_SUCCESS != (ret = generateCharacterTable(pg->type, pg->n, pg->order, pg->sops, &pg->ct))) goto err;
    }
    *ct = ctx->pg->ct;
    
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

msym_error_t msymSetCenterOfMass(msym_context ctx, double cm[3]){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    vcopy(cm,ctx->cm);
err:
    return ret;
}

msym_error_t msymGetGeometry(msym_context ctx, msym_geometry_t *geometry){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->elements == NULL) {ret = MSYM_INVALID_ELEMENTS;goto err;}
    if(ctx->geometry == MSYM_GEOMETRY_UNKNOWN) {ret = MSYM_INVALID_GEOMETRY;goto err;}
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
    for(int i = 0;i < ctx->elementsl;i++){
        double abs = vabs(ctx->elements[i].v);
        r = r > abs ? r : abs;
    }
    *radius = r;

err:
    return ret;
}

msym_error_t msymGetSymmetryOperations(msym_context ctx, int *sopsl, const msym_symmetry_operation_t **sops){
    msym_error_t ret = MSYM_SUCCESS;
    msym_symmetry_operation_t *rsops = NULL;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->pg == NULL || ctx->pg->sops == NULL) {ret = MSYM_INVALID_POINT_GROUP;goto err;}
    
    *sops = ctx->pg->sops;
    *sopsl = ctx->pg->order;
    return ret;
err:
    free(rsops);
    *sops = NULL;
    *sopsl = 0;
    return ret;
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

msym_error_t ctxGetThresholds(msym_context ctx, msym_thresholds_t **thresholds){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->thresholds == NULL) ret = MSYM_INVALID_THRESHOLD;
    msym_thresholds_t *t = ctx->thresholds;
    if(t->angle < 1.0 && !signbit(t->angle) &&
       t->equivalence < 1.0 && !signbit(t->equivalence) &&
       t->geometry < 1.0 && !signbit(t->geometry) &&
       !signbit(t->eigfact) &&
       !signbit(t->orthogonalization) &&
       !signbit(t->zero) &&
       !signbit(t->permutation)){
        *thresholds = t;
    } else {
        ret = MSYM_INVALID_THRESHOLD;
        goto err;
    }
err:
    return ret;
}

msym_error_t ctxGetElements(msym_context ctx, int *l, msym_element_t **elements){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->elements == NULL) {ret = MSYM_INVALID_ELEMENTS; goto err;}
    *elements = (msym_element_t *) ctx->elements;
    *l = ctx->elementsl;
err:
    return ret;
}

msym_error_t ctxGetExternalElements(msym_context ctx, int *l, msym_element_t **elements){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->ext.elements == NULL) {ret = MSYM_INVALID_ELEMENTS; goto err;}
    *elements = (msym_element_t *) ctx->ext.elements;
    *l = ctx->elementsl;
err:
    return ret;
}


msym_error_t ctxUpdateExternalElementCoordinates(msym_context ctx){
    msym_error_t ret = MSYM_SUCCESS;
    if(NULL == ctx) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(NULL == ctx->elements || NULL == ctx->ext.elements) {ret = MSYM_INVALID_ELEMENTS; goto err;}
    
    msym_element_t *internal = ctx->elements, *external = ctx->ext.elements;
    
    for(int i = 0; i < ctx->elementsl;i++){
        vadd(internal[i].v, ctx->cm, external[i].v);
    }
    
err:
    return ret;
}


msym_error_t ctxGetInternalElement(msym_context ctx, msym_element_t *ext, msym_element_t **element){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->ext.elements == NULL) {ret = MSYM_INVALID_ELEMENTS;goto err;}
    if(ext < ctx->ext.elements || ext >= ctx->ext.elements+ctx->elementsl){
        msymSetErrorDetails("Element pointer (%p) outside memory block (%p -> %p)", ext, ctx->ext.elements, ctx->ext.elements + ctx->elementsl);
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


msym_error_t ctxGetElementPtrs(msym_context ctx, int *l, msym_element_t ***pelements){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT;goto err;}
    if(ctx->pelements == NULL) {ret = MSYM_INVALID_ELEMENTS; goto err;}
    *pelements = (msym_element_t **) ctx->pelements;
    *l = ctx->elementsl;
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
    int el = 0;
    msym_element_t *elements = NULL;
    
    if(MSYM_SUCCESS != (ret = ctxDestroyEquivalcenceSets(ctx))) goto err;
    if(MSYM_SUCCESS != (ret = msymGetElements(ctx, &el, &elements))) goto err; //a bit lazy
    if(MSYM_SUCCESS != (ret = copyEquivalenceSets(esl,es,&ctx->ext.es))) goto err;
    
    for(int i = 0; i < esl;i++){
        for(int j = 0;j < es[i].length;j++){
            ctx->ext.es[i].elements[j] = ctx->ext.es[i].elements[j] - ctx->elements + elements;
        }
    }
    
    ctx->ext.eesmap = calloc(ctx->elementsl, sizeof(msym_equivalence_set_t *));
    for(int i = 0;i < esl;i++){
        for(int j = 0;j < ctx->ext.es[i].length;j++){
            ctx->ext.eesmap[ctx->ext.es[i].elements[j] - ctx->ext.elements] = &ctx->ext.es[i];
        }
    }
    for(int i = 0; i < ctx->elementsl;i++){
        if(!ctx->ext.eesmap[i]){
            ret = MSYM_INVALID_EQUIVALENCE_SET;
            msymSetErrorDetails("Element %d does not map to any equivalence set",i);
            goto err;
        }
    }
    ctx->es = es;
    ctx->esl = esl;
    return ret;
err:
    free(ctx->ext.es);
    free(ctx->ext.eesmap);
    ctx->ext.es = NULL;
    ctx->ext.eesmap = NULL;
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

msym_error_t ctxGetExternalEquivalenceSets(msym_context ctx, int *esl, msym_equivalence_set_t **es){
    msym_error_t ret = MSYM_SUCCESS;
    
    if(ctx->ext.es == NULL){ret = MSYM_INVALID_EQUIVALENCE_SET;goto err;}
    
    *es = ctx->ext.es;
    *esl = ctx->esl;
    
err:
    return ret;
}

msym_error_t ctxGetExternalElementEquivalenceSetMap(msym_context ctx, msym_equivalence_set_t ***eesmap){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    if(ctx->ext.eesmap == NULL) {ret = MSYM_INVALID_EQUIVALENCE_SET;goto err;}
    *eesmap = ctx->ext.eesmap;
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

msym_error_t ctxGetSALCSubspaces(msym_context ctx, int *ssl, msym_subspace_t **ss, int **span){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    if(ctx->salc_ss == NULL) {ret = MSYM_INVALID_SUBSPACE; goto err;}
    *ssl = ctx->salc_ssl;
    *ss = ctx->salc_ss;
    *span = ctx->salc_span;
err:
    return ret;
}

msym_error_t ctxSetSALCSubspaces(msym_context ctx, int ssl, msym_subspace_t *ss, int *span){
    msym_error_t ret = MSYM_SUCCESS;
    if(MSYM_SUCCESS != (ret = ctxDestroySALCSubspaces(ctx))) goto err;
    ctx->salc_ssl = ssl;
    ctx->salc_ss = ss;
    ctx->salc_span = span;
err:
    return ret;
}

msym_error_t ctxGetGeometry(msym_context ctx, msym_geometry_t *g, double eigval[3], double eigvec[3][3]){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    if(ctx->geometry == MSYM_GEOMETRY_UNKNOWN) {ret = MSYM_INVALID_GEOMETRY; goto err;}
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
    ctxDestroySALCSubspaces(ctx);
    ctxDestroyBasisFunctions(ctx);
    free(ctx->elements);
    free(ctx->pelements);
    free(ctx->ext.eesmap);
    free(ctx->ext.elements);
    
    ctx->ext.set_elements_ptr = NULL;
    ctx->elements = NULL;
    ctx->pelements = NULL;
    ctx->ext.eesmap = NULL;
    ctx->ext.elements = NULL;
    ctx->elementsl = 0;
    ctx->geometry = MSYM_GEOMETRY_UNKNOWN;
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
    free(ctx->ext.eesmap);
    free(ctx->es);
    free(ctx->ext.es);
    ctx->ext.eesmap = NULL;
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
   
    free(ctx->pg->perm);
    free(ctx->pg->ct);
    free(ctx->pg->sops);
    free(ctx->pg);
    
    ctx->pg = NULL;
err:
    return ret;
}

msym_error_t ctxDestroySubgroups(msym_context ctx){
    msym_error_t ret = MSYM_SUCCESS;
    free(ctx->sg);
    ctx->sg = NULL;
    ctx->sgl = 0;
err:
    return ret;
}

msym_error_t ctxDestroyBasisFunctions(msym_context ctx){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    ctxDestroySALCSubspaces(ctx);
    free(ctx->basis);
    ctx->basis = NULL;
    ctx->basisl = 0;
err:
    return ret;
}

msym_error_t ctxDestroySALCSubspaces(msym_context ctx){
    msym_error_t ret = MSYM_SUCCESS;
    if(ctx == NULL) {ret = MSYM_INVALID_CONTEXT; goto err;}
    freeSALCSubspaces(ctx->salc_ssl, ctx->salc_ss);
    free(ctx->salc_span);
    ctx->salc_ss = NULL;
    ctx->salc_span = NULL;
    ctx->salc_ssl = 0;
err:
    return ret;
}
