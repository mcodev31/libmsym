//
//  msym.c
//  libmsym
//
//  Created by Marcus Johansson on 30/01/15.
//  Copyright (c) 2015 Marcus Johansson. 
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "msym.h"
#include "context.h"
#include "symmetry.h"
#include "equivalence_set.h"
#include "point_group.h"
#include "symmetrize.h"
#include "linalg.h"
#include "subspace.h"

#include "debug.h"

msym_error_t msymFindSymmetry(msym_context ctx){
    msym_error_t ret = MSYM_SUCCESS;
    int elementsl = 0, esl = 0;
    msym_element_t *elements = NULL;
    msym_thresholds_t *t = NULL;
    msym_equivalence_set_t *es = NULL, *des = NULL;;
    msym_point_group_t *pg = NULL;
    int sopsl = 0;
    msym_symmetry_operation_t *sops = NULL;
    msym_equivalence_set_t *ses = NULL;
    int sesl = 0;
    msym_point_group_t *fpg = NULL;
    
    if(MSYM_SUCCESS != (ret = ctxGetElements(ctx, &elementsl, &elements))) goto err;
    
    if(MSYM_SUCCESS != (ret = ctxGetThresholds(ctx, &t))) goto err;
    
    if(MSYM_SUCCESS != (ret = ctxGetEquivalenceSets(ctx, &esl, &des))){
        if(MSYM_SUCCESS != (ret = msymFindEquivalenceSets(ctx))) goto err;
    }
    
    if(MSYM_SUCCESS != (ret = ctxGetEquivalenceSets(ctx, &esl, &es))) goto err;
    if(MSYM_SUCCESS != (ret = ctxGetPointGroup(ctx, &pg))){
        if(MSYM_SUCCESS != (ret = findSymmetryOperations(esl,es,t,&sopsl,&sops))) goto err;
        if(MSYM_SUCCESS != (ret = findPointGroup(sopsl, sops, t, &fpg))) goto err;
        pg = fpg;
        if(MSYM_SUCCESS != (ret = ctxSetPointGroup(ctx, pg))) {
            free(pg);
            goto err;
        }
    }
    
    if(NULL != fpg || isLinearPointGroup(pg)){
        // Rebuild equivalence sets after determining point group in case they are very similar
        if(MSYM_SUCCESS != (ret = ctxReduceLinearPointGroup(ctx))) goto err;
        
        if(MSYM_SUCCESS != (ret = splitPointGroupEquivalenceSets(pg, esl, es, &sesl, &ses, t))) goto err;
        if(MSYM_SUCCESS != (ret = ctxSetEquivalenceSets(ctx, sesl, ses))) goto err;
        ses = NULL; sesl = 0;
        if(MSYM_SUCCESS != (ret = ctxGetEquivalenceSets(ctx, &esl, &es))) goto err;
    }
    
    if(MSYM_SUCCESS != (ret = msymFindEquivalenceSetPermutations(ctx))) goto err;
    
    if(MSYM_SUCCESS != (ret = ctxGetEquivalenceSets(ctx, &esl, &es))) goto err; //This is only for printing, since permutation may regenerate sets
    
    free(sops);
    return ret;
    
err:
    free(ses);
    free(sops);
    if(des == NULL) {
        ctxDestroyEquivalcenceSets(ctx);
    }
    return ret;
}

msym_error_t msymSetPointGroupByName(msym_context ctx, const char *name){
    msym_error_t ret = MSYM_SUCCESS;
    msym_point_group_t *pg = NULL, *ppg = NULL;
    msym_thresholds_t *t = NULL;
    
    
    if(MSYM_SUCCESS != (ret = ctxGetThresholds(ctx, &t))) goto err;
    
    if(MSYM_SUCCESS != (ret = ctxGetPointGroup(ctx, &ppg))){
        double transform[3][3];
        mleye(3,transform);
        if(MSYM_SUCCESS != (ret = generatePointGroupFromName(name, transform, t, &pg))) goto err;
    } else if(MSYM_SUCCESS != (ret = generatePointGroupFromName(name, ppg->transform, t, &pg))) goto err;
    
    if(MSYM_SUCCESS != (ret = ctxSetPointGroup(ctx, pg))) goto err;
    
    return ret;
    
err:
    free(pg);
    return ret;
}

msym_error_t msymSetPointGroupByType(msym_context ctx, msym_point_group_type_t type, int n){
    msym_error_t ret = MSYM_SUCCESS;
    msym_point_group_t *pg = NULL, *ppg = NULL;
    msym_thresholds_t *t = NULL;
    
    
    if(MSYM_SUCCESS != (ret = ctxGetThresholds(ctx, &t))) goto err;
    
    if(MSYM_SUCCESS != (ret = ctxGetPointGroup(ctx, &ppg))){
        double transform[3][3];
        mleye(3,transform);
        if(MSYM_SUCCESS != (ret = generatePointGroupFromType(type, n, transform, t, &pg))) goto err;
    } else if(MSYM_SUCCESS != (ret = generatePointGroupFromType(type, n, ppg->transform, t, &pg))) goto err;
    
    if(MSYM_SUCCESS != (ret = ctxSetPointGroup(ctx, pg))) goto err;
    
    return ret;
    
err:
    free(pg);
    return ret;
}

msym_error_t msymGenerateElements(msym_context ctx, int length, msym_element_t elements[length]){
    msym_error_t ret = MSYM_SUCCESS;
    msym_point_group_t *pg = NULL;
    msym_thresholds_t *t = NULL;
    msym_element_t *gelements = NULL;
    msym_equivalence_set_t *es = NULL;
    msym_element_t **pelements = NULL;
    double err = 0.0;
    double cm[3];
    
    int glength = 0, plength = 0, esl = 0;
    if(MSYM_SUCCESS != (ret = ctxGetThresholds(ctx, &t))) goto err;
    if(MSYM_SUCCESS != (ret = msymGetCenterOfMass(ctx, cm))) goto err;

    if(MSYM_SUCCESS != (ret = ctxGetPointGroup(ctx, &pg))) goto err;
    if(MSYM_SUCCESS != (ret = generateEquivalenceSet(pg, length, elements, cm, &glength, &gelements, &esl, &es,t))) goto err;
    if(MSYM_SUCCESS != (ret = ctxSetElements(ctx, glength, gelements))) goto err;
    if(MSYM_SUCCESS != (ret = ctxGetElementPtrs(ctx, &plength, &pelements))) goto err;
    if(plength != glength){
        ret = MSYM_INVALID_ELEMENTS;
        msymSetErrorDetails("Inconsistency detected when setting elements");
        goto err;
    }
    for(int i = 0;i < esl;i++){
        for(int j = 0;j < es[i].length;j++){
            long int index = es[i].elements[j] - gelements;
            es[i].elements[j] = pelements[index];
        }
    }
    if(MSYM_SUCCESS != (ret = ctxSetEquivalenceSets(ctx, esl, es))) goto err;
    es = NULL; esl = 0;
    if(MSYM_SUCCESS != (ret = msymFindEquivalenceSetPermutations(ctx))) goto err;
    if(MSYM_SUCCESS != (ret = msymSymmetrizeElements(ctx, &err))) goto err;
    if(MSYM_SUCCESS != (ret = msymSetCenterOfMass(ctx, cm))) goto err;
    free(gelements);
    return ret;
    
err:
    free(gelements);
    free(es);
    return ret;
}

msym_error_t msymFindEquivalenceSets(msym_context ctx){
    msym_error_t ret = MSYM_SUCCESS;
    int pelementsl = 0;
    msym_element_t **pelements = NULL;
    msym_thresholds_t *t = NULL;
    msym_point_group_t *pg = NULL;
    msym_geometry_t g = MSYM_GEOMETRY_UNKNOWN;
    double eigvec[3][3];
    double eigval[3];
    int esl = 0;
    msym_equivalence_set_t *es;
    
    if(MSYM_SUCCESS != (ret = ctxGetElementPtrs(ctx, &pelementsl, &pelements))) goto err;
    if(MSYM_SUCCESS != (ret = ctxGetThresholds(ctx, &t))) goto err;
    if(MSYM_SUCCESS != (ret = ctxGetPointGroup(ctx, &pg))) {
        if(MSYM_SUCCESS != (ret = ctxGetGeometry(ctx, &g, eigval, eigvec))) goto err;
        if(MSYM_SUCCESS != (ret = findEquivalenceSets(pelementsl, pelements, g, &esl, &es, t))) goto err;
    } else {
        if(MSYM_SUCCESS != (ret = findPointGroupEquivalenceSets(pg, pelementsl, pelements, &esl, &es, t))) goto err;
    }
    if(MSYM_SUCCESS != (ret = ctxSetEquivalenceSets(ctx, esl, es))) goto err;
err:
    return ret;
}


msym_error_t msymAlignAxes(msym_context ctx){
    
    msym_error_t ret = MSYM_SUCCESS;
    msym_element_t *elements = NULL;
    msym_point_group_t *pg;
    int elementsl = 0;
    double zero[3] = {0,0,0};
    
    if(MSYM_SUCCESS != (ret = ctxGetElements(ctx, &elementsl, &elements))) goto err;
    if(MSYM_SUCCESS != (ret = ctxGetPointGroup(ctx, &pg))) goto err;
    
    
    if(pg->sops == NULL || pg->order == 0){
        msymSetErrorDetails("No symmetry operations in point group");
        ret = MSYM_INVALID_POINT_GROUP;
        goto err;
    }
    
    if(MSYM_SUCCESS != (ret = msymSetCenterOfMass(ctx, zero))) goto err;
    
    for(int i = 0; i < elementsl;i++) mvmul(elements[i].v, pg->transform, elements[i].v);
    for(int i = 0; i < pg->order;i++) mvmul(pg->sops[i].v, pg->transform, pg->sops[i].v);
    mleye(3,pg->transform);
    if(MSYM_SUCCESS != (ret = ctxUpdateExternalElementCoordinates(ctx))) goto err;
    
err:
    return ret;
}

msym_error_t msymGetAlignmentAxes(msym_context ctx, double primary[3], double secondary[3]){
    msym_error_t ret = MSYM_SUCCESS;
    msym_point_group_t *pg;
    
    if(MSYM_SUCCESS != (ret = ctxGetPointGroup(ctx, &pg))) goto err;
    
    double m[3][3], x[3] = {1,0,0}, z[3] = {0,0,1};
    minv(pg->transform,m);
    mvmul(z, m, primary);
    mvmul(x, m, secondary);
    
err:
    return ret;

}


msym_error_t msymGetAlignmentTransform(msym_context ctx, double transform[3][3]){
    msym_error_t ret = MSYM_SUCCESS;
    msym_point_group_t *pg;
    
    if(MSYM_SUCCESS != (ret = ctxGetPointGroup(ctx, &pg))) goto err;
    
    mcopy(pg->transform, transform);
    
err:
    return ret;
    
}

msym_error_t msymSetAlignmentTransform(msym_context ctx, double transform[3][3]){
    msym_error_t ret = MSYM_SUCCESS;
    msym_point_group_t *pg;
    msym_element_t *elements = NULL;
    msym_thresholds_t *t = NULL;
    msym_equivalence_set_t *es = NULL;
    int elementsl = 0, esl = 0;
    double m[3][3];
    
    if(MSYM_SUCCESS != (ret = ctxGetThresholds(ctx, &t))) goto err;
    if(MSYM_SUCCESS != (ret = ctxGetElements(ctx, &elementsl, &elements))){
        elements = NULL;
        elementsl = 0;
    }
    
    if(MSYM_SUCCESS != (ret = ctxGetEquivalenceSets(ctx, &esl, &es))){
        es = NULL;
        esl = 0;
    }
    
    if(MSYM_SUCCESS != (ret = ctxGetPointGroup(ctx, &pg))) goto err;
    
    if(pg->sops == NULL || pg->order == 0){
        msymSetErrorDetails("No symmetry operations in point group for setting alignment transform");
        ret = MSYM_INVALID_POINT_GROUP;
        goto err;
    }
    /* Don't transform elements if we don't have an equivalence set the current pg is set manually */
    if(NULL != es){
        for(int i = 0; i < elementsl;i++) mvmul(elements[i].v, pg->transform, elements[i].v);
    }
    for(int i = 0; i < pg->order;i++) mvmul(pg->sops[i].v, pg->transform, pg->sops[i].v);
    
    minv(transform,m);
    mcopy(transform, pg->transform);
    if(NULL != es){
        for(int i = 0; i < elementsl;i++) mvmul(elements[i].v, m, elements[i].v);
    }
    for(int i = 0; i < pg->order;i++) mvmul(pg->sops[i].v, m, pg->sops[i].v);
    
err:
    return ret;
}

msym_error_t msymSetAlignmentAxes(msym_context ctx, double primary[3], double secondary[3]){
    
    msym_error_t ret = MSYM_SUCCESS;
    msym_point_group_t *pg;
    msym_element_t *elements = NULL;
    msym_thresholds_t *t = NULL;
    msym_equivalence_set_t *es = NULL;
    int elementsl = 0, esl = 0;
    double x[3] = {1,0,0}, z[3] = {0,0,1}, m[3][3], p[3], s[3];
    
    vnorm2(primary, p);
    vnorm2(secondary,s);
    
    if(MSYM_SUCCESS != (ret = ctxGetThresholds(ctx, &t))) goto err;
    if(MSYM_SUCCESS != (ret = ctxGetElements(ctx, &elementsl, &elements))){
        elements = NULL;
        elementsl = 0;
    }
    
    if(MSYM_SUCCESS != (ret = ctxGetEquivalenceSets(ctx, &esl, &es))){
        es = NULL;
        esl = 0;
    }
    
    if(MSYM_SUCCESS != (ret = ctxGetPointGroup(ctx, &pg))) goto err;
    
    if(pg->sops == NULL || pg->order == 0){
        msymSetErrorDetails("No symmetry operations in point group for setting alignment axes");
        ret = MSYM_INVALID_POINT_GROUP;
        goto err;
    }
    
    if(!vperpendicular(primary, secondary, t->angle)) {
        msymSetErrorDetails("Alignment axes are not orthogonal");
        ret = MSYM_INVALID_AXES;
        goto err;
    }
    
    /* Don't transform elements if we don't have an equivalence set the current pg is set manually */
    if(NULL != es){
        for(int i = 0; i < elementsl;i++) mvmul(elements[i].v, pg->transform, elements[i].v);
    }
    for(int i = 0; i < pg->order;i++) mvmul(pg->sops[i].v, pg->transform, pg->sops[i].v);
    
    vproj_plane(s, p, s);
    malign(p,z,pg->transform);
    mvmul(s, pg->transform, s);
    malign(s,x,m);
    mmmul(m,pg->transform,pg->transform);
    minv(pg->transform,m);
    
    if(NULL != es){
        for(int i = 0; i < elementsl;i++) mvmul(elements[i].v, m, elements[i].v);
    }
    for(int i = 0; i < pg->order;i++) mvmul(pg->sops[i].v, m, pg->sops[i].v);
    
    
err:
    return ret;
}

msym_error_t msymSelectSubgroup(msym_context ctx, const msym_subgroup_t *sg){
    msym_error_t ret = MSYM_SUCCESS;
    msym_subgroup_t *sgs;
    msym_point_group_t *pg;
    msym_thresholds_t *t = NULL;
    int sgl = 0;
    
    if(MSYM_SUCCESS != (ret = ctxGetSubgroups(ctx, &sgl, &sgs))) goto err;
    if(sg < sgs || sg >= sgs + sgl){
        msymSetErrorDetails("Subgroup not available in current context");
        ret = MSYM_INVALID_SUBGROUPS;
        goto err;
    }
    if(MSYM_SUCCESS != (ret = ctxGetThresholds(ctx, &t))) goto err;
    if(MSYM_SUCCESS != (ret = pointGroupFromSubgroup(sg, t, &pg))) goto err;
    if(MSYM_SUCCESS != (ret = ctxSetPointGroup(ctx, pg))) goto err;
    if(MSYM_SUCCESS != (ret = msymFindEquivalenceSets(ctx))) goto err;
    if(MSYM_SUCCESS != (ret = msymFindEquivalenceSetPermutations(ctx))) goto err;

err:
    return ret;
}

msym_error_t msymSymmetrizeElements(msym_context ctx, double *oerr){
    msym_error_t ret = MSYM_SUCCESS;
    
    msym_point_group_t *pg = NULL;
    msym_equivalence_set_t *es = NULL;
    msym_element_t *elements = NULL;
    
    msym_permutation_t **perm = NULL;
    msym_thresholds_t *t = NULL;
    double error = 0.0;
    int perml = 0, esl = 0, elementsl = 0, sopsl = 0;
    
    if(MSYM_SUCCESS != (ret = ctxGetThresholds(ctx, &t))) goto err;
    if(MSYM_SUCCESS != (ret = ctxGetElements(ctx, &elementsl, &elements))) goto err;
    if(MSYM_SUCCESS != (ret = ctxGetPointGroup(ctx, &pg))) goto err;
    if(MSYM_SUCCESS != (ret = ctxGetEquivalenceSets(ctx, &esl, &es))){
        if(MSYM_SUCCESS != (ret = msymFindEquivalenceSets(ctx))) goto err;
        if(MSYM_SUCCESS != (ret = msymFindEquivalenceSetPermutations(ctx))) goto err;
        if(MSYM_SUCCESS != (ret = ctxGetEquivalenceSets(ctx, &esl, &es))) goto err;
    }
    if(MSYM_SUCCESS != (ret = ctxGetEquivalenceSetPermutations(ctx, &perml, &sopsl, &perm))) goto err;
    if(sopsl != pg->order || perml != esl) {
        msymSetErrorDetails("Detected inconsistency between point group, equivalence sets and permutations");
        ret = MSYM_INVALID_PERMUTATION;
        goto err;
    }
    
    if(MSYM_SUCCESS != (ret = symmetrizeElements(pg, esl, es, perm, t, &error))) goto err;
    
    if(MSYM_SUCCESS != (ret = ctxUpdateGeometry(ctx))) goto err;
    
    if(MSYM_SUCCESS != (ret = ctxUpdateExternalElementCoordinates(ctx))) goto err;
    
    *oerr = error;
err:
    return ret;
}

msym_error_t msymApplyTranslation(msym_context ctx, msym_element_t *ext, double v[3]){
    msym_error_t ret = MSYM_SUCCESS;
    
    msym_point_group_t *pg;
    msym_equivalence_set_t *es, *ees;
    msym_element_t *eelements;
    msym_equivalence_set_t **eesmap = NULL;
    msym_permutation_t **perm;
    msym_thresholds_t *t = NULL;
    int perml = 0, esl = 0, eesl = 0, eelementsl = 0, sopsl = 0;
    
    if(MSYM_SUCCESS != (ret = ctxGetThresholds(ctx, &t))) goto err;
    if(MSYM_SUCCESS != (ret = ctxGetPointGroup(ctx, &pg))) goto err;
    if(MSYM_SUCCESS != (ret = ctxGetExternalElements(ctx, &eelementsl, &eelements))) goto err;
    if(MSYM_SUCCESS != (ret = ctxGetEquivalenceSets(ctx, &esl, &es))){
        if(MSYM_SUCCESS != (ret = msymFindEquivalenceSets(ctx))) goto err;
        if(MSYM_SUCCESS != (ret = msymFindEquivalenceSetPermutations(ctx))) goto err;
        if(MSYM_SUCCESS != (ret = ctxGetEquivalenceSets(ctx, &esl, &es))) goto err;
    }
    if(MSYM_SUCCESS != (ret = ctxGetExternalEquivalenceSets(ctx, &eesl, &ees))) goto err;
    if(MSYM_SUCCESS != (ret = ctxGetExternalElementEquivalenceSetMap(ctx, &eesmap))) goto err;
    
    if(MSYM_SUCCESS != (ret = ctxGetEquivalenceSetPermutations(ctx, &perml, &sopsl, &perm))) goto err;
    if(sopsl != pg->order || perml != esl) {
        msymSetErrorDetails("Detected inconsistency between point group, equivalence sets and permutations");
        ret = MSYM_INVALID_PERMUTATION;
        goto err;
    }
    
    int esmi = (int)(ext - eelements);
    
    if(esmi > eelementsl) {
        msymSetErrorDetails("Element outside of memory block of external elements");
        ret = MSYM_INVALID_ELEMENTS;
        goto err;
    }
    
    int fesi = (int)(eesmap[esmi] - ees);
    msym_equivalence_set_t *fes = eesmap[esmi];
    int fi = 0;
    for(fi = 0;fi < fes->length;fi++){
        if(fes->elements[fi] == ext) break;
    }
    
    if(fi >= fes->length){
        msymSetErrorDetails("Could not find index of element %s in equivalence set %d", ext->name, fesi);
        ret = MSYM_INVALID_ELEMENTS;
        goto err;
    }
    
    if(MSYM_SUCCESS != (ret = symmetrizeTranslation(pg, &es[fesi], perm[fesi], fi, v))) goto err;
    
    if(MSYM_SUCCESS != (ret = ctxUpdateExternalElementCoordinates(ctx))) goto err;
    
    return ret;
err:
    return ret;
}

msym_error_t msymGenerateSubrepresentationSpaces(msym_context ctx){
    msym_error_t ret = MSYM_SUCCESS;
    
    msym_point_group_t *pg = NULL;
    msym_basis_function_t *basis = NULL;
    msym_equivalence_set_t *es = NULL;
    msym_equivalence_set_t **eesmap = NULL;
    msym_permutation_t **perm = NULL;
    msym_thresholds_t *t = NULL;
    msym_subrepresentation_space_t *srs = NULL;
    msym_basis_function_t **srsbf = NULL;
    msym_element_t *elements = NULL;
    const msym_subgroup_t *sg = NULL;
    int *span = NULL;
    
    int basisl = 0, esl = 0, perml = 0, sopsl = 0, srsl = 0, elementsl = 0, sgl = 0;
    
    if(MSYM_SUCCESS != (ret = ctxGetThresholds(ctx, &t))) goto err;
    if(MSYM_SUCCESS != (ret = ctxGetExternalElements(ctx, &elementsl, &elements))) goto err;
    if(MSYM_SUCCESS != (ret = ctxGetPointGroup(ctx, &pg))) goto err;
    if(pg->ct == NULL){
        if(MSYM_SUCCESS != (ret = generateCharacterTable(pg->type, pg->n, pg->order, pg->sops, &pg->ct))) goto err;
    }
    if(MSYM_SUCCESS != (ret = ctxGetExternalEquivalenceSets(ctx, &esl, &es))) goto err;
    if(MSYM_SUCCESS != (ret = ctxGetExternalElementEquivalenceSetMap(ctx, &eesmap))) goto err;
    if(MSYM_SUCCESS != (ret = ctxGetBasisFunctions(ctx, &basisl, &basis))) goto err;
    if(MSYM_SUCCESS != (ret = ctxGetEquivalenceSetPermutations(ctx, &perml, &sopsl, &perm))) goto err;
    if(sopsl != pg->order || perml != esl) {ret = MSYM_INVALID_PERMUTATION; goto err;}
    
    if(MSYM_SUCCESS != (ret = msymGetSubgroups(ctx, &sgl, &sg))) goto err;
    
    if(MSYM_SUCCESS != (ret = generateSubrepresentationSpaces(pg, sgl, sg, esl, es, perm, basisl, basis, elements, eesmap, t, &srsl, &srs, &srsbf, &span))) goto err;
    
    if(MSYM_SUCCESS != (ret = ctxSetSubrepresentationSpaces(ctx,srsl,srs,srsbf,span))) goto err;
    
    return ret;
err:
    freeSubrepresentationSpaces(srsl, srs);
    free(srs);
    free(span);
    return ret;
}

msym_error_t msymGetSALCs(msym_context ctx, int l, double c[l][l], int species[l], msym_partner_function_t pf[l]){
    msym_error_t ret = MSYM_SUCCESS;
    
    
    msym_subrepresentation_space_t *srs = NULL;
    msym_basis_function_t *basis = NULL;
    
    int *span = NULL;
    
    int srsl = 0, basisl = 0;
    
    if(MSYM_SUCCESS != (ret = ctxGetBasisFunctions(ctx, &basisl, &basis))) goto err;
    
    if(MSYM_SUCCESS != (ret = ctxGetSubrepresentationSpaces(ctx, &srsl, &srs, &span))){
        if(MSYM_SUCCESS != (ret = msymGenerateSubrepresentationSpaces(ctx))) goto err;
        if(MSYM_SUCCESS != (ret = ctxGetSubrepresentationSpaces(ctx, &srsl, &srs, &span))) goto err;
    }
    
    if(l != basisl){
        ret = MSYM_INVALID_INPUT;
        msymSetErrorDetails("Supplied SALC matrix size (%dx%d) does not match number of basis functions (%d)",l,l,basisl);
        goto err;
    }
    
    memset(c,0,sizeof(double[l][l]));
    int wf = 0;
    for(int i = 0;i < srsl;i++){
        int s = srs[i].s;
        for(int j = 0;j < srs[i].salcl;j++){
            int pwf = wf;
            double (*mpf)[srs[i].salc[j].fl] = srs[i].salc[j].pf;
            for(int d = 0;d < srs[i].salc[j].d;d++){
                if(wf >= basisl){
                    ret = MSYM_INVALID_SUBSPACE;
                    msymSetErrorDetails("Generated more SALCs than the number of basis functions (%d)", basisl);
                    goto err;
                }
                for(int f = 0;f < srs[i].salc[j].fl;f++){
                    int index = (int)(srs[i].salc[j].f[f] - basis);
                    c[wf][index] = mpf[d][f];
                }
                if(NULL != pf){
                    pf[wf].i = pwf;
                    pf[wf].d = d;
                }
                if(NULL != species) species[wf] = s;
                wf++;
            }
        }
    }
    
    if(wf != basisl){
        ret = MSYM_INVALID_BASIS_FUNCTIONS;
        msymSetErrorDetails("Number of generated SALC wavefunctions (%d) does not match orbital basis (%d)",wf,basisl);
        goto err;
    }
    
err:
    return ret;
    
}

msym_error_t msymSymmetrySpeciesComponents(msym_context ctx, int wfl, double *wf, int sl, double *s){
    msym_error_t ret = MSYM_SUCCESS;
    
    msym_point_group_t *pg = NULL;
    msym_subrepresentation_space_t *srs = NULL;
    msym_basis_function_t *basis = NULL;
    int *span = NULL;
    
    int srsl = 0, basisl = 0;
    
    if(MSYM_SUCCESS != (ret = ctxGetPointGroup(ctx, &pg))) goto err;
    if(pg->ct == NULL){
        if(MSYM_SUCCESS != (ret = generateCharacterTable(pg->type, pg->n, pg->order, pg->sops, &pg->ct))) goto err;
    }
    
    if(MSYM_SUCCESS != (ret = ctxGetBasisFunctions(ctx, &basisl, &basis))) goto err;
    
    if(basisl != wfl) {
        ret = MSYM_INVALID_INPUT;
        msymSetErrorDetails("Supplied coefficient vector size (%d) does not match number of basis functions (%d)",wfl,basisl);
        goto err;
    }
    
    if(sl != pg->ct->d) {
        ret = MSYM_INVALID_INPUT;
        msymSetErrorDetails("Supplied symmetry species vector size (%d) does not match character table (%d)",sl,pg->ct->d);
        goto err;
    }
    
    if(MSYM_SUCCESS != (ret = ctxGetSubrepresentationSpaces(ctx, &srsl, &srs, &span))){
        if(MSYM_SUCCESS != (ret = msymGenerateSubrepresentationSpaces(ctx))) goto err;
        if(MSYM_SUCCESS != (ret = ctxGetSubrepresentationSpaces(ctx, &srsl, &srs, &span))) goto err;
    }
    
    if(MSYM_SUCCESS != (ret = symmetrySpeciesComponents(pg, srsl, srs, basisl, basis, wf, s))) goto err;
    
err:
    return ret;
    
}

msym_error_t msymSymmetrizeWavefunctions(msym_context ctx, int l, double c[l][l], int species[l], msym_partner_function_t pf[l]){
    msym_error_t ret = MSYM_SUCCESS;
    msym_point_group_t *pg = NULL;
    msym_subrepresentation_space_t *srs = NULL;
    msym_basis_function_t *basis = NULL;
    int *span = NULL;
    
    int srsl = 0, basisl = 0;
    
    if(MSYM_SUCCESS != (ret = ctxGetPointGroup(ctx, &pg))) goto err;
    if(pg->ct == NULL){
        if(MSYM_SUCCESS != (ret = generateCharacterTable(pg->type, pg->n, pg->order, pg->sops, &pg->ct))) goto err;
    }
    
    
    if(MSYM_SUCCESS != (ret = ctxGetBasisFunctions(ctx, &basisl, &basis))) goto err;
    
    if(basisl != l) {
        ret = MSYM_INVALID_INPUT;
        msymSetErrorDetails("Supplied wavefunction matrix size (%d) does not match number of basis functions (%d)",l,basisl);
        goto err;
    }
    
    if(MSYM_SUCCESS != (ret = ctxGetSubrepresentationSpaces(ctx, &srsl, &srs, &span))){
        if(MSYM_SUCCESS != (ret = msymGenerateSubrepresentationSpaces(ctx))) goto err;
        if(MSYM_SUCCESS != (ret = ctxGetSubrepresentationSpaces(ctx, &srsl, &srs, &span))) goto err;
    }
    
    if(MSYM_SUCCESS != (ret = symmetrizeWavefunctions(pg, srsl, srs, span, basisl, basis, c , c, species, pf))) goto err;
    
err:
    return ret;
}

msym_error_t msymFindEquivalenceSetPermutations(msym_context ctx) {
    msym_error_t ret = MSYM_SUCCESS;
    //We can't allocate this as a double[][] unless we typecast it every time, since the compiler doesn't have the indexing information in the context
    msym_permutation_t **perm = NULL;
    msym_permutation_t *bperm = NULL;
    msym_point_group_t *pg = NULL;
    msym_equivalence_set_t *es = NULL;
    msym_thresholds_t *t = NULL;
    double (**esv)[3] = NULL;
    int esl = 0;
    
    if(MSYM_SUCCESS != (ret = ctxGetThresholds(ctx, &t))) goto err;
    if(MSYM_SUCCESS != (ret = ctxGetPointGroup(ctx, &pg))) goto err;
    if(MSYM_SUCCESS != (ret = ctxGetEquivalenceSets(ctx, &esl, &es))) goto err;
    
    perm = (msym_permutation_t**)malloc(esl*sizeof(msym_permutation_t*) + esl*pg->order*sizeof(msym_permutation_t));
    bperm = (msym_permutation_t*)(perm + esl);
    memset(bperm,0,esl*pg->order*sizeof(msym_permutation_t));
    
    
    
    for(int i = 0; i < esl;i++){
        perm[i] = bperm + i*pg->order;
        if(es[i].length > pg->order){
            msymSetErrorDetails("Equivalence set has more elements (%d) than the order of the point group %s (%d)",es[i].length,pg->name,pg->order);
            ret = MSYM_INVALID_EQUIVALENCE_SET;
            goto err;
        }
    }
    /*
    if(perm == NULL){
        perm = (msym_permutation_t**)malloc(esl*sizeof(msym_permutation_t*) + esl*pg->sopsl*sizeof(msym_permutation_t));
        bperm = (msym_permutation_t*)(perm + esl);
        memset(bperm,0,esl*pg->sopsl*sizeof(msym_permutation_t));
        for(int i = 0; i < esl;i++){ //This really shouldn't happen
            perm[i] = bperm + i*pg->sopsl;
            if(es[i].length > pg->order){
                msymSetErrorDetails("Equivalence set has more elements (%d) than the order of the point group %s (%d)",es[i].length,pg->name,pg->order);
                ret = MSYM_INVALID_EQUIVALENCE_SET;
                goto err;
            }
        }
    }*/
    
    esv = malloc(sizeof(double (*[pg->order])[3]));
    for(int i = 0; i < esl;i++){
        for(int j = 0; j < es[i].length;j++){
            esv[j] = &es[i].elements[j]->v;
        }
        
        for(int j = 0; j < pg->order;j++){
            if(MSYM_SUCCESS != (ret = findPermutation(&pg->sops[j], es[i].length, esv, t, &perm[i][j]))) goto err;
        }
    }
        
    if(MSYM_SUCCESS != (ret = ctxSetEquivalenceSetPermutations(ctx, esl, pg->order, perm))) goto err;
    
    free(esv);
    return ret;
    
err:
    free(esv);
    free(perm);
    return ret;
}

