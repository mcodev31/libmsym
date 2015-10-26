//
//  example.c
//  libmsym
//
//  Created by Marcus Johansson on 24/04/15.
//  Copyright (c) 2015 Marcus Johansson. 
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "msym.h"
#include "example.h"

int read_xyz(const char *name, msym_element_t **ratoms);

int example(const char* in_file, msym_thresholds_t *thresholds){
    msym_error_t ret = MSYM_SUCCESS;
    msym_element_t *elements = NULL;
    
    const char *error = NULL;
    char point_group[6];
    double cm[3], radius = 0.0, symerr = 0.0;
    
    /* Do not free these variables */
    msym_element_t *melements = NULL;
    msym_basis_function_t *mbfs = NULL;
    /* these are not mutable */
    const msym_symmetry_operation_t *msops = NULL;
    const msym_subgroup_t *msg = NULL;
    const msym_subrepresentation_space_t *msrs = NULL;
    const msym_character_table_t *mct = NULL;
    double *irrep = NULL;
    
    
    
    
    msym_basis_function_t *bfs = NULL;
    
    int msgl = 0, msopsl = 0, mlength = 0, msrsl = 0, mbfsl = 0;
    int orbitalsl = 0, bfsl = 0;
    
    //char *orbitals[15] = {"2px", "2py", "2pz", "3d2-", "3d1-", "3d0", "3d1+", "3d2+","4f3-", "4f2-", "4f1-", "4f0", "4f1+", "4f2+", "4f3+"};
    
    char *orbitals[13] = {"7i6-", "7i5-", "7i4-", "7i3-", "7i2-", "7i1-", "7i0", "7i1+", "7i2+", "7i3+", "7i4+", "7i5+", "7i6+"};
    
    /* This function reads xyz files.
     * It initializes an array of msym_element_t to 0,
     * then sets the coordinates and name of the elements */
    int length = read_xyz(in_file, &elements);
    if(length <= 0) return -1;
    
    orbitalsl = sizeof(orbitals)/sizeof(char*);
    bfsl = orbitalsl*length;
    bfs = calloc(bfsl, sizeof(msym_basis_function_t));
    double (*salcs)[bfsl] = calloc(bfsl, sizeof(*salcs)); // SALCs in matrix form, and input for symmetrization
    double *cmem = calloc(bfsl, sizeof(double)); // Some temporary memory
    int *species = calloc(bfsl, sizeof(*species));
    msym_partner_function_t *pf = calloc(bfsl, sizeof(*pf));
    
    /* Create a context */
    msym_context ctx = msymCreateContext();
    
    if(NULL != thresholds){
        if(MSYM_SUCCESS != (ret = msymSetThresholds(ctx, thresholds))) goto err;
    }
    
    /* Use default thresholds otherwise call:
     * msymSetThresholds(msym_context ctx, msym_thresholds_t *thresholds); */
    
    /* Set elements */
    if(MSYM_SUCCESS != (ret = msymSetElements(ctx, length, elements))) goto err;
    
    /* Get elements msym elements */
    if(MSYM_SUCCESS != (ret = msymGetElements(ctx, &mlength, &melements))) goto err;
    
    for(int i = 0, k = 0;i < length;i++){
        for(int j = 0; j < orbitalsl; j++){
            snprintf(bfs[k].name,sizeof(bfs[i].name),"%s",orbitals[j]);
            bfs[k].element = &melements[i];
            bfs[k].type = MSYM_BASIS_TYPE_REAL_SPHERICAL_HARMONIC;
            k++;
        }
    }
    
    printf("Created %d basis functions\n",bfsl);
    
    /* Get elements msym elements */
    if(MSYM_SUCCESS != (ret = msymSetBasisFunctions(ctx, bfsl, bfs))) goto err;
    
    /* These are no longer needed, internal versions of these are kept in the context,
     * They are indexed in the same way that they have been allocated.
     * I.e. during orbital symmetrization or when getting the symmetrized LCAO,
     * the coefficients will correspond to the same indexing as "orbitals",
     * this is the main reason for the two levels of indirection */
    free(elements);  elements = NULL;
    
    /* Some trivial information */
    if(MSYM_SUCCESS != (ret = msymGetCenterOfMass(ctx,cm))) goto err;
    if(MSYM_SUCCESS != (ret = msymGetRadius(ctx,&radius))) goto err;
    
    printf("Molecule has center of mass [%lf; %lf; %lf] "
           "and a radius of %lf\n",cm[0],cm[1],cm[2],radius);
    
    /* Find molecular symmetry */
    if(MSYM_SUCCESS != (ret = msymFindSymmetry(ctx))) goto err;
    
    /* Get the point group name */
    if(MSYM_SUCCESS != (ret = msymGetPointGroupName(ctx, sizeof(char[6]), point_group))) goto err;
    if(MSYM_SUCCESS != (ret = msymGetSubgroups(ctx, &msgl, &msg))) goto err;
    printf("Found point group [0] %s select subgroup\n",point_group);
    for(int i = 0; i < msgl;i++) printf("\t [%d] %s\n",i+1,msg[i].name);
    int ssg = 0;
    printf("\nEnter point group[0-%d]:",msgl);
    while(scanf("%d", &ssg) <= 0 && ssg >= 0 && ssg <= msgl) printf("\nEnter point group[0-%d]:",msgl);
    if(ssg > 0){
        ssg--;
        printf("Selected point group %s\n",msg[ssg].name);
        if(MSYM_SUCCESS != (ret = msymSelectSubgroup(ctx, &msg[ssg]))) goto err;
        /* Retreive the symmetry operations again.
         * Everything has been rebuilt, and the old msops is no longer valid
         * Neither are the equivalence sets
         */
        if(MSYM_SUCCESS != (ret = msymGetSymmetryOperations(ctx, &msopsl, &msops))) goto err;
        
    }
    
    /* Set pointgroup to the C3v subgroup if it has XXX symmetry
     * using the same alignment as the original.
     * If specific axes are wanted the alignment axes/transform can be set using:
     * msym Get/Set Alignment Transform/Axes */
    if(0 == strncmp(point_group, "XXX", 3) && ssg == 0){
        //double transform[3][3];
        printf("Changing pointgroup from XXX -> C3v\n");
        //if(MSYM_SUCCESS != (ret = msymGetAlignmentTransform(ctx, transform))) goto err;
        if(MSYM_SUCCESS != (ret = msymSetPointGroupByType(ctx, MSYM_POINT_GROUP_TYPE_Cnv,3))) goto err;
        //if(MSYM_SUCCESS != (ret = msymSetAlignmentTransform(ctx, transform))) goto err;
        if(MSYM_SUCCESS != (ret = msymFindSymmetry(ctx))) goto err;
        if(MSYM_SUCCESS != (ret = msymGetPointGroupName(ctx, sizeof(char[6]), point_group))) goto err;
    }
    
    /* Retreive the symmetry operations */
    if(MSYM_SUCCESS != (ret = msymGetSymmetryOperations(ctx, &msopsl, &msops))) goto err;
    
    for(int i = 0; i < msopsl;i++){
        if(msops[i].type == MSYM_SYMMETRY_OPERATION_TYPE_PROPER_ROTATION && msops[i].order == 3 && msops[i].power == 1){
            printf("Found a C3^1 axis, YEY!\n");
        }
    }
    
    /* Aligning axes prior to orbital symmetrization will
     * change the orientation of orbitals with l >= 1 */
    if(MSYM_SUCCESS != (ret = msymAlignAxes(ctx))) goto err;
    
    /* Symmetrize the molecule.
     * You can do this before orbital symmetrization as well,
     * but the permutations are already built, so you don't need to */
    if(MSYM_SUCCESS != (ret = msymSymmetrizeElements(ctx, &symerr))) goto err;
    
    printf("Molecule has been symmetrized to point group %s "
           "with an error of %lf\n",point_group, symerr);
    
    if(MSYM_SUCCESS != (ret = msymGetElements(ctx, &mlength, &melements))) goto err;
    if(mlength != length){ printf("Not possible!\n"); goto err;}
    
    
    
    printf("New element coordinates:\n%d\n\n",mlength);
    for(int i = 0;i < mlength;i++){
        printf("%s %12.9lf %12.9lf %12.9lf\n",
               melements[i].name,
               melements[i].v[0],
               melements[i].v[1],
               melements[i].v[2]);
    }

    if(MSYM_SUCCESS != (ret = msymGetBasisFunctions(ctx, &mbfsl, &mbfs))) goto err;
    if(MSYM_SUCCESS != (ret = msymGetSubrepresentationSpaces(ctx, &msrsl, &msrs))) goto err;
    if(MSYM_SUCCESS != (ret = msymGetCharacterTable(ctx, &mct))) goto err;
    
    irrep = calloc(mct->d, sizeof(*irrep));
    
    for(int i = 0; i < msrsl;i++){
        printf("Got %d SALCs with %d partner functions of symmetry species %s\n",msrs[i].salcl,mct->s[msrs[i].s].d, mct->s[msrs[i].s].name);
        for(int j = 0;j < msrs[i].salcl;j++){
            char *type = "";
            msym_salc_t *salc = &msrs[i].salc[j];
            msym_basis_function_t *bf = salc->f[0];
            if(bf->type == MSYM_BASIS_TYPE_REAL_SPHERICAL_HARMONIC) type = "real spherical harmonic ";
            printf("\tSALC %d was constructed from %d %sbasis functions on %s with quantum numbers n=%d and l=%d\n",j,salc->fl,type,bf->element->name,bf->f.sh.n,bf->f.sh.l);
        }
    }
    
    if(MSYM_SUCCESS != (ret = msymGetSALCs(ctx, bfsl, salcs, species, pf))) goto err;
    
    printf("Species: ");
    for(int i = 0;i < bfsl;i++){
        printf("%s",mct->s[species[i]].name);
        if(i == bfsl - 1) printf("\n");
        else printf(", ");
    }

    
    /* Reorder the SALCs */
    for(int i = 0;i < bfsl;i++){
        int ni = i*i % bfsl;
        if(ni != i){
            memcpy(cmem, salcs[i], sizeof(double[bfsl]));
            memcpy(salcs[i], salcs[ni], sizeof(double[bfsl]));
            memcpy(salcs[ni], cmem, sizeof(double[bfsl]));
        }
    }
    
    /* Add some noise */
    srand((unsigned)time(NULL));
    for(int i = 0;i < bfsl;i++){
        for(int j = 0;j < bfsl;j++){
            double r = ((double) (rand() - (RAND_MAX >> 1)))/RAND_MAX;
            double x = i % 2 == 0 ? -1.0 : 1.0;
            salcs[i][j] += r*1.0e-5;
            salcs[i][j] *= x;
        }
    }
    
    /* Symmetrize wavefunctions */
    if(MSYM_SUCCESS != (ret = msymSymmetrizeWavefunctions(ctx, bfsl, salcs, species, pf))) goto err;
    
    for(int i = 0; i < bfsl;i++) cmem[i] = 1;
    if(MSYM_SUCCESS != (ret = msymSymmetrySpeciesComponents(ctx, bfsl, cmem, mct->d, irrep))) goto err;
    printf("Test wavefunction 1,1... components\n");
    printf("%lf%s", irrep[0], mct->s[0].name);
    for(int j = 1;j < mct->d;j++){
        printf(" + %lf%s", irrep[j], mct->s[j].name);
    }
    printf("\n");
    
    printf("Wave function symmetrization returned new linear combinations:\n");
    for(int i = 0;i < bfsl;i++){
        int s = species[i];
        printf("\t wf:%d is of symmetry species %s",i,mct->s[s].name);
        if(mct->s[s].d > 1){
            printf(" partner function %d to wf:%d",pf[i].d, pf[i].i);
        }
        printf("\n");
    }
    
    /* Make a new element with the same type as the first one we read */
    msym_element_t myelement;
    memset(&myelement,0,sizeof(msym_element_t));
    myelement.n = melements[0].n;
    myelement.v[0] = melements[0].v[0];
    myelement.v[1] = melements[0].v[1];
    myelement.v[2] = melements[0].v[2];
    
    /* Generate some new elements of the same point group */
    if(MSYM_SUCCESS != (ret = msymGenerateElements(ctx,1,&myelement))) goto err;
    
    /* This is not a memory leak, context keeps track of melements,
     * and it should never be freed, msymReleaseContext does this. */
    if(MSYM_SUCCESS != (ret = msymGetElements(ctx, &mlength, &melements))) goto err;
    
    printf("Generated element coordinates:\n%d\n\n",mlength);
    for(int i = 0;i < mlength;i++){
        printf("%s %lf %lf %lf\n",
               melements[i].name,
               melements[i].v[0],
               melements[i].v[1],
               melements[i].v[2]);
    }
    
    msymReleaseContext(ctx);
    printf("We're done!\n");
    
    free(salcs);
    free(cmem);
    free(species);
    free(pf);
    free(bfs);
    free(elements);
    free(irrep);
    
    return ret;
err:
    free(salcs);
    free(cmem);
    free(species);
    free(pf);
    free(bfs);
    free(elements);
    free(irrep);
    error = msymErrorString(ret);
    fprintf(stderr,"Error %s: ",error);
    error = msymGetErrorDetails();
    fprintf(stderr,"%s\n",error);
    return ret;
}

int read_xyz(const char *name, msym_element_t **ratoms) {
    FILE *fp = fopen(name,"r");
    msym_element_t *a;
    int l;
    char buf[1024];
    if(NULL == fp){
        fprintf(stderr, "could not open file: %s\n",name);
        return -1;
    }
    if (NULL == fgets(buf, sizeof(buf), fp) || sscanf(buf," %d ",&l) != 1){
        fprintf(stderr,"Unable to read file %s\n",name);
        fclose(fp);
        return -1;
    }
    if(l < 300000) {
        a = malloc(l*sizeof(msym_element_t));
        memset(a,0,l*sizeof(msym_element_t));
    } else {
        fprintf(stderr, "Too many elements in file %d\n",l);
        fclose(fp);
        return -1;
    }
    
    //char * fgets ( comment, sizeof(comment), fp );
    if(NULL != fgets(buf, sizeof(buf), fp)){
        printf("Comment: %.*s", (int)sizeof(buf), buf);
    }
    
    for (int i = 0; i < l && fgets(buf, sizeof(buf), fp) && sscanf(buf, "%s %lf %lf %lf", a[i].name, &(a[i].v[0]),  &(a[i].v[1]),  &(a[i].v[2])) == 4 && i < l; i++) {}
    *ratoms = a;
    fclose(fp);
    return l;
    
}
