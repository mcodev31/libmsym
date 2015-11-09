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
#include "example.h"

int read_xyz(const char *name, msym_element_t **ratoms);

void printSALC(msym_salc_t *salc, msym_element_t *melements){

    
    double (*space)[salc->fl] = (double (*)[salc->fl]) salc->pf;
    for(int d = 0;d < salc->d;d++){
        if(salc->d > 1) printf("Component %d:\n",d+1);
        for(int line = 0; line < salc->fl; line+=6){
            for(int i = line;i < line + 6 && i < salc->fl;i++){
                msym_basis_function_t *bf = salc->f[i];
                printf(" %d%s %-8s\t",(int)(bf->element-melements)+1, bf->element->name,bf->name);
            }
            printf("\n");
            
            for(int i = line;i < line + 6 && i < salc->fl;i++){
                printf("%10.7lf\t", space[d][i]);
                
            }
            printf("\n\n");
        }
        printf("\n");
    }

}

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
    const msym_equivalence_set_t *mes = NULL;
    int mesl = 0;
    double *irrep = NULL;
    
    msym_basis_function_t *bfs = NULL;
    
    int msgl = 0, msopsl = 0, mlength = 0, msrsl = 0, mbfsl = 0, bfsl = 0;
        
    /* This function reads xyz files.
     * It initializes an array of msym_element_t to 0,
     * then sets the coordinates and name of the elements */
    int length = read_xyz(in_file, &elements);
    if(length <= 0) return -1;
    
    
    double (*psalcs)[bfsl] = NULL; // SALCs in matrix form, and input for symmetrization
    double *pcmem = NULL; // Some temporary memory
    int *pspecies = NULL;
    msym_partner_function_t *ppf = NULL;
    
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
    
    printf("Found point group [0] %s with %d subgroups:\n",point_group, msgl);
    printf("\t [0] %s\n\t -------\n", point_group);
    for(int i = 0; i < msgl;i++) printf("\t [%d] %s\n",i+1,msg[i].name);
    int ssg = 0;
    
    do {printf("\nChoose point group to use [0-%d]:",msgl);} while(scanf(" %d", &ssg) <= 0 || ssg < 0 || ssg > msgl);
    if(ssg > 0){
        ssg--;
        printf("Selected point group %s\n",msg[ssg].name);
        if(MSYM_SUCCESS != (ret = msymSelectSubgroup(ctx, &msg[ssg]))) goto err;
        if(MSYM_SUCCESS != (ret = msymGetPointGroupName(ctx, sizeof(char[6]), point_group))) goto err;
    }
    
    
    char yn = 'n';
    printf("\nWould you like to add basis functions? [y/N]:");
    if(scanf(" %c", &yn) > 0 && (yn | 0x60) == 'y'){
        int sel_es = 0;
        if(MSYM_SUCCESS != (ret = msymGetEquivalenceSets(ctx, &mesl, &mes))) goto err;
        do {
            
            int sel_n = 0, sel_l = 0, ele_bfsl = 0;
            printf("Basis functions can be added to entire molecule [0] or any of %d symmetrically equivalent sets:\n",mesl);
            printf("\t [0] Entire molecule\n");
            for(int i = 0; i < mesl;i++) printf("\t [%d] %d%s\n",i+1,mes[i].length, mes[i].elements[0]->name);
            
            do {printf("\nChoose set of elements [0-%d]:",mesl);} while(scanf(" %d", &sel_es) <= 0 || sel_es < 0 || sel_es > mesl);
            
            do {
                //z functions should be enough for an example
                printf("\nSelect principal quantum number (n) [1-21]:");
            } while (scanf(" %d", &sel_n) <= 0 || sel_n < 1 || sel_n > 21);
            
            do {
                printf("\nSelect angular momentum quantum number (l) [0-%d]:",sel_n-1);
            } while (scanf(" %d", &sel_l) <= 0 || sel_l < 0 || sel_l >= sel_n);
            
            
            ele_bfsl = (2*sel_l+1);
            
            if(sel_es > 0){
                const msym_equivalence_set_t *smes = &mes[sel_es-1];
                int nbfsl = smes->length*ele_bfsl;
                bfs = realloc(bfs, bfsl*sizeof(*bfs) + nbfsl*sizeof(*bfs));
                memset(&bfs[bfsl], 0, nbfsl*sizeof(*bfs));
                for(int i = 0;i < smes->length;i++){
                    for(int m = -sel_l;m <= sel_l;m++){
                        bfs[bfsl].element = smes->elements[i];
                        bfs[bfsl].type = MSYM_BASIS_TYPE_REAL_SPHERICAL_HARMONIC;
                        bfs[bfsl].f.rsh.n = sel_n;
                        bfs[bfsl].f.rsh.l = sel_l;
                        bfs[bfsl].f.rsh.m = m;
                        bfsl++;
                    }
                }
                printf("Will add %d real spherical harmonics basis functions with n=%d l=%d m=[%d,%d]to equivalence set %d\n", nbfsl, sel_n, sel_l, -sel_l, sel_l, sel_es);
            } else {
                int nbfsl = mlength*ele_bfsl;
                bfs = realloc(bfs, bfsl*sizeof(*bfs) + nbfsl*sizeof(*bfs));
                memset(&bfs[bfsl], 0, nbfsl*sizeof(*bfs));
                for(int i = 0;i < mlength;i++){
                    for(int m = -sel_l;m <= sel_l;m++){
                        // You can use either melements or your initial elements here, but melements is "more correct"
                        bfs[bfsl].element = &melements[i];
                        bfs[bfsl].type = MSYM_BASIS_TYPE_REAL_SPHERICAL_HARMONIC;
                        bfs[bfsl].f.rsh.n = sel_n;
                        bfs[bfsl].f.rsh.l = sel_l;
                        bfs[bfsl].f.rsh.m = m;
                        bfsl++;
                    }
                }
                printf("Will add %d real spherical harmonics basis functions with n=%d l=%d m=[%d,%d] to entire molecule\n", nbfsl, sel_n, sel_l, -sel_l, sel_l);
            }

            printf("\nWould you like to add more basis functions? [y/N]:");
        } while(scanf(" %c", &yn) > 0 && (yn | 0x60) == 'y');
        
        printf("Adding %d basis functions\n",bfsl);
        if(MSYM_SUCCESS != (ret = msymSetBasisFunctions(ctx, bfsl, bfs))) goto err;
        
        msym_point_group_type_t mtype;
        int mn;
        if(MSYM_SUCCESS != (ret = msymGetPointGroupType(ctx, &mtype, &mn))) goto err;
        
        /* Equivalence sets, subgroups and symmetry elements are updated when setting basis functions on linear groups.
         * Take care to update them (you can always update after setting basis functions, it's just a pointer update)
         */
        if((MSYM_POINT_GROUP_TYPE_Dnh == mtype || MSYM_POINT_GROUP_TYPE_Cnv == mtype) && 0 == mn){
            
            if(MSYM_SUCCESS != (ret = msymGetSubgroups(ctx, &msgl, &msg))) goto err;
            
            printf("Your selecton of basis functions resulted in new relevant subgroups\n");
            printf("Can now use [0] %s or any of %d subgroups:\n",point_group, msgl);
            printf("\t [0] %s\n\t -------\n", point_group);
            for(int i = 0; i < msgl;i++) printf("\t [%d] %s\n",i+1,msg[i].name);
            int ssg = 0;
            
            do {printf("\nChoose point group to use [0-%d]:",msgl);} while(scanf(" %d", &ssg) <= 0 || ssg < 0 || ssg > msgl);
            if(ssg > 0){
                ssg--;
                printf("Selected point group %s\n",msg[ssg].name);
                if(MSYM_SUCCESS != (ret = msymSelectSubgroup(ctx, &msg[ssg]))) goto err;
                if(MSYM_SUCCESS != (ret = msymGetPointGroupName(ctx, sizeof(char[6]), point_group))) goto err;
            }
            
            if(MSYM_SUCCESS != (ret = msymGetEquivalenceSets(ctx, &mesl, &mes))) goto err;
        }
        
        
    }
    
    /* Get elements msym elements */
    if(MSYM_SUCCESS != (ret = msymGetSymmetryOperations(ctx, &msopsl, &msops))) goto err;
    
    
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
    
    
    printf("\nWould you like to print the symmetry elements? [y/N]:");
    if(scanf(" %c", &yn) > 0 && (yn | 0x60) == 'y'){
        
        for(int i = 0; i < msopsl;i++){
            const msym_symmetry_operation_t *sop = &msops[i];
            char *rn = "";
            char *cn = "";
            switch(sop->orientation){
                case MSYM_SYMMETRY_OPERATION_ORIENTATION_HORIZONTAL : rn = "h"; break;
                case MSYM_SYMMETRY_OPERATION_ORIENTATION_VERTICAL   : rn = "v"; cn = "'"; break;
                case MSYM_SYMMETRY_OPERATION_ORIENTATION_DIHEDRAL   : rn = "d"; cn = "''"; break;
                default: break;
            }
            switch (sop->type) {
                case MSYM_SYMMETRY_OPERATION_TYPE_PROPER_ROTATION :
                    if(sop->order == 2) printf("C%d%s",sop->order,cn);
                    else printf("C%d%s^%d",sop->order,cn,sop->power);
                    printf(" around [%lf %lf %lf]\n", sop->v[0],sop->v[1],sop->v[2]);
                    break;
                case MSYM_SYMMETRY_OPERATION_TYPE_IMPROPER_ROTATION :
                    printf("S%d^%d",sop->order,sop->power);
                    printf(" around [%lf %lf %lf]\n", sop->v[0],sop->v[1],sop->v[2]);
                    break;
                case MSYM_SYMMETRY_OPERATION_TYPE_REFLECTION :
                    printf("R%s",rn);
                    printf(" with normal vector [%lf %lf %lf]\n", sop->v[0],sop->v[1],sop->v[2]);
                    break;
                case MSYM_SYMMETRY_OPERATION_TYPE_INVERSION :
                    printf("i\n");
                    break;
                case MSYM_SYMMETRY_OPERATION_TYPE_IDENTITY :
                    printf("E\n");
                    break;
                default  :
                    printf("?"); break;
            }
            
        }
    }
    
    if(bfsl > 0){
        // Just keeping memory in context of error handling
        double (*salcs)[bfsl] = psalcs = calloc(bfsl, sizeof(*salcs)); // SALCs in matrix form, and input for symmetrization
        double *cmem = pcmem = calloc(bfsl, sizeof(*cmem)); // Some temporary memory
        int *species = pspecies = calloc(bfsl, sizeof(*species));
        msym_partner_function_t *pf = ppf = calloc(bfsl, sizeof(*pf));
        
        if(MSYM_SUCCESS != (ret = msymGetBasisFunctions(ctx, &mbfsl, &mbfs))) goto err;
        if(MSYM_SUCCESS != (ret = msymGetSubrepresentationSpaces(ctx, &msrsl, &msrs))) goto err;
        if(MSYM_SUCCESS != (ret = msymGetCharacterTable(ctx, &mct))) goto err;
        
        irrep = calloc(mct->d, sizeof(*irrep));
        
        printf("\nGenerated SALCs from %d basis functions of %d symmetry species.\nWould you like to view them? [y/N]:", mbfsl, mct->d);
        if(scanf(" %c", &yn) > 0 && (yn | 0x60) == 'y'){
            do{
                
                int sel_ss = 0, sel_salc = 0;
                for(int i = 0; i < msrsl;i++){
                    if(msrs[i].salcl > 0)
                        printf("\t [%d] %s (%d SALCs with %d partner functions each)\n",i, mct->s[msrs[i].s].name, msrs[i].salcl, mct->s[msrs[i].s].d);
                    else
                        printf("\t [-] %s (no SALCs of this symmetry species)\n",mct->s[msrs[i].s].name);
                }
                do {printf("\nChoose subspace [0-%d]:",msrsl-1);} while(scanf(" %d", &sel_ss) <= 0 || sel_ss < 0 || sel_ss >= msrsl || msrs[sel_ss].salcl <= 0);
                
                int salcl = msrs[sel_ss].salcl;
                
                for(int i = 0; i < salcl;i++){
                    char *type = "";
                    msym_salc_t *salc = &msrs[sel_ss].salc[i];
                    msym_basis_function_t *bf = salc->f[0];
                    const msym_equivalence_set_t *salces = NULL;
                    if(MSYM_SUCCESS != (ret = msymGetEquivalenceSetByElement(ctx, bf->element, &salces))) goto err;
                    if(bf->type == MSYM_BASIS_TYPE_REAL_SPHERICAL_HARMONIC) type = "real spherical harmonic ";
                    printf("\t [%d] SALC constructed from %d %sbasis on equivalence set %d with n=%d and l=%d\n",i,salc->fl,type,(int) (salces - mes),bf->f.rsh.n,bf->f.rsh.l);
                        
                }
                
                do {printf("\nChoose SALC [0-%d]:",salcl-1);} while(scanf(" %d", &sel_salc) <= 0 || sel_salc < 0 || sel_salc > salcl);
                
                msym_salc_t *salc = &msrs[sel_ss].salc[sel_salc];
                
                printSALC(salc, melements);
                
                
                printf("\nWould you like to view more SALCs? [y/N]:");
            } while(scanf(" %c", &yn) > 0 && (yn | 0x60) == 'y');
            
        }
        
        printf("\nWould you like to test wavefunction symmetrization and component determination? [y/N]:");
        if(scanf(" %c", &yn) > 0 && (yn | 0x60) == 'y'){
            printf("Scrambling and adding noise to SALCS to use for symmetrization\n");
            if(MSYM_SUCCESS != (ret = msymGetSALCs(ctx, bfsl, salcs, species, pf))) goto err;
            
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
            
            printf("Wave function symmetrization returned new linear combinations:\n");
            for(int i = 0;i < bfsl;i++){
                int s = species[i];
                printf("\t wf:%d is of symmetry species %s",i,mct->s[s].name);
                if(mct->s[s].d > 1){
                    printf(" partner function %d to wf:%d",pf[i].d, pf[i].i);
                }
                printf("\n");
            }
            
            for(int i = 0; i < bfsl;i++) cmem[i] = 1;
            
            if(MSYM_SUCCESS != (ret = msymSymmetrySpeciesComponents(ctx, bfsl, cmem, mct->d, irrep))) goto err;
            
            printf("Test wavefunction 1,1... components\n");
            printf("%lf%s", irrep[0], mct->s[0].name);
            for(int j = 1;j < mct->d;j++){
                printf(" + %lf%s", irrep[j], mct->s[j].name);
            }
            printf("\n");
            
        }
        
        
    }
    int symprint = 0;
    printf("\nWould you like to symmetrize the molecule? [y/N]:");
    if(scanf(" %c", &yn) > 0 && (yn | 0x60) == 'y'){
        symprint = 1;
        /* Symmetrize the molecule.
         * You can do this before orbital symmetrization as well,
         * but the permutations are already built, so you don't need to */
        if(MSYM_SUCCESS != (ret = msymSymmetrizeElements(ctx, &symerr))) goto err;
        printf("Molecule has been symmetrized to point group %s "
               "with an error of %lf\n",point_group, symerr);
    }
    
    printf("\nWould you like to align it to the xyz axis? [y/N]:");
    if(scanf(" %c", &yn) > 0 && (yn | 0x60) == 'y'){
        symprint = 1;
        /* Aligning axes prior to orbital symmetrization will
         * change the orientation of orbitals with l >= 1
         * relative to the molecular orientation */
        if(MSYM_SUCCESS != (ret = msymAlignAxes(ctx))) goto err;
        printf("Molecule has been aligned to the xyz axes\n");
    }
    
    if(symprint){
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
    }
    
    /* Make a new element with the same type as the first one we read */
    msym_element_t myelement;
    memset(&myelement,0,sizeof(msym_element_t));
    myelement.n = melements[0].n;
    myelement.v[0] = melements[0].v[0];
    myelement.v[1] = melements[0].v[1];
    myelement.v[2] = melements[0].v[2];
    
    printf("\nWould you like generate a new molecule based on %s at [%lf %lf %lf]? [y/N]:", melements[0].name, melements[0].v[0], melements[0].v[1], melements[0].v[2]);
    if(scanf(" %c", &yn) > 0 && (yn | 0x60) == 'y'){
        
        /* Generate some new elements of the same point group */
        if(MSYM_SUCCESS != (ret = msymGenerateElements(ctx,1,&myelement))) goto err;
        
        /* This is not a memory leak, context keeps track of melements,
         * and it should never be freed, msymReleaseContext does this. */
        if(MSYM_SUCCESS != (ret = msymGetElements(ctx, &mlength, &melements))) goto err;
        
        printf("\n\nGenerated element coordinates based on %s at [%lf %lf %lf]:\n%d\n\n",melements[0].name, melements[0].v[0], melements[0].v[1], melements[0].v[2], mlength);
        for(int i = 0;i < mlength;i++){
            printf("%s %12.9lf %12.9lf %12.9lf\n",
                   melements[i].name,
                   melements[i].v[0],
                   melements[i].v[1],
                   melements[i].v[2]);
        }
    }
    
    msymReleaseContext(ctx);
    printf("We're done!\n");
    
    free(psalcs);
    free(pcmem);
    free(pspecies);
    free(ppf);
    free(bfs);
    free(elements);
    free(irrep);
    
    return ret;
err:
    free(psalcs);
    free(pcmem);
    free(pspecies);
    free(ppf);
    free(bfs);
    free(elements);
    free(irrep);
    error = msymErrorString(ret);
    fprintf(stderr,"Error %s: ",error);
    error = msymGetErrorDetails();
    fprintf(stderr,"%s\n",error);
    msymReleaseContext(ctx);
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
