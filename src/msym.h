//
//  msym.h
//  libmsym
//
//  Created by Marcus Johansson on 30/01/15.
//  Copyright (c) 2015 Marcus Johansson. 
//
//  Distributed under the MIT License ( See LICENSE file or copy at http://opensource.org/licenses/MIT )
//

#ifndef __MSYM_H
#define __MSYM_H

#ifdef __cplusplus
extern "C" {
#endif
    
#include "msym_error.h"
    
    typedef struct _msym_context * msym_context;

    typedef enum _msym_geometry {
        GEOMETRY_UNKNOWN = -1,
        SPHERICAL,
        LINEAR,
        PLANAR_REGULAR,
        PLANAR_IRREGULAR,
        POLYHEDRAL_PROLATE,
        POLYHEDRAL_OBLATE,
        ASSYMETRIC
    } msym_geometry_t;
    
    typedef enum _msym_basis_type {
        ATOMIC_ORBITAL,
        MASS_WEIGHTED_COORDINATES
    } msym_basis_type_t;
    
    typedef struct _msym_symmetry_operation {
        enum _msym_symmetry_operation_type {
            IDENTITY = 0,
            PROPER_ROTATION = 1,
            IMPROPER_ROTATION = 2,
            REFLECTION = 3,
            INVERSION = 4
        } type;
        int order;                              // Order of proper/improper rotation
        int power;                              // Power (e.g. C2^2 = I)
        enum _msym_symmetry_operation_orientation {
            NONE = -10,
            HORIZONTAL = -11,
            VERTICAL = -12,
            DIHEDRAL = -13
        } orientation;
        double v[3];                            // Proper/improper rotation vector or reflection plane normal
        int cla;                                // Class of symmetry operation (point group dependant)
    } msym_symmetry_operation_t;
    
    typedef enum _msym_point_group_type {
        POINT_GROUP_Ci,
        POINT_GROUP_Cs,
        POINT_GROUP_Cn,
        POINT_GROUP_Cnh,
        POINT_GROUP_Cnv,
        POINT_GROUP_Dn,
        POINT_GROUP_Dnh,
        POINT_GROUP_Dnd,
        POINT_GROUP_Sn,
        POINT_GROUP_T,
        POINT_GROUP_Td,
        POINT_GROUP_Th,
        POINT_GROUP_O,
        POINT_GROUP_Oh,
        POINT_GROUP_I,
        POINT_GROUP_Ih,
        POINT_GROUP_K,
        POINT_GROUP_Kh
    } msym_point_group_type_t;
    
    typedef struct _msym_subgroup {
        msym_point_group_type_t type;
        int n;
        int order;
        msym_symmetry_operation_t *primary;
        msym_symmetry_operation_t **sops;
        struct _msym_subgroup *subgroup[2];
        char name[6];
    } msym_subgroup_t;
    
    typedef struct _msym_thresholds {
        double zero;                            // For determining if something is zero (e.g. vectors close to center of mass)
        double geometry;                        // For translating inertial tensor eigenvalues to geometric structures
        double angle;                           // For determining angles, (e.g. if vectors are parallel)
        double equivalence;                     // Equivalence test threshold
        double eigfact;                         // Jacobi eigenvalue algorithm threshold
        double permutation;                     // Equality test when determining permutation for symmetry operation
        double orthogonalization;               // For orthogonalizing orbital subspaces
    } msym_thresholds_t;
    
    
    typedef struct _msym_element {
        void *id;                               // custom identifier
        double m;                               // Mass
        double v[3];                            // Position
        int n;                                  // Nuclear charge
        char name[4];                           // Name
    } msym_element_t;
    
    typedef struct _msym_equivalence_set {
        msym_element_t **elements;              // Pointers to elements
        double err;                             // Maximum error when detecting this equivalence set
        int length;                             // Number of elements
    } msym_equivalence_set_t ;
    
    typedef struct _msym_spherical_harmonic {
        int n;                                  // Principal
        int l;                                  // Azimuthal
        int m;                                  // Liniear combination of magnetic quantum number (e.g. 2pz = 0, 2px = 1, 2py = -1)
    } msym_spherical_harmonic_t;
    
    typedef struct _msym_basis_function {
        void *id;                               // custom identifier
        enum _msym_basis_type_t {
            REAL_SPHERICAL_HARMONIC,
            CARTESIAN
        } type;
        msym_element_t *element;
        union {
            msym_spherical_harmonic_t sh;       // Atomic orbital basis
        } f;
        char name[8];
    
    } msym_basis_function_t;
    
    typedef struct _msym_salc {
        void *pf;         // partner functions double[d][fl]
        int fl;             // number of basis functions
        int d;              // dimension of space (same as msym_character_table_t.s[msym_subspace_t.s].d)
        msym_basis_function_t **f;
    } msym_salc_t;
    
    typedef struct _msym_subspace_2 {
        int s, salcl; //symmetry species, irreps
        msym_salc_t *salc;
    } msym_subspace_2_t;
    
    //rename representation? we need info on wheather or not it's irreducible
    typedef struct _msym_symmetry_species {
        int d;
        char name[8];
    } msym_symmetry_species_t;
    
    typedef struct _msym_character_table {
        int d;
        int *classc;
        msym_symmetry_operation_t **sops;
        msym_symmetry_species_t *s;
        void *table;  //double[d][d]
    } msym_character_table_t;

    //There are better ways of representing a permutation (lika a Lehmer code) but I'll leave that for later
    typedef struct _msym_permutation_cycle_t {
        int l;
        int s;
    } msym_permutation_cycle_t;
    
    typedef struct _msym_permutation {
        int *p;
        int p_length;
        msym_permutation_cycle_t *c;
        int c_length;
    } msym_permutation_t;
    
    typedef struct {
        msym_point_group_type_t type;
        int n;
        int order;
        msym_symmetry_operation_t *primary;
        msym_symmetry_operation_t *sops;
        msym_permutation_t *perm;
        double transform[3][3];
        msym_character_table_t *ct;
        char name[8];
    } msym_point_group_t;
    
    
    msym_context msymCreateContext();
    msym_error_t msymReleaseContext(msym_context ctx);
    
    msym_error_t msymSetThresholds(msym_context ctx, msym_thresholds_t *thresholds);
    msym_error_t msymGetThresholds(msym_context ctx, msym_thresholds_t **thresholds);
    msym_error_t msymSetElements(msym_context ctx, int length, msym_element_t *elements);
    msym_error_t msymGetElements(msym_context ctx, int *length, msym_element_t **elements);
    msym_error_t msymSetBasisFunctions(msym_context ctx, int length, msym_basis_function_t *basis);
    msym_error_t msymGetBasisFunctions(msym_context ctx, int *length, msym_basis_function_t **basis);
    msym_error_t msymGetPointGroup(msym_context ctx, msym_point_group_t **pg);
    msym_error_t msymSetPointGroupByName(msym_context ctx, const char *name);
    msym_error_t msymGetPointGroupName(msym_context ctx, int l, char buf[l]);
    msym_error_t msymGetSubgroups(msym_context ctx, int *l, msym_subgroup_t **subgroups);
    msym_error_t msymSelectSubgroup(msym_context ctx, msym_subgroup_t *subgroup);
    msym_error_t msymGetSymmetryOperations(msym_context ctx, int *sopsl, msym_symmetry_operation_t **sops);
    msym_error_t msymGetEquivalenceSets(msym_context ctx, int *l, msym_equivalence_set_t **es);
    msym_error_t msymGetSALCSubspaces(msym_context ctx, int *l, msym_subspace_2_t **ss);
    msym_error_t msymGetCharacterTable(msym_context ctx, msym_character_table_t **ct);
    
    msym_error_t msymFindEquivalenceSets(msym_context ctx);
    msym_error_t msymFindEquivalenceSetPermutations(msym_context ctx);
    msym_error_t msymFindSymmetry(msym_context ctx);
    msym_error_t msymSymmetrizeMolecule(msym_context context, double *err);
    msym_error_t msymApplyTranslation(msym_context ctx, msym_element_t *element, double v[3]);
    msym_error_t msymSymmetrizeWavefunctions(msym_context ctx, int l, double c[l][l]);
    msym_error_t msymGenerateElements(msym_context ctx, int length, msym_element_t elements[length]);
    msym_error_t msymGenerateSALCSubspaces(msym_context ctx);
    msym_error_t msymAlignAxes(msym_context ctx);
    
    msym_error_t msymGetCenterOfMass(msym_context ctx, double v[3]);
    msym_error_t msymSetCenterOfMass(msym_context ctx, double v[3]);
    msym_error_t msymGetRadius(msym_context ctx, double *radius);
    msym_error_t msymGetGeometry(msym_context ctx, msym_geometry_t *geometry);
    msym_error_t msymGetPrincipalMoments(msym_context ctx, double eigval[3]);
    msym_error_t msymGetPrincipalAxes(msym_context ctx, double eigvec[3][3]);
    msym_error_t msymGetAlignmentAxes(msym_context ctx, double primary[3], double secondary[3]);
    msym_error_t msymSetAlignmentAxes(msym_context ctx, double primary[3], double secondary[3]);
    msym_error_t msymGetAlignmentTransform(msym_context ctx, double transform[3][3]);
    msym_error_t msymSetAlignmentTransform(msym_context ctx, double transform[3][3]);
    
#ifdef __cplusplus
}
#endif
    
#endif /* defined(__MSYM_H) */
