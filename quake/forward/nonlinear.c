/* -*- C -*- */

/* @copyright_notice_start
 *
 * This file is part of the CMU Hercules ground motion simulator developed
 * by the CMU Quake project.
 *
 * Copyright (C) Carnegie Mellon University. All rights reserved.
 *
 * This program is covered by the terms described in the 'LICENSE.txt' file
 * included with this software package.
 *
 * This program comes WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * 'LICENSE.txt' file for more details.
 *
 *  @copyright_notice_end
 */

#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <gsl/gsl_poly.h>
#include <stdio.h>

#include "geometrics.h"
#include "nonlinear.h"
#include "octor.h"
#include "psolve.h"
#include "quake_util.h"
#include "util.h"

#define  QC  qc = 0.577350269189 /* sqrt(3.0)/3.0; */

#define MAX(a, b) ((a)>(b)?(a):(b))
#define MIN(a, b) ((a)<(b)?(a):(b))

#define  XI  xi[3][8] = { {-1,  1, -1,  1, -1,  1, -1, 1} , \
                          {-1, -1,  1,  1, -1, -1,  1, 1} , \
                          {-1, -1, -1, -1,  1,  1,  1, 1} }

static int superflag = 0;

/* -------------------------------------------------------------------------- */
/*                             Global Variables                               */
/* -------------------------------------------------------------------------- */

static nlsolver_t           *myNonlinSolver;
static int32_t               thePropertiesCount;
static materialmodel_t       theMaterialModel;
static plasticitytype_t      thePlasticityModel;
static noyesflag_t           theApproxGeoState  = NO;
static noyesflag_t           theTensionCutoff   = NO;
static double               *theVsLimits;
static double               *theAlphaCohes;
static double               *theKayPhis;
static double               *theStrainRates;
static double               *theSensitivities;
static double               *theHardeningModulus;
static double               *theBetaDilatancy;
static double               *theGamma0;
static double               *thePsi;
static double               *theM;
static double               *theTheta1;
static double               *theTheta2;
static double               *theTheta3;
static double               *theTheta4;
static double               *theTheta5;
static double               *theTau_y;
static double                theGeostaticLoadingT = 0;
static double                theErrorTol          = 1E-03;
static double                theGeostaticCushionT = 0;
static int                   theGeostaticFinalStep;
static int                   theNoSubsteps=1000;
static int32_t              *myStationsElementIndices;
//static nlstation_t          *myNonlinStations;
static int32_t              *myNonlinStationsMapping;
static int32_t               myNumberOfNonlinStations;
static int32_t               myNonlinElementsCount = 0;
static int32_t              *myNonlinElementsMapping;
static int32_t               myBottomElementsCount = 0;
static bottomelement_t      *myBottomElements;
static int32_t               theNonlinearFlag = 0;

static double totalWeight = 0;

/* -------------------------------------------------------------------------- */
/*                                 Utilities                                  */
/* -------------------------------------------------------------------------- */

double get_geostatic_total_time() {
    return theGeostaticLoadingT + theGeostaticCushionT;
}

/*
 * Return YES if an element is to be considered nonlinear, NO otherwise.
 */
int isThisElementNonLinear(mesh_t *myMesh, int32_t eindex) {

	elem_t  *elemp;
	edata_t *edata;

	if ( theNonlinearFlag == 0 )
		return NO;

	elemp = &myMesh->elemTable[eindex];
	edata = (edata_t *)elemp->data;

	if ( ( edata->Vs <=  theVsLimits[thePropertiesCount-1] ) && ( edata->Vs >=  theVsLimits[0] ) )
		return YES;

	return NO;
}

/*
 * Performs linear interpolation between to given vectors for a given value
 * and gives back the corresponding pair for the requested Vs.
 *
 */
double interpolate_property_value(double vsRequest, double *propVector)
{
    int i;
    double vs0, vs1;
    double prop0, prop1;

    /* below floor case */
    if ( vsRequest <= theVsLimits[0] ) {
        return propVector[0];
    }

    /* above ceiling case */
    if ( vsRequest > theVsLimits[thePropertiesCount-1] ) {
        return propVector[thePropertiesCount-1];
    }

    /* in between cases */

    for ( i = 0; i < thePropertiesCount-1; i++ ) {

        if ( ( vsRequest >  theVsLimits[i]   ) &&
             ( vsRequest <= theVsLimits[i+1] ) ) {

            vs0 = theVsLimits[i];
            vs1 = theVsLimits[i+1];

            prop0 = propVector[i];
            prop1 = propVector[i+1];
        }
    }

    double result = prop0 + (vsRequest - vs0)*( (prop1 - prop0)/(vs1 - vs0) );

    return result;
}

/*
 * Returns the value of the constant alpha for Drucker-Prager's material
 * model
 */
double get_alpha(double vs, double phi) {

    double alpha;
    alpha = 2. * sin(phi) / ( sqrt(3.0) * ( 3. - sin(phi) ) );
    if ( alpha > 0.40 ) {
    	fprintf(stderr,"Illegal alpha= %f "
    			"Friction angle larger that 50deg\n", alpha);
    	MPI_Abort(MPI_COMM_WORLD, ERROR);
    	exit(1);
    }

    return alpha;
}


double get_gamma(double vs, double phi) {

    double gamma;

    gamma     = 6. * cos(phi) / ( sqrt(3.0) * ( 3. - sin(phi) ) );

    return gamma;
}

/*
 * Returns the value of the constant phi  in Drucker-Prager's material
 * model
 */
double get_phi(double vs) {

	double phi;

	phi   = interpolate_property_value(vs, theKayPhis) * PI / 180.0;

	return phi;
}

/*
 * Returns the value of the constant phi  in Drucker-Prager's material
 * model
 */
double get_dilatancy(double vs) {

	double dilt;

	dilt   = interpolate_property_value(vs, theBetaDilatancy) * PI / 180.0;

	return dilt;
}

/*
 * Returns the value of the constant beta (dilatancy angle)  in Drucker-Prager's material
 * model.
 */
double get_beta(double vs) {

	double  beta, dil;

	dil   = interpolate_property_value(vs, theBetaDilatancy) * PI / 180.0;
	beta  = 2. * sin(dil) / ( sqrt(3.0) * ( 3. - sin(dil) ) );

	return beta;
}


double get_cohesion(double vs) {

	double  coh;
	coh   = interpolate_property_value(vs, theAlphaCohes);
	return coh;
}

double get_hardmod(double vs) {

	double  hrd;
	hrd   = interpolate_property_value(vs, theHardeningModulus);

	return hrd;
}

double get_gamma_eff (double vs30, double zo)  {

	double gamma_eff=0;

	if ( ( vs30 >= 760 ) && ( vs30 <= 1500 ) )  { /* Site class B "Rock"  */
		if ( ( zo >= 0 ) && ( zo <= 6.096 ) ) {
			gamma_eff = pow(10.0,-1.73) / 100.0;
		} else if ( ( zo > 6.096 ) && ( zo <= 15.24 ) )   {
			gamma_eff = pow(10.0,-1.66) / 100.0;
		} else if ( ( zo > 15.24 ) && ( zo <= 36.576 ) )  {
			gamma_eff = pow(10.0,-1.58) / 100.0;
		} else if ( ( zo > 36.576 ) && ( zo <= 76.20 ) )  {
			gamma_eff = pow(10.0,-1.45) / 100.0;
		} else if ( ( zo > 76.20 ) && ( zo <= 152.40 ) )  {
			gamma_eff = pow(10.0,-1.39) / 100.0;
		} else if ( ( zo > 152.40 ) && ( zo <= 304.80 ) ) {
			gamma_eff = pow(10.0,-1.317) / 100.0;
		} else { /* This should not happen */
	    	fprintf(stderr, "nonlinear error for site class B. Attempting to get gamma_eff for vs30= %f at z= %f \n", vs30, zo);
	    	MPI_Abort(MPI_COMM_WORLD, ERROR);
	    	exit(1);
	    }
	} else if ( ( vs30 >= 200 ) && ( vs30 < 760 ) ) { /* Site class C & D "Sand"  */

		if ( ( zo >= 0 ) && ( zo <= 6.096 ) ) {
			gamma_eff = pow(10.0,-1.49) / 100.0;
		} else if ( ( zo > 6.096 ) && ( zo <= 15.24 ) )   {
			gamma_eff = pow(10.0,-1.29) / 100.0;
		} else if ( ( zo > 15.24 ) && ( zo <= 36.576 ) )  {
			gamma_eff = pow(10.0,-1.14) / 100.0;
		} else if ( ( zo > 36.576 ) && ( zo <= 76.20 ) )  {
			gamma_eff = pow(10.0,-1.00) / 100.0;
		} else if ( ( zo > 76.20 ) && ( zo <= 152.40 ) )  {
			gamma_eff = pow(10.0,-0.87) / 100.0;
		} else if ( ( zo > 152.40 ) && ( zo <= 304.80 ) ) {
			gamma_eff = pow(10.0,-0.7) / 100.0;
		} else { /* This should not happen */
	    	fprintf(stderr, "nonlinear error for site class C and D.  Attempting to get gamma_eff for vs30= %f at z= %f \n", vs30, zo);
	    	MPI_Abort(MPI_COMM_WORLD, ERROR);
	    	exit(1);
	    }
	} else {

			fprintf(stderr, "nonlinear error attempting to get gamma_eff - vs30: %f or z= %f are out of range \n", vs30, zo);
			MPI_Abort(MPI_COMM_WORLD, ERROR);
			exit(1);
	}

return gamma_eff;
}

/* -------------------------------------------------------------------------- */
/*       Initialization of parameters, structures and memory allocations      */
/* -------------------------------------------------------------------------- */



void nonlinear_init( int32_t     myID,
                     const char *parametersin,
                     double      theDeltaT,
                     double      theEndT )
{
    double  double_message[3];
    int     int_message[8];

    /* Capturing data from file --- only done by PE0 */
    if (myID == 0) {
        if (nonlinear_initparameters(parametersin, theDeltaT, theEndT) != 0) {
            fprintf(stderr,"Thread 0: nonlinear_local_init: "
                    "nonlinear_initparameters error\n");
            MPI_Abort(MPI_COMM_WORLD, ERROR);
            exit(1);
        }
    }

    /* Broadcasting data */
    double_message[0] = theGeostaticLoadingT;
    double_message[1] = theGeostaticCushionT;
    double_message[2] = theErrorTol;

    int_message[0] = (int)theMaterialModel;
    int_message[1] = thePropertiesCount;
    int_message[2] = theGeostaticFinalStep;
    int_message[3] = (int)thePlasticityModel;
    int_message[4] = (int)theApproxGeoState;
    int_message[5] = (int)theNonlinearFlag;
    int_message[6] = (int)theTensionCutoff;
    int_message[7] = (int)theNoSubsteps;

    MPI_Bcast(double_message, 3, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(int_message,    8, MPI_INT,    0, comm_solver);

    theGeostaticLoadingT  = double_message[0];
    theGeostaticCushionT  = double_message[1];
    theErrorTol           = double_message[2];

    theMaterialModel      = int_message[0];
    thePropertiesCount    = int_message[1];
    theGeostaticFinalStep = int_message[2];
    thePlasticityModel    = int_message[3];
    theApproxGeoState     = int_message[4];
    theNonlinearFlag      = int_message[5];
    theTensionCutoff      = int_message[6];
    theNoSubsteps         = int_message[7];

    /* allocate table of properties for all other PEs */

    if (myID != 0) {
        theVsLimits         = (double*)malloc(sizeof(double) * thePropertiesCount);
        theAlphaCohes       = (double*)malloc(sizeof(double) * thePropertiesCount);
        theKayPhis          = (double*)malloc(sizeof(double) * thePropertiesCount);
        theStrainRates      = (double*)malloc(sizeof(double) * thePropertiesCount);
        theSensitivities    = (double*)malloc(sizeof(double) * thePropertiesCount);
        theHardeningModulus = (double*)malloc(sizeof(double) * thePropertiesCount);
        theBetaDilatancy    = (double*)malloc(sizeof(double) * thePropertiesCount);
        theGamma0           = (double*)malloc(sizeof(double) * thePropertiesCount);
        thePsi              = (double*)malloc(sizeof(double) * thePropertiesCount);
        theM                = (double*)malloc(sizeof(double) * thePropertiesCount);
        theTheta1	 	    = (double*)malloc( sizeof(double) * thePropertiesCount);
        theTheta2	 	    = (double*)malloc( sizeof(double) * thePropertiesCount);
        theTheta3	 	    = (double*)malloc( sizeof(double) * thePropertiesCount);
        theTheta4	 	    = (double*)malloc( sizeof(double) * thePropertiesCount);
        theTheta5	 	    = (double*)malloc( sizeof(double) * thePropertiesCount);
        theTau_y 	 	    = (double*)malloc( sizeof(double) * thePropertiesCount);
    }

    /* Broadcast table of properties */
    MPI_Bcast(theVsLimits,         thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(theAlphaCohes,       thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(theKayPhis,          thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(theStrainRates,      thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(theSensitivities,    thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(theHardeningModulus, thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(theBetaDilatancy,    thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(theGamma0,           thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(thePsi,              thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(theM,                thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(theTheta1,           thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(theTheta2,           thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(theTheta3,           thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(theTheta4,           thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(theTheta5,           thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(theTau_y,            thePropertiesCount, MPI_DOUBLE, 0, comm_solver);
}

/*
 * Reads from parameters.in and stores in PE0 globals.
 */
int32_t nonlinear_initparameters ( const char *parametersin,
                                   double      theDeltaT,
                                   double      theEndT )
{
    FILE    *fp;
    int32_t  properties_count, no_substeps;
    int      row;
    double   geostatic_loading_t, geostatic_cushion_t, errorTol,
            *auxiliar;
    char     material_model[64],
             plasticity_type[64], approx_geostatic_state[64], tension_cutoff[64];

    materialmodel_t      materialmodel;
    plasticitytype_t     plasticitytype;
    noyesflag_t          approxgeostatic = -1;
    noyesflag_t          tensioncutoff = -1;

    /* Opens numericalin file */
    if ((fp = fopen(parametersin, "r")) == NULL) {
        fprintf(stderr, "Error opening %s\n at nl_init_parameters", parametersin);
        return -1;
    }

    /* Parses parameters.in to capture nonlinear single-value parameters */

    if ( (parsetext(fp, "geostatic_loading_time_sec",   'd', &geostatic_loading_t    ) != 0) ||
         (parsetext(fp, "geostatic_cushion_time_sec",   'd', &geostatic_cushion_t    ) != 0) ||
         (parsetext(fp, "error_tolerance",              'd', &errorTol               ) != 0) ||
         (parsetext(fp, "material_model",               's', &material_model         ) != 0) ||
         (parsetext(fp, "approximate_geostatic_state",  's', &approx_geostatic_state ) != 0) ||
         (parsetext(fp, "material_plasticity_type",     's', &plasticity_type        ) != 0) ||
         (parsetext(fp, "material_properties_count",    'i', &properties_count       ) != 0) ||
         (parsetext(fp, "no_substeps",                  'i', &no_substeps            ) != 0) ||
         (parsetext(fp, "tension_cutoff",               's', &tension_cutoff         ) != 0) )
    {
        fprintf(stderr, "Error parsing nonlinear parameters from %s\n", parametersin);
        return -1;
    }

    /* Performs sanity checks */
    if ( (geostatic_loading_t < 0) || (geostatic_cushion_t < 0) ||
         (geostatic_loading_t + geostatic_cushion_t > theEndT) ) {
        fprintf(stderr, "Illegal geostatic loading/cushion time %f/%f\n",
                geostatic_loading_t, geostatic_cushion_t);
        return -1;
    }

    if ( strcasecmp(material_model, "linear") == 0 ) {
        materialmodel = LINEAR;
    } else if ( strcasecmp(material_model, "vonMises_ep") == 0 )   {
        materialmodel = VONMISES_EP;
    } else if ( strcasecmp(material_model, "vonMises_fa") == 0 )   {
        materialmodel = VONMISES_FA;
    } else if ( strcasecmp(material_model, "vonMises_faM") == 0 )  {
        materialmodel = VONMISES_FAM;
    }  else if ( strcasecmp(material_model, "vonMises_baE") == 0 ) {
        materialmodel = VONMISES_BAE;
    }  else if ( strcasecmp(material_model, "vonMises_baH") == 0 ) {
        materialmodel = VONMISES_BAH;
    } else if ( strcasecmp(material_model, "vonMises_GQH") == 0 )  {
        materialmodel = VONMISES_GQH;
    } else if ( strcasecmp(material_model, "vonMises_MKZ") == 0 )  {
        materialmodel = VONMISES_MKZ;
    } else if ( strcasecmp(material_model, "vonMises_RO") == 0 )  {
        materialmodel = VONMISES_RO;
    } else if ( strcasecmp(material_model, "MohrCoulomb") == 0 )   {
        materialmodel = MOHR_COULOMB;
    } else if ( strcasecmp(material_model, "DruckerPrager") == 0 ) {
        materialmodel = DRUCKERPRAGER;
    }
    else {
        fprintf(stderr,
                "Illegal material model for nonlinear analysis \n"
                "(linear, vonMises_ep (Elasto-plastic), vonMises_FA (Frederick-Armstrong), vonMises_FAM (Frederick-Armstrong modified), \n "
                "vonMises_BAE (Borja-Aimes exponential), vonMises_BAH (Borja-Aimes hyperbolic), vonMises_GQH (Groholski et al GQH model), \n"
        		"vonMises_MKZ (Matasovic 1994), vonMises_RO (Ramberg-Osgood), DruckerPrager, MohrCoulomb): %s\n", material_model);
        return -1;
    }

    if (properties_count < 1) {
        fprintf(stderr, "Illegal material properties count %d\n", properties_count);
        return -1;
    }

    if ( strcasecmp(plasticity_type, "rate_dependant") == 0 ) {
        plasticitytype = RATE_DEPENDANT;
    } else if ( strcasecmp(plasticity_type, "rate_independant") == 0 ) {
        plasticitytype = RATE_INDEPENDANT;
    } else {
        fprintf(stderr,
                "Illegal material plasticity type for nonlinear "
                "analysis (rate_dependant, rate_independant): %s\n",
                plasticity_type);
        return -1;
    }

    if ( strcasecmp(approx_geostatic_state, "yes") == 0 ) {
        approxgeostatic = YES;
    } else if ( strcasecmp(approx_geostatic_state, "no") == 0 ) {
    	approxgeostatic = NO;
    } else {
        fprintf(stderr,
        		":Unknown response for considering an "
                "approximate geostatic state (yes or no): %s\n",
                approx_geostatic_state );
        return -1;
    }

    if ( ( (geostatic_loading_t > 0) || (geostatic_cushion_t > 0) ) &&
         ( approxgeostatic == YES ) ) {
        fprintf(stderr, "Approximate geostatic-state must be set to (no) when geostatic loading/cushion time > 0. %s\n",
        		approx_geostatic_state);
        return -1;
    }

    if ( ( strcasecmp(tension_cutoff, "yes") == 0 ) && ( ( materialmodel == DRUCKERPRAGER ) || ( materialmodel == MOHR_COULOMB ) ) ) {
    	tensioncutoff = YES;
    } else if ( ( strcasecmp(tension_cutoff, "yes") == 0 ) && ( ( materialmodel != DRUCKERPRAGER ) || ( materialmodel != MOHR_COULOMB ) ) ) {
    	fprintf(stderr,
    			":Tension cutoff option available "
    			"only for Mohr-Coulomb or Drucker-Prager models: %s\n",
    			material_model );
    	return -1;
    } else if ( strcasecmp(tension_cutoff, "no") == 0 ) {
    	tensioncutoff = NO;
    } else {
    	fprintf(stderr,
    			":Unknown response for considering "
    			"tension cutoff (yes or no): %s\n",
    			tension_cutoff );
    	return -1;
    }


    /* Initialize the static global variables */
    theGeostaticLoadingT  = geostatic_loading_t;
    theGeostaticCushionT  = geostatic_cushion_t;
    theErrorTol           = errorTol;
    theGeostaticFinalStep = (int)( (geostatic_loading_t + geostatic_cushion_t) / theDeltaT );
    theMaterialModel      = materialmodel;
    thePropertiesCount    = properties_count;
    thePlasticityModel    = plasticitytype;
    theApproxGeoState     = approxgeostatic;
    theTensionCutoff      = tensioncutoff;
    theNoSubsteps         = no_substeps;

    auxiliar             = (double*)malloc( sizeof(double) * thePropertiesCount * 16 );
    theVsLimits          = (double*)malloc( sizeof(double) * thePropertiesCount );
    theAlphaCohes        = (double*)malloc( sizeof(double) * thePropertiesCount );
    theKayPhis           = (double*)malloc( sizeof(double) * thePropertiesCount );
    theStrainRates       = (double*)malloc( sizeof(double) * thePropertiesCount );
    theSensitivities     = (double*)malloc( sizeof(double) * thePropertiesCount );
    theHardeningModulus  = (double*)malloc( sizeof(double) * thePropertiesCount );
    theBetaDilatancy     = (double*)malloc( sizeof(double) * thePropertiesCount );
    theGamma0            = (double*)malloc( sizeof(double) * thePropertiesCount );
    thePsi 				 = (double*)malloc( sizeof(double) * thePropertiesCount );
    theM 				 = (double*)malloc( sizeof(double) * thePropertiesCount );
    theTheta1	 	     = (double*)malloc( sizeof(double) * thePropertiesCount );
    theTheta2	 	     = (double*)malloc( sizeof(double) * thePropertiesCount );
    theTheta3	 	     = (double*)malloc( sizeof(double) * thePropertiesCount );
    theTheta4	 	     = (double*)malloc( sizeof(double) * thePropertiesCount );
    theTheta5	 	     = (double*)malloc( sizeof(double) * thePropertiesCount );
    theTau_y	 	     = (double*)malloc( sizeof(double) * thePropertiesCount );


    if ( parsedarray( fp, "material_properties_list", thePropertiesCount * 16, auxiliar ) != 0) {
        fprintf(stderr, "Error parsing nonlinear material properties list from %s\n", parametersin);
        return -1;
    }

    for ( row = 0; row < thePropertiesCount; row++) {
        theVsLimits[row]          = auxiliar[ row * 16     ];
        theAlphaCohes[row]        = auxiliar[ row * 16 + 1 ];
        theKayPhis[row]           = auxiliar[ row * 16 + 2 ];
        theStrainRates[row]       = auxiliar[ row * 16 + 3 ];
        theSensitivities[row]     = auxiliar[ row * 16 + 4 ];
        theHardeningModulus[row]  = auxiliar[ row * 16 + 5 ];
        theBetaDilatancy[row]     = auxiliar[ row * 16 + 6 ];
        theGamma0[row]            = auxiliar[ row * 16 + 7 ];
        thePsi[row]               = auxiliar[ row * 16 + 8 ];
        theM[row]                 = auxiliar[ row * 16 + 9 ];
        theTheta1[row]            = auxiliar[ row * 16 + 10 ];
        theTheta2[row]            = auxiliar[ row * 16 + 11 ];
        theTheta3[row]            = auxiliar[ row * 16 + 12 ];
        theTheta4[row]            = auxiliar[ row * 16 + 13 ];
        theTheta5[row]            = auxiliar[ row * 16 + 14 ];
        theTau_y [row]            = auxiliar[ row * 16 + 15 ];

    }

    theNonlinearFlag = 1;

    return 0;
}

/*
 * Counts the number of nonlinear elements in my local mesh
 */
void nonlinear_elements_count(int32_t myID, mesh_t *myMesh) {

    int32_t eindex;
    int32_t count = 0;

    for (eindex = 0; eindex < myMesh->lenum; eindex++) {

        if ( isThisElementNonLinear(myMesh, eindex) == YES ) {
            count++;
        }
    }

    if ( count > myMesh-> lenum ) {
        fprintf(stderr,"Thread %d: nl_elements_count: "
                "more elements than expected\n", myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }

    myNonlinElementsCount = count;

    return;
}

/*
 * Re-counts and stores the nonlinear element indices to a static local array
 * that will serve as mapping tool to the local mesh elements table.
 */
void nonlinear_elements_mapping(int32_t myID, mesh_t *myMesh) {

    int32_t eindex;
    int32_t count = 0;

    XMALLOC_VAR_N(myNonlinElementsMapping, int32_t, myNonlinElementsCount);

    for (eindex = 0; eindex < myMesh->lenum; eindex++) {

        if ( isThisElementNonLinear(myMesh, eindex) == YES ) {
            myNonlinElementsMapping[count] = eindex;
            count++;
        }
    }

    if ( count != myNonlinElementsCount ) {
        fprintf(stderr,"Thread %d: nl_elements_mapping: "
                "more elements than the count\n", myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }

    return;
}

noyesflag_t isThisElementsAtTheBottom( mesh_t  *myMesh,
                                       int32_t  eindex,
                                       double   depth )
{
    elem_t  *elemp;
    int32_t  nindex;
    double   z_m;

    /* Capture the element's last node at the bottom */
    elemp  = &myMesh->elemTable[eindex];
    nindex = elemp->lnid[7];

    z_m = (myMesh->ticksize)*(double)myMesh->nodeTable[nindex].z;

    if ( z_m == depth ) {
        return YES;
    }

    return NO;
}

void bottom_elements_count(int32_t myID, mesh_t *myMesh, double depth ) {

    int32_t eindex;
    int32_t count = 0;

    for (eindex = 0; eindex < myMesh->lenum; eindex++) {

        if ( isThisElementsAtTheBottom(myMesh, eindex, depth) == YES ) {
            count++;
        }
    }

    if ( count > myMesh-> lenum ) {
        fprintf(stderr,"Thread %d: bottom_elements_count: "
                "more elements than expected\n", myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }

    myBottomElementsCount = count;

    return;
}

void bottom_elements_mapping(int32_t myID, mesh_t *myMesh, double depth) {

    int32_t eindex;
    int32_t count = 0;

    XMALLOC_VAR_N(myBottomElements, bottomelement_t, myBottomElementsCount);

    for (eindex = 0; eindex < myMesh->lenum; eindex++) {

        if ( isThisElementsAtTheBottom(myMesh, eindex, depth) == YES ) {
            myBottomElements[count].element_id = eindex;
            count++;
        }
    }

    if ( count != myBottomElementsCount ) {
        fprintf(stderr,"Thread %d: bottom_elements_mapping: "
                "more elements than expected\n", myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }

    return;
}

/*
 * Prints statistics about number of nonlinear elements and stations in nonlinear
 * elements.
 */
void nonlinear_print_stats(int32_t *nonlinElementsCount,
                           int32_t *nonlinStationsCount,
                           int32_t *bottomElementsCount,
                           int32_t  theGroupSize)
{

    int pid;
    global_id_t totalElements = 0;
    global_id_t totalStations = 0;
    global_id_t totalBottom   = 0;

    FILE *fp = hu_fopen( "stat-nonlin.txt", "w" );

    fputs( "\n"
           "# ---------------------------------------- \n"
           "# Nonlinear elements and stations count:   \n"
           "# ---------------------------------------- \n"
           "# Rank    Elements    Stations      Bottom \n"
           "# ---------------------------------------- \n", fp );

    for ( pid = 0; pid < theGroupSize; pid++ ) {

        fprintf( fp, "%06d %11d %11d %11d\n", pid,
                 nonlinElementsCount[pid],
                 nonlinStationsCount[pid],
                 bottomElementsCount[pid] );

        totalElements += nonlinElementsCount[pid];
        totalStations += nonlinStationsCount[pid];
        totalBottom   += bottomElementsCount[pid];
    }

    fprintf( fp,
             "# ---------------------------------------- \n"
             "# Total%11"INT64_FMT" %11"INT64_FMT" %11"INT64_FMT" \n"
             "# ---------------------------------------- \n\n",
             totalElements, totalStations, totalBottom );

    hu_fclosep( &fp );

    /* output aggregate information to the monitor file / stdout */
    fprintf( stdout,
             "\nNonlinear mesh information\n"
             "Total number of nonlinear elements: %11"INT64_FMT"\n"
             "Total number of nonlinear stations: %11"INT64_FMT"\n"
             "Total number of bottom elements:    %11"INT64_FMT"\n\n",
             totalElements, totalStations, totalBottom );

    fflush( stdout );

}

/*
 * Collects statistics about number of nonlinear elements and stations in
 * nonlinear elements.
 */
void nonlinear_stats(int32_t myID, int32_t theGroupSize) {

    int32_t *nonlinElementsCount = NULL;
    int32_t *nonlinStationsCount = NULL;
    int32_t *bottomElementsCount = NULL;

    if ( myID == 0 ) {
        XMALLOC_VAR_N( nonlinElementsCount, int32_t, theGroupSize);
        XMALLOC_VAR_N( nonlinStationsCount, int32_t, theGroupSize);
        XMALLOC_VAR_N( bottomElementsCount, int32_t, theGroupSize);
    }

    MPI_Gather( &myNonlinElementsCount,    1, MPI_INT,
                nonlinElementsCount,       1, MPI_INT, 0, comm_solver );
    MPI_Gather( &myNumberOfNonlinStations, 1, MPI_INT,
                nonlinStationsCount,       1, MPI_INT, 0, comm_solver );
    MPI_Gather( &myBottomElementsCount,    1, MPI_INT,
                bottomElementsCount,       1, MPI_INT, 0, comm_solver );

    if ( myID == 0 ) {

        nonlinear_print_stats( nonlinElementsCount, nonlinStationsCount,
                               bottomElementsCount, theGroupSize);

        xfree_int32_t( &nonlinElementsCount );
    }

    return;
}

/*
 * nonlinear_solver_init: Initialize all the structures needed for nonlinear
 *                        analysis and the material/element constants.
 */
void nonlinear_solver_init(int32_t myID, mesh_t *myMesh, double depth) {

    int32_t eindex, nl_eindex;

    nonlinear_elements_count(myID, myMesh);
    nonlinear_elements_mapping(myID, myMesh);

    if ( theGeostaticLoadingT > 0 ) {
        bottom_elements_count(myID, myMesh, depth);
        bottom_elements_mapping(myID, myMesh, depth);
    }

    /* Memory allocation for mother structure */

    myNonlinSolver = (nlsolver_t *)malloc(sizeof(nlsolver_t));

    if (myNonlinSolver == NULL) {
        fprintf(stderr, "Thread %d: nonlinear_init: out of memory\n", myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }

    /* Memory allocation for internal structures */

    myNonlinSolver->constants =
        (nlconstants_t *)calloc(myNonlinElementsCount, sizeof(nlconstants_t));
    myNonlinSolver->stresses =
        (qptensors_t *)calloc(myNonlinElementsCount, sizeof(qptensors_t));
    myNonlinSolver->strains =
        (qptensors_t *)calloc(myNonlinElementsCount, sizeof(qptensors_t));
    myNonlinSolver->strains1 =
        (qptensors_t *)calloc(myNonlinElementsCount, sizeof(qptensors_t));
    myNonlinSolver->pstrains1 =
        (qptensors_t *)calloc(myNonlinElementsCount, sizeof(qptensors_t));
    myNonlinSolver->pstrains2 =
        (qptensors_t *)calloc(myNonlinElementsCount, sizeof(qptensors_t));
    myNonlinSolver->alphastress1 =
        (qptensors_t *)calloc(myNonlinElementsCount, sizeof(qptensors_t));
    myNonlinSolver->alphastress2 =
        (qptensors_t *)calloc(myNonlinElementsCount, sizeof(qptensors_t));
    myNonlinSolver->ep1 =
        (qpvectors_t *)calloc(myNonlinElementsCount, sizeof(qpvectors_t));
    myNonlinSolver->ep2 =
        (qpvectors_t *)calloc(myNonlinElementsCount, sizeof(qpvectors_t));
    myNonlinSolver->LoUnlo_n =
        (qpvectors_t *)calloc(myNonlinElementsCount, sizeof(qpvectors_t));
    myNonlinSolver->Sv_max =
        (qpvectors_t *)calloc(myNonlinElementsCount, sizeof(qpvectors_t));
    myNonlinSolver->Sv_n =
        (qpvectors_t *)calloc(myNonlinElementsCount, sizeof(qpvectors_t));
    myNonlinSolver->psi_n =
        (qpvectors_t *)calloc(myNonlinElementsCount, sizeof(qpvectors_t));
    myNonlinSolver->kappa =
        (qpvectors_t *)calloc(myNonlinElementsCount, sizeof(qpvectors_t));
    myNonlinSolver->Sref =
        (qptensors_t *)calloc(myNonlinElementsCount, sizeof(qptensors_t));


    if ( (myNonlinSolver->constants           == NULL) ||
         (myNonlinSolver->stresses            == NULL) ||
         (myNonlinSolver->strains             == NULL) ||
         (myNonlinSolver->strains1            == NULL) ||
         (myNonlinSolver->ep1                 == NULL) ||
         (myNonlinSolver->ep2                 == NULL) ||
         (myNonlinSolver->alphastress1        == NULL) ||
         (myNonlinSolver->alphastress2        == NULL) ||
         (myNonlinSolver->pstrains1           == NULL) ||
         (myNonlinSolver->pstrains2           == NULL) ||
         (myNonlinSolver->LoUnlo_n            == NULL) ||
         (myNonlinSolver->Sv_max              == NULL) ||
         (myNonlinSolver->Sv_n                == NULL) ||
         (myNonlinSolver->psi_n               == NULL) ||
         (myNonlinSolver->kappa               == NULL) ||
         (myNonlinSolver->Sref                == NULL) ) {

        fprintf(stderr, "Thread %d: nonlinear_init: out of memory\n", myID);
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }

    /* Initialization of element constants
     * Tensors have been initialized to 0 by calloc
     */

    for (nl_eindex = 0; nl_eindex < myNonlinElementsCount; nl_eindex++) {

        elem_t     *elemp;
        edata_t    *edata;
        nlconstants_t *ecp;
        double      mu, lambda;
        double      elementVs, elementVp;

        eindex = myNonlinElementsMapping[nl_eindex];

        elemp = &myMesh->elemTable[eindex];
        edata = (edata_t *)elemp->data;
        ecp   = myNonlinSolver->constants + nl_eindex;

        int32_t lnid0 = elemp->lnid[0];
        double  zo    = myMesh->ticksize * myMesh->nodeTable[lnid0].z;

        /* get element Vs */

        elementVs   = (double)edata->Vs;
        elementVp   = (double)edata->Vp;

        /* Calculate the lame constants and store in element */

        mu_and_lambda(&mu, &lambda, edata, eindex);
        ecp->lambda = lambda;
        ecp->mu     = mu;

        /* Calculate the vertical stress as a homogeneous half-space */
        if ( theApproxGeoState == YES )
        	ecp->sigmaZ_st = edata->rho * 9.80 * ( zo + edata->edgesize / 2.0 );

        /* Calculate the yield function constants */
        switch (theMaterialModel) {

            case LINEAR:
                ecp->alpha    = 0.0;
                ecp->phi      = 0.0;
                ecp->beta     = 0.0;
                ecp->h        = 0.0;
                ecp->Sstrain0 = 0.0;
                break;

            case VONMISES_FA:
            	ecp->c         = get_cohesion(elementVs);
            	ecp->phi       = 0.0;
            	ecp->dil_angle = 0.0;

            	ecp->alpha     = 0.0;
            	ecp->beta      = 0.0;
            	ecp->gamma     = 0.0;

            	ecp->Sstrain0  = interpolate_property_value(elementVs, theGamma0);

            	ecp->h         = 0.0;
            	ecp->psi0      = interpolate_property_value(elementVs, thePsi);
            	break;

            case VONMISES_FAM:
            	ecp->c         = get_cohesion(elementVs);
            	ecp->phi       = 0.0;
            	ecp->dil_angle = 0.0;

            	ecp->alpha     = 0.0;
            	ecp->beta      = 0.0;
            	ecp->gamma     = 0.0;

            	ecp->Sstrain0  = 0.0;

            	ecp->h         = 0.0;
            	ecp->psi0      = interpolate_property_value(elementVs, thePsi);
            	break;

            case VONMISES_BAE:
            	ecp->c         = get_cohesion(elementVs);
            	ecp->phi       = 0.0;
            	ecp->dil_angle = 0.0;

            	ecp->alpha     = 0.0;
            	ecp->beta      = 0.0;
            	ecp->gamma     = 0.0;

            	ecp->Sstrain0  = 0.0;
            	ecp->m         = interpolate_property_value(elementVs, theM);

            	ecp->h         = 0.0;
            	ecp->psi0      = interpolate_property_value(elementVs, thePsi);
            	break;

            case VONMISES_BAH:
            	ecp->c         = get_cohesion(elementVs);
            	ecp->phi       = 0.0;
            	ecp->dil_angle = 0.0;

            	ecp->alpha     = 0.0;
            	ecp->beta      = 0.0;
            	ecp->gamma     = 0.0;

            	ecp->Sstrain0  = 0.0;
            	ecp->m         = 0.0;

            	ecp->h         = 0.0;
            	ecp->psi0      = 0.0;
            	break;

            case VONMISES_GQH:
            	ecp->c         = get_cohesion(elementVs);
            	ecp->phi       = 0.0;
            	ecp->dil_angle = 0.0;

            	ecp->alpha     = 0.0;
            	ecp->beta      = 0.0;
            	ecp->gamma     = 0.0;

            	ecp->Sstrain0  = 0.0;
            	ecp->m         = 0.0;

            	ecp->h         = 0.0;
            	ecp->psi0      = 0.0;

            	ecp->thetaGQH[0] = interpolate_property_value(elementVs, theTheta1);
            	ecp->thetaGQH[1] = interpolate_property_value(elementVs, theTheta2);
            	ecp->thetaGQH[2] = interpolate_property_value(elementVs, theTheta3);
            	ecp->thetaGQH[3] = interpolate_property_value(elementVs, theTheta4);
            	ecp->thetaGQH[4] = interpolate_property_value(elementVs, theTheta5);

            	break;

            case VONMISES_MKZ:
            	ecp->c           = get_cohesion(elementVs);

            	ecp->beta_MKZ    = interpolate_property_value(elementVs, thePsi);
            	ecp->s_MKZ       = interpolate_property_value(elementVs, theM);

            	break;

            case VONMISES_RO:
            	ecp->c           = get_cohesion(elementVs);

            	ecp->alpha_RO    = interpolate_property_value(elementVs, thePsi);
            	ecp->eta_RO      = interpolate_property_value(elementVs, theM);
            	ecp->tauy_RO     = interpolate_property_value(elementVs, theTau_y);

            	break;


            case VONMISES_EP:
            	ecp->c         = get_cohesion(elementVs);
            	ecp->phi       = 0.0;
            	ecp->dil_angle = 0.0;


            	ecp->alpha     = 0.0;
            	ecp->beta      = 0.0;
            	ecp->gamma     = 1.0;

            	ecp->Sstrain0  = 0.0;

            	ecp->h         = 0; /*  no isotropic hardening  in von Mises model */
            	break;

            case DRUCKERPRAGER:
                ecp->c         = get_cohesion(elementVs);
                ecp->phi       = get_phi(elementVs);
                ecp->dil_angle = get_dilatancy(elementVs);
                ecp->h         = get_hardmod(elementVs);

                ecp->alpha     =  2.0 * sin(ecp->phi)       / ( sqrt(3.0) * ( 3.0 - sin(ecp->phi) ) );
                ecp->beta      =  2.0 * sin(ecp->dil_angle) / ( sqrt(3.0) * ( 3.0 - sin(ecp->dil_angle) ) );
                ecp->gamma     =  6.0 * cos(ecp->phi)       / ( sqrt(3.0) * ( 3.0 - sin(ecp->phi) ) );

                ecp->Sstrain0  = 0.0;

            	break;
            case MOHR_COULOMB:
                ecp->c     = get_cohesion(elementVs);
                ecp->phi   = get_phi(elementVs);
                ecp->dil_angle = get_dilatancy(elementVs);

                ecp->alpha = 0.;
                ecp->beta  = 0.;
                ecp->gamma = 0.;

                ecp->h         = get_hardmod(elementVs);
                ecp->Sstrain0  = 0.0;
            	break;

            default:
                fprintf(stderr, "Thread %d: nonlinear_solver_init:\n"
                        "\tUnexpected error with the material model\n", myID);
                MPI_Abort(MPI_COMM_WORLD, ERROR);
                exit(1);
                break;
        }


        ecp->strainrate  =
        		interpolate_property_value(elementVs, theStrainRates  );
        ecp->sensitivity =
        		interpolate_property_value(elementVs, theSensitivities );

    } /* for all elements */

}

/* -------------------------------------------------------------------------- */
/*                   Auxiliary tensor manipulation methods                    */
/* -------------------------------------------------------------------------- */

/*
 * Returns the isotropic tensor B = lambda*I .
 */
tensor_t isotropic_tensor(double lambda) {

    tensor_t B;

    B.xx = lambda;
    B.yy = lambda;
    B.zz = lambda;
    B.xy = 0.0;
    B.yz = 0.0;
    B.xz = 0.0;

    return B;
}


/*
 * Returns the scaled tensor B = lambda*A .
 */
tensor_t scaled_tensor(tensor_t A, double lambda) {

    tensor_t B;

    B.xx = lambda*A.xx;
    B.yy = lambda*A.yy;
    B.zz = lambda*A.zz;
    B.xy = lambda*A.xy;
    B.yz = lambda*A.yz;
    B.xz = lambda*A.xz;

    return B;
}

/*
 * tensor_I1: Returns the invariant I1 of a tensor.
 */
double tensor_I1(tensor_t tensor) {

    return tensor.xx + tensor.yy + tensor.zz;
}

/*
 * Returns the octahedral of a tensor given its first invariant I1
 */
double tensor_octahedral(double I1) {

    return I1/3.0;
}

/*
 * Returns the deviator of a tensor given the tensor and its octahedral.
 */
tensor_t tensor_deviator(tensor_t tensor, double oct) {

    tensor_t deviator;

    deviator.xx = tensor.xx - oct;
    deviator.yy = tensor.yy - oct;
    deviator.zz = tensor.zz - oct;
    deviator.xy = tensor.xy;
    deviator.yz = tensor.yz;
    deviator.xz = tensor.xz;

    return deviator;
}

/*
 * tensor_J2: Returns the second invariant of a tensor given its deviator.
 */
double tensor_J2(tensor_t dev) {

    return ( (dev.xx * dev.xx) + (dev.yy * dev.yy) + (dev.zz * dev.zz) ) * 0.5
           + (dev.xy * dev.xy) + (dev.yz * dev.yz) + (dev.xz * dev.xz);
}

/*
 * combtensor_J2: Returns the combined second invariant J2_comb = (A:B)/2
 * A y B symmetric tensors
 */
double combtensor_J2(tensor_t A, tensor_t B) {

    return ( (A.xx * B.xx) + (A.yy * B.yy) + (A.zz * B.zz) ) * 0.5
           + (A.xy * B.xy) + (A.yz * B.yz) + (A.xz * B.xz);
}

/*
 * tensos_Det: Returns the determinan of a tensor.
 */
double tensor_J3(tensor_t dev) {

    return ( dev.xx * ( dev.yy * dev.zz - dev.yz * dev.yz ) - dev.xy * ( dev.xy * dev.zz - dev.xz * dev.yz ) + dev.xz * ( dev.xy * dev.yz - dev.xz * dev.yy ) );
}

/*
 * Computes the contribution of the i-th node to the three derivatives
 * dx, dy, dz of the shape functions evaluated at local coordinates
 * lx, ly, lz.
 */
void point_dxi ( double *dx, double *dy, double *dz,
                 double  lx, double  ly, double  lz, double h, int i )
{
    double XI;
    double Jij = 0.25 / h; /* Jacobian 1/(4h) */

    *dx = Jij * (       xi[0][i]      )
              * ( 1.0 + xi[1][i] * ly )
              * ( 1.0 + xi[2][i] * lz );

    *dy = Jij * ( 1.0 + xi[0][i] * lx )
              * (       xi[1][i]      )
              * ( 1.0 + xi[2][i] * lz );

    *dz = Jij * ( 1.0 + xi[0][i] * lx )
              * ( 1.0 + xi[1][i] * ly )
              * (       xi[2][i]      );

    return;
}

/*
 * Computes the three derivatives of the shape functions evaluated at
 * coordinates j-th quadrature point.
 */
void compute_qp_dxi (double *dx, double *dy, double *dz, int i, int j, double h)
{
    double XI, QC;

    /* quadrature point local coordinates */
    double lx = xi[0][j] * qc ;
    double ly = xi[1][j] * qc ;
    double lz = xi[2][j] * qc ;

    point_dxi(dx, dy, dz, lx, ly, lz, h, i);

    return;
}

/*
 * Resets a tensor to zero in all its components.
 */
tensor_t init_tensor() {

    tensor_t tensor;

    tensor.xx = 0.0;
    tensor.yy = 0.0;
    tensor.zz = 0.0;
    tensor.xy = 0.0;
    tensor.yz = 0.0;
    tensor.xz = 0.0;

    return tensor;
}

void init_tensorptr(tensor_t *tensor) {

    tensor->xx = 0.0;
    tensor->yy = 0.0;
    tensor->zz = 0.0;
    tensor->xy = 0.0;
    tensor->yz = 0.0;
    tensor->xz = 0.0;

    return;
}

/*
 * Compute strain tensor of a given point in the element.
 */
tensor_t point_strain (fvector_t *u, double lx, double ly, double lz, double h) {

    int i;

    tensor_t strain = init_tensor();

    /* Contribution of each node */
    for (i = 0; i < 8; i++) {

        double dx, dy, dz;

        point_dxi(&dx, &dy, &dz, lx, ly, lz, h, i);

        strain.xx += dx * u[i].f[0];
        strain.yy += dy * u[i].f[1];
        strain.zz += dz * u[i].f[2];

        strain.xy += 0.5 * ( dy * u[i].f[0] + dx * u[i].f[1] );
        strain.yz += 0.5 * ( dz * u[i].f[1] + dy * u[i].f[2] );
        strain.xz += 0.5 * ( dz * u[i].f[0] + dx * u[i].f[2] );

    } /* nodes contribution */

    return strain;
}

/*
 * Computes the stress tensor in a given point in the element from the point's
 * strain tensor and the element properties according to the linear elastic
 * stress-strain relationship.
 */
tensor_t point_stress (tensor_t strain, double mu, double lambda) {

    double mu2, lkk;
    double strain_kk;
    tensor_t stress;

    /* calculate strtain_kk */
    strain_kk = tensor_I1(strain);

    mu2 = 2.0 * mu;
    lkk = lambda * strain_kk;

    stress.xx = mu2 * strain.xx + lkk;
    stress.yy = mu2 * strain.yy + lkk;
    stress.zz = mu2 * strain.zz + lkk;

    stress.xy = mu2 * strain.xy;
    stress.yz = mu2 * strain.yz;
    stress.xz = mu2 * strain.xz;

    return stress;
}

/*
 * Computes the strain tensor given the stress tensor and the material properties.
 */
tensor_t elastic_strains (tensor_t stress, double mu, double kappa) {

    double Skk, e_kk;
    tensor_t strain, strain_dev;

    Skk = tensor_I1(stress);
    e_kk = Skk / ( 9.0 * kappa);

    strain_dev =  tensor_deviator( stress, Skk / 3.0 );

    strain.xx = 1./(2. * mu) * strain_dev.xx + e_kk;
    strain.yy = 1./(2. * mu) * strain_dev.yy + e_kk;
    strain.zz = 1./(2. * mu) * strain_dev.zz + e_kk;
    strain.xy = 1./(2. * mu) * strain_dev.xy;
    strain.xz = 1./(2. * mu) * strain_dev.xz;
    strain.yz = 1./(2. * mu) * strain_dev.yz;

    return strain;
}

/*
 * Computes the stress tensors for the quadrature points of an element given
 * the strain tensors.
 */
qptensors_t compute_qp_stresses (qptensors_t strains, double mu, double lambda)
{
    int i;
    qptensors_t stresses;

    /* Loop over the quadrature points */
    for ( i = 0; i < 8; i++ ) {
        stresses.qp[i] = point_stress(strains.qp[i], mu, lambda);
    }

    return stresses;
}

/*
 * Returns the subtraction between two tensors.
 */
tensor_t subtrac_tensors(tensor_t A, tensor_t B) {

    tensor_t C;

    C.xx = A.xx - B.xx;
    C.yy = A.yy - B.yy;
    C.zz = A.zz - B.zz;
    C.xy = A.xy - B.xy;
    C.yz = A.yz - B.yz;
    C.xz = A.xz - B.xz;

    return C;
}


/*
 * Returns the double product between two symetric tensors.
 */
double ddot_tensors(tensor_t A, tensor_t B) {

    double C;

    C = A.xx * B.xx + A.yy * B.yy + A.zz * B.zz + 2 * (A.xy * B.xy + A.yz * B.yz + A.xz * B.xz);

    return C;
}

/*
 * Returns the summation of two tensors.
 */
tensor_t add_tensors(tensor_t A, tensor_t B) {

    tensor_t C;

    C.xx = A.xx + B.xx;
    C.yy = A.yy + B.yy;
    C.zz = A.zz + B.zz;
    C.xy = A.xy + B.xy;
    C.yz = A.yz + B.yz;
    C.xz = A.xz + B.xz;

    return C;
}

/*
 * Returns the ZERO tensor.
 */
tensor_t zero_tensor() {

    tensor_t C;

    C.xx = 0.0;
    C.yy = 0.0;
    C.zz = 0.0;
    C.xy = 0.0;
    C.yz = 0.0;
    C.xz = 0.0;

    return C;
}

/*
 * Returns the approximated self-weight tensor.
 */
tensor_t ApproxGravity_tensor(double Szz, double phi, double h, double lz, double rho) {

	double Ko = 1 - sin(phi);
	double Sigma = -( Szz + 9.8 * rho * h * lz * 0.5);

    tensor_t C;

    C.xx = Ko * Sigma;
    C.yy = Ko * Sigma;
    C.zz = Sigma;
    C.xy = 0.0;
    C.yz = 0.0;
    C.xz = 0.0;

    return C;
}

/*
 * Returns the subtraction between two tensors for all quadrature points.
 */
qptensors_t subtrac_qptensors(qptensors_t A, qptensors_t B) {

    int i;
    qptensors_t C;

    /* Loop over quadrature points */
    for (i = 0; i < 8; i++) {
        C.qp[i] = subtrac_tensors(A.qp[i], B.qp[i]);
    }

    return C;
}

tensor_t copy_tensor (tensor_t original) {

    tensor_t copy;

    copy.xx = original.xx;
    copy.yy = original.yy;
    copy.zz = original.zz;
    copy.xy = original.xy;
    copy.yz = original.yz;
    copy.xz = original.xz;

    return copy;
}

double compute_yield_surface_stateII ( double J3, double J2, double I1, double alpha, double phi, tensor_t Sigma ) {

	double Yf=0., p, q, r, teta, Rmc, s1, s3;

	if ( ( theMaterialModel == VONMISES_EP )  || ( theMaterialModel == DRUCKERPRAGER ) ||
	     ( theMaterialModel == VONMISES_FAM ) || ( theMaterialModel == VONMISES_FA   ) ||
	     ( theMaterialModel == VONMISES_BAE ) || ( theMaterialModel == VONMISES_BAH  ) ||
	     ( theMaterialModel == VONMISES_GQH   )  ) {
		if ( theMaterialModel == DRUCKERPRAGER )
			Yf = alpha * I1 + sqrt( J2 );
		else
			Yf = sqrt( J2 ); // alpha must be zero in any vonMises model
	} else {
		p = (1./3.) * I1;
		q = 2.0 * pow( J2, 1.5 );
		r = -3.0 * sqrt(3.0) * (J3);

		if ( ( r/q <= 1.00000001 ) && ( r/q >= 0.99999999 ) )
			teta = PI / 6.0;
		else if ( ( r/q >= -1.00000001 ) && ( r/q <= -0.99999999 ) )
			teta = -PI / 6.0;
		else
			teta = 1./3. * asin(r/q);

		if ( ( teta > PI / 6.0 ) || ( teta < -PI / 6.0 ) ) {

			vect1_t n1, n2, n3, eig_vals;

		    specDecomp(Sigma, &n1, &n2, &n3, &eig_vals);

		    s1 = eig_vals.x;
		    s3 = eig_vals.z;

		    Yf = s1 - s3 + ( s1 + s3 )*sin(phi);

		} else {
			Rmc = -1./sqrt(3.0) * ( sin(phi)*sin(teta) ) + cos(teta);
			Yf = 2. * ( Rmc * sqrt(J2) + p * sin(phi) );
		}

	}

	return Yf;

}

double compute_hardening ( double gamma, double c, double Sy, double h, double ep_bar, double phi, double psi ) {
	double H=0.;

	if ( theMaterialModel == VONMISES_EP ) {
		H = c * 2.0 / sqrt(3.0);   // c=Su and tao_max = 2Su/sqrt(3)
	} else if ( theMaterialModel == VONMISES_FA ) {
		H = Sy; // Since Sy comes from the G/Gmax it does not require the constant 2/sqrt(3)
	} else if ( theMaterialModel == VONMISES_FAM || theMaterialModel == VONMISES_BAE ||
			    theMaterialModel == VONMISES_BAH || theMaterialModel == VONMISES_GQH ) { // no elastic region in vonMisesKinHard_Modified
		H = 0.0;
	} else if ( theMaterialModel == DRUCKERPRAGER ) {
		H = gamma * ( c + h * ep_bar);
	}	else {
		H = 2.0 * ( c + h * ep_bar) * cos(phi);
	}


	return H;

}

/*===============================================================*/
/*===============================================================*/
/*   Material update function for material models based on (1994) Borja & Amies approach    */
void MatUpd_vMGeneral ( nlconstants_t el_cnt, double *kappa,
		                tensor_t e_n, tensor_t e_n1, tensor_t *sigma_ref, tensor_t *sigma,
		                int *FlagTolSubSteps, int *FlagNoSubSteps, double *ErrMax ) {

/*	 INPUTS:
    * el_cnt            : Material constants

 	* e_n           	: Total strain tensor.
 	* e_n1         		: Total strain tensor at t-1
    * sigma_ref     	: reference stress
    * substepTol, BoundSurfTol 	: substep Tolerance, Bounding surface Tolerance

 	* OUTPUTS:
 	* sigma         : Updated stress tensor
 	* kappa         : Updated hardening variable
    * Sref          : Updated reference deviator stress tensor          */

	double   Dt=1.0, T=0.0, Dtmin, Dt_sup, kappa_n, load_unload, Den1, Den2, kappa_up,
			 ErrB, ErrS, xi, xi_sup, kappa_o, K,
			 G=el_cnt.mu, Lambda = el_cnt.lambda, popo;
	tensor_t sigma_n, sigma_up, Num, Sdev;
	int cnt=0;

	Dtmin = Dt/theNoSubsteps;
	K     = Lambda + 2.0 * G / 3.0;

	/* At  this point *sigma and *kappa have the information at t-1 */
	kappa_n = *kappa;
	sigma_n = copy_tensor(*sigma);

	/* deviatoric stress at t-1. At  this point *sigma has the information at t-1  */
	tensor_t Sdev_n1   = tensor_deviator( *sigma, tensor_octahedral ( tensor_I1 ( *sigma ) ) );


	/* total strain increment and deviatoric strain increment */
	tensor_t De       = subtrac_tensors ( e_n, e_n1 );
	double   De_vol   = tensor_I1 ( De );
	tensor_t De_dev   = tensor_deviator( De, tensor_octahedral ( De_vol ) );

	Den1 = ddot_tensors(Sdev_n1, subtrac_tensors (Sdev_n1 , *sigma_ref));
	Den2 = kappa_n * ( ddot_tensors(subtrac_tensors (Sdev_n1 , *sigma_ref), subtrac_tensors (Sdev_n1 , *sigma_ref)) );
	Num  = add_tensors ( scaled_tensor( Sdev_n1, (1.0+kappa_n) ), scaled_tensor( (subtrac_tensors (Sdev_n1 , *sigma_ref) ) ,kappa_n*(1.0+kappa_n) ) );

	load_unload = -ddot_tensors(Num,De_dev) / (Den1 + Den2);

	if ( load_unload > 0 ) {

		*kappa = get_kappaUnLoading_II(  el_cnt, Sdev_n1,  De_dev, ErrMax );
	    *sigma_ref = copy_tensor( Sdev_n1 );


		/* get sigma_n deviatoric */
		double H_n      = getHardening( el_cnt, *kappa);
		double xi1      = 2.0 * G / ( 1.0 + 3.0 * G / H_n );
		tensor_t DSdev  = scaled_tensor( De_dev, xi1 );

		if ( isnan( tensor_J2(DSdev) ) || isinf( tensor_J2(DSdev) ) ) {
			fprintf(stdout," NAN at unloading: %f.  \n",tensor_J2(DSdev));
	        MPI_Abort(MPI_COMM_WORLD, ERROR);
	        exit(1);
		}

		*sigma          = add_tensors (  add_tensors( sigma_n, isotropic_tensor(K * De_vol) ),  DSdev  );
		return;

	}

	EvalSubStep ( el_cnt,  sigma_n,  De,  De_dev,  De_vol, Dt,  sigma_ref,  &sigma_up,  kappa_n, &kappa_up,  &ErrB,  &ErrS);


	double Emax     = 0;
	int    step_Emax = -1, i;

	if ( ErrB > theErrorTol ) { // begin sub-stepping
		xi_sup  = 0.0;
		kappa_o = kappa_n;

	    for (i = 0; i < theNoSubsteps ; i++) {

	    	while ( ErrB > theErrorTol ) {

	    		Dt_sup = MIN( xi_sup * Dt, 1-T );
	    		xi     = MAX( 0.9 * sqrt(theErrorTol / ErrB), 0.10 );
	    		Dt     = MAX( xi * Dt, Dtmin );
	    		Dt     = MIN( Dt, 1 - T );

	    		/*  compute state for xi_sup (xi_sup is an extrapolated value of xi)  */
	    		if ( Dt_sup > Dt ) {
	    			EvalSubStep (  el_cnt, sigma_n, De, De_dev,  De_vol, Dt_sup,  sigma_ref,  &sigma_up, kappa_n, &kappa_up, &ErrB, &ErrS );
	    		}

	    		if ( ErrB > theErrorTol ) {
	    			EvalSubStep ( el_cnt, sigma_n, De, De_dev, De_vol, Dt, sigma_ref, &sigma_up, kappa_n, &kappa_up,  &ErrB, &ErrS);
	    			xi_sup = 0.0;  // forget previous xi_sup
	    		} else
	    			Dt = Dt_sup;

	    		if ( Dt == Dtmin )
	    			break;

	    		cnt = cnt + 1;

	    		if (cnt > theNoSubsteps)
	    			break;
	    	}

	        if ( (Dt == Dtmin) && (ErrB > theErrorTol) ) {
	            if (ErrB > Emax) {
	            	*ErrMax = ErrB;
	                step_Emax = i;
	            }
	        }

	        /* Update initial values  */
	        sigma_n = copy_tensor(sigma_up);
	        kappa_n = kappa_up;
	        Sdev    = tensor_deviator( sigma_n, tensor_octahedral ( tensor_I1 ( sigma_n ) ) );
	        *kappa  = get_kappa(  el_cnt, Sdev,  *sigma_ref,  kappa_o  );
	        *sigma  = copy_tensor(sigma_up);
	        T       = T + Dt;

	        if ( T == 1 ) {
	        	if ( step_Emax != -1 ) {
	        		*FlagTolSubSteps = 1;
	        	}

	    		if ( isnan( tensor_J2(Sdev) ) || isinf( tensor_J2(Sdev) ) ) {
	    			fprintf(stdout," NAN or INF at T=1. J2= %f.  \n",tensor_J2(Sdev));
	    	        MPI_Abort(MPI_COMM_WORLD, ERROR);
	    	        exit(1);
	    		}
	    		*ErrMax = ErrB;
	        	return;
	        }
	        xi_sup = MIN(0.9*sqrt(theErrorTol/ErrB),1.1);
	        ErrB   = 1E10;
	    }

	    if ( i == theNoSubsteps - 1 ) { // reached maximum sub-steps
	    		*FlagNoSubSteps = 1.0;

	    		if ( isnan( tensor_J2(Sdev) ) || isinf( tensor_J2(Sdev) ) ) {
	    			fprintf(stdout," NAN or INF when reaching Maxsubsteps. J2=%f.  \n",tensor_J2(Sdev));
	    	        MPI_Abort(MPI_COMM_WORLD, ERROR);
	    	        exit(1);
	    		}
	    		*ErrMax = ErrB;
	    		return;
	    }
	} else {
		*kappa = kappa_up;
		*sigma  = copy_tensor(sigma_up);
		*ErrMax = ErrB;

		if ( isnan(tensor_J2(*sigma)) || isinf(tensor_J2(*sigma)) ) {
			fprintf(stdout," NAN without substepping. J2=%f, kappa=%f,  Sxx=%f, Syy=%f, Szz=%f."
					"\n",tensor_J2(*sigma), kappa_up,sigma->xx,sigma->yy, sigma->zz);
	        MPI_Abort(MPI_COMM_WORLD, ERROR);
	        exit(1);
		}
	}

}


void EvalSubStep (nlconstants_t el_cnt, tensor_t  sigma_n, tensor_t De, tensor_t De_dev, double De_vol,
		          double Dt, tensor_t *sigma_ref, tensor_t *sigma_up, double kappa_n,
		          double *kappa_up, double *ErrB, double *ErrS) {

	tensor_t Sdev_0, DSdev1, DSdev2, Sdev1, Sdev2, Dsigma1, Dsigma2, Dss;
	double   H_n, H_n2, xi1, xi2, K, kappa1, kappa2, Su=el_cnt.c, Lambda=el_cnt.lambda, G=el_cnt.mu;

	K       = Lambda + 2.0 * G / 3.0;
	De 		= scaled_tensor(De,Dt);
	De_dev 	= scaled_tensor(De_dev,Dt);
	De_vol 	= De_vol * Dt;

	/* get sigma_n deviatoric */
	Sdev_0   = tensor_deviator( sigma_n, tensor_octahedral ( tensor_I1 ( sigma_n ) ) );
	H_n      = getHardening( el_cnt, kappa_n );
	xi1      = 2.0 * G / ( 1.0 + 3.0 * G / H_n );

	DSdev1   = scaled_tensor( De_dev,xi1 );
	Sdev1    = add_tensors( Sdev_0, DSdev1 );
	Dsigma1  = add_tensors( isotropic_tensor(K * De_vol), Sdev1 );
	kappa1   = get_kappa( el_cnt, Sdev1, *sigma_ref, kappa_n );

	/* get second set of values */
	H_n2     = getHardening( el_cnt, kappa1 );
	xi2      = 2.0 * G / ( 1.0 + 3.0 * G / H_n2 );

	DSdev2   = scaled_tensor( De_dev,xi2 );
	Sdev2    = add_tensors( Sdev_0, DSdev2 );
	kappa2   = get_kappa( el_cnt, Sdev2, *sigma_ref, kappa_n );
	Dsigma2  = add_tensors( isotropic_tensor(K * De_vol), Sdev2 );

	*sigma_up = add_tensors (  add_tensors( sigma_n, isotropic_tensor(K*De_vol) ),
			                   add_tensors( scaled_tensor(DSdev1,0.5), scaled_tensor(DSdev2,0.5) ) );
	*kappa_up = (kappa1+kappa2)/2;

	/* compute errors */
	Dss       =  subtrac_tensors(Dsigma2,Dsigma1);
	*ErrS     =  sqrt(2.0 *  tensor_J2 ( Dss ) ) / sqrt(2.0 *  tensor_J2 ( *sigma_up ) );

	double R  = Su * sqrt(8.0/3.0);

	tensor_t SmSo = subtrac_tensors( add_tensors( scaled_tensor(Sdev1,0.5), scaled_tensor(Sdev2,0.5) ), *sigma_ref );
	tensor_t S1   = add_tensors    ( add_tensors( scaled_tensor(Sdev1,0.5), scaled_tensor(Sdev2,0.5) ), scaled_tensor(SmSo,*kappa_up) );

	*ErrB         = fabs( sqrt( ddot_tensors(S1,S1) ) - R ) / R;

}


double getH_GQHmodel (nlconstants_t el_cnt, double kappa) {

double H = 0.0, tao_bar, Theta1, Theta2, Theta3, Theta4, Theta5,
	   A1, B1, C1, gamma_baro, gamma_bar, theta, Dergamma, Eo, Jac, A, B, C;

int    cnt=0, cnt_max=200;

if (kappa == 0.0)
		return H;

Theta1 = el_cnt.thetaGQH[0];
Theta2 = el_cnt.thetaGQH[1];
Theta3 = el_cnt.thetaGQH[2];
Theta4 = el_cnt.thetaGQH[3];
Theta5 = el_cnt.thetaGQH[4];

tao_bar = 1.0/(1.0 + kappa);

// approximate initial gamma_bar assuming Theta4=Theta5=1
A1 = 1.0 - tao_bar;
B1 = Theta3 * A1 - tao_bar + ( Theta1 + Theta2 ) * tao_bar * tao_bar;
C1 = tao_bar * Theta3 * ( Theta1 * tao_bar - 1.0 );
gamma_baro = (-B1 + sqrt( B1*B1 -4.0 * A1 * C1 ) ) / ( 2.0 * A1 );

// get theta initial
if (gamma_baro == 0) {
	theta    = Theta1;
	Dergamma = 0.0;
	H = FLT_MAX;
	return H;
} else {
	theta    = Theta1 + Theta2 * ( Theta4 * pow(gamma_baro,Theta5) ) / ( pow(Theta3,Theta5) + Theta4 * pow(gamma_baro,Theta5) );
	Dergamma = Theta2 * ( pow(Theta3,Theta5) ) * Theta4 * Theta5 * pow(gamma_baro,(Theta5 - 1.0)) / (pow(( pow(Theta3,Theta5) + Theta4 * pow(gamma_baro,Theta5) ),2.0 ));
}


Eo  = gamma_baro * ( 1.0 - tao_bar ) - tao_bar + tao_bar * tao_bar * theta;
gamma_bar = gamma_baro;

while ( fabs(Eo) > 1E-10 ) {
	Jac       = (1.0 - tao_bar) + tao_bar * tao_bar * Dergamma;
	gamma_bar = gamma_bar - Eo/Jac;

	if (gamma_bar == 0) {
		theta    = Theta1;
		Dergamma = 0.0;
	} else {
		theta     = Theta1 + Theta2 * ( Theta4 * pow(gamma_bar,Theta5) ) / ( pow(Theta3,Theta5) + Theta4 * pow(gamma_bar,Theta5) );
		Dergamma  = Theta2 * ( pow(Theta3,Theta5) ) * Theta4 * Theta5 * pow(gamma_bar,(Theta5 - 1.0)) / (pow(( pow(Theta3,Theta5) + Theta4 * pow(gamma_bar,Theta5) ), 2.0 ));
	}

	Eo        = gamma_bar * ( 1.0 - tao_bar ) - tao_bar + tao_bar * tao_bar * theta;
	cnt++;
	if (cnt == cnt_max)
		break;
}

/*  Sanity check.   */
if ( cnt == cnt_max ) {
    fprintf(stderr,"Material update warning: "
            "could not find GammaBar for GQ/H model. GammaBar=%f, Error=%f\n", gamma_bar, Eo);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
}

A =  1.0 + gamma_bar + sqrt( (1.0 + gamma_bar) * (1.0 + gamma_bar) - 4.0 * theta * gamma_bar );
C =  1.0 + ( 1+ gamma_bar - 2.0 * (theta + gamma_bar * Dergamma) )/sqrt( (1.0 + gamma_bar) * (1.0 + gamma_bar) - 4.0 * theta * gamma_bar );
B =  A - gamma_bar * C;

H = el_cnt.mu * 4.0 * B / ( A * A - 2.0 * B ) * 3.0 / 2.0; // the 3.0/2.0 scaling factor converts Kp to H as needed in (1994) Borja & Amies technique

return H ;

}


/*  Plastic modulus */
double getHardening(nlconstants_t el_cnt, double kappa) {

	double H = 0, psi=el_cnt.psi0, m=el_cnt.m, G=el_cnt.mu;

	if ( kappa == FLT_MAX ) {
		H = FLT_MAX;
		return H;
	}
	if ( theMaterialModel == VONMISES_BAE )
		H = ( psi * G ) * pow( kappa, m );
	else {
		if ( theMaterialModel == VONMISES_BAH )
			H = 3.0 * G * pow(kappa,2.0) / ( 1.0 + 2.0 * kappa );
		else {
			if ( theMaterialModel == VONMISES_GQH )
				H = getH_GQHmodel ( el_cnt,  kappa);
		}
	}

	return H;
}

/*

 Derivative of the plastic modulus
double getDerHardening(nlconstants_t el_cnt, double kappa) {

	double H = 0, psi=el_cnt.psi0, m=el_cnt.m, G=el_cnt.mu;

	if ( kappa == FLT_MAX ) {
		if ( theMaterialModel == VONMISES_BAE )
			H = FLT_MAX;
		else
			H = 1.50 * G;
		return H;
	}

	if ( theMaterialModel == VONMISES_BAE )
		H = ( psi * G * m ) * pow( kappa, m - 1.0 );
	else {
		if ( theMaterialModel == VONMISES_BAH )
			H = 6.0 * kappa * G * ( 1.0 + kappa ) / ( ( 1.0 + 2.0 * kappa ) * ( 1.0 + 2.0 * kappa ) );
		else{
			if ( theMaterialModel == VONMISES_GQH )
				H = ( getH_GQHmodel ( el_cnt,  kappa * 1.00001) - getH_GQHmodel ( el_cnt,  kappa) ) / 0.00001; // todo: get the exact expression for the derivative of the hardening function;
		}
	}

	return H;
}

*/


double get_kappa( nlconstants_t el_cnt, tensor_t Sdev, tensor_t Sref, double kn ) {

	double R, Fk, Dk, Jk, kappa;
	int    cnt=0, cnt_max=200;
	tensor_t SmSo, S1;

	double Su=el_cnt.c;

	kappa = kn;
	R = sqrt(8.0/3.0) * Su;

	SmSo = subtrac_tensors(Sdev,Sref);
	S1   = add_tensors(Sdev, scaled_tensor(SmSo,kappa));

	Fk   = sqrt(ddot_tensors(S1,S1)) - R;

	while ( fabs(Fk) > theErrorTol ) {
		Jk     = ddot_tensors(SmSo,S1)/(sqrt(ddot_tensors(S1,S1)));
		Dk     = -Fk / Jk;
		kappa  = kappa + Dk;
		S1     = add_tensors(Sdev, scaled_tensor(SmSo,kappa));
		Fk     = sqrt(ddot_tensors(S1,S1)) - R;
		cnt = cnt + 1;
		if (cnt == cnt_max)
			break;
	}

	if ( kappa < 0)
		kappa = kn;

	return kappa;

}


double Pegasus(double beta, nlconstants_t el_cnt) {
	// (1973) King, R. An Improved Pegasus method for root finding

	double k0 = 0.0, k1 = 1.0, k2,  f0, f1, f2,  G=el_cnt.mu, tmp;
	int cnt1=1, cnt2=1, cntMax = 200;

	f0 = ( 1.0 + k0 - beta ) * getHardening(el_cnt, k0) / G - 3.0 * beta;
	f1 = ( 1.0 + k1 - beta ) * getHardening(el_cnt, k1) / G - 3.0 * beta;

	// get initial range for k
	while ( f0 * f1 > 0 && cnt1 < cntMax ) {
		k0 = k1;
		k1 = 2.0 * k1;

		f0 = ( 1.0 + k0 - beta ) * getHardening(el_cnt, k0) / G - 3.0 * beta;
		f1 = ( 1.0 + k1 - beta ) * getHardening(el_cnt, k1) / G - 3.0 * beta;
		cnt1++;
	}

	if (cnt1 == cntMax) {
		fprintf(stdout,"Cannot obtain initial range for kappa at unloading: k0=%f, k1=%f. \n", k0, k1 );
		MPI_Abort(MPI_COMM_WORLD, ERROR);
		exit(1);
	}

	cnt1=1;

	k2 = k1 - f1 * ( k1 - k0 ) / ( f1 - f0 );
	f2 = ( 1.0 + k2 - beta ) * getHardening(el_cnt, k2) / G - 3.0 * beta;

	while ( fabs(f2) > theErrorTol && cnt1 < cntMax ) {

		if ( f1 * f2 < 0 ) {
			tmp = k0;
			k0  = k1;
			k1  = tmp;
			tmp = f0;
			f0  = f1;
			f1  = tmp;
		}

		while( f1 * f2 > 0 && cnt2 < cntMax) {
			f0 = f0 * f1 / ( f1 + f2 );

			k1 = k2;
			f1 = f2;

			k2 = k1 - f1 * ( k1 - k0 ) / ( f1 - f0 );
			f2 = ( 1.0 + k2 - beta ) * getHardening(el_cnt, k2) / G - 3.0 * beta;

			cnt2++;
		}

		k0 = k1;
		f0 = f1;

		k1 = k2;
		f1 = f2;

		if ( k0 == k1 ) {
			k2 = k1;
			break;
		}

		k2 = k1 - f1 * ( k1 - k0 ) / ( f1 - f0 );
		f2 = ( 1.0 + k2 - beta ) * getHardening(el_cnt, k2) / G - 3.0 * beta;

		cnt1++;

	}

	if ( cnt1 == cntMax || cnt2 == cntMax )
		fprintf(stdout,"Increase the number of steps for root finding in Pegasus method. k=%f, beta=%f, Error=%f, ErroTol=%f \n", k2, beta, fabs(f2), theErrorTol );

	return k2;

}


double get_kappaUnLoading_II( nlconstants_t el_cnt, tensor_t Sn, tensor_t De, double *Err ) {

	double R, A, B, C, kappa1, kappa2, beta, phi, G=el_cnt.mu, Su=el_cnt.c, kn;

	R     = sqrt(8.0/3.0) * Su;

	A = ddot_tensors(Sn,Sn);
	B = 2.0 * ddot_tensors(Sn,De);
	C = ddot_tensors(De,De);

	/* ========== Improve kn value =========== */


	if ( ( B*B - 4.0*C*(A-R*R) ) > 0 ) {

		kappa1 = (-B + sqrt( B*B - 4.0*C*(A-R*R) ) ) / (2 * C);
		kappa2 = (-B - sqrt( B*B - 4.0*C*(A-R*R) ) ) / (2 * C);

		phi  = MAX(kappa1,kappa2);

		/* Sanity check. Should not get here !!!  */
		if ( phi < 0.0 ) {
			fprintf(stderr,"Material update error: "
					"negative phi at unloading:%f \n",phi);
			MPI_Abort(MPI_COMM_WORLD, ERROR);
			exit(1);
		}

		beta = phi * 0.50 / G;
		kn   = Pegasus( beta,  el_cnt);  // this is a check
		*Err    = ( 1.0 + kn - beta ) * getHardening(el_cnt, kn) / G - 3.0 * beta;

		if ( fabs(*Err) > theErrorTol  ) {
			fprintf(stdout," Warning --- UnloadingError/ErrorTolerance= %f/%f \n", fabs(*Err), theErrorTol );
		}

		if ( kn < 0  ) {
			fprintf(stdout," =*=*=*=* CHECK FOR UNSTABLE BEHAVIOR =*=*=*=* \n"
					"Found negative kappa at unloading.  k=%f  \n", kn );
			MPI_Abort(MPI_COMM_WORLD, ERROR);
			exit(1);
		}

	} else {
		fprintf(stdout," =*=*=*=* CHECK FOR UNSTABLE BEHAVIOR =*=*=*=* \n"
				"Cannot compute kappa at unloading. \n" );
		MPI_Abort(MPI_COMM_WORLD, ERROR);
		exit(1);
	}

	return kn;
}


/*===============================================================*/
/*===============================================================*/



/*   Material update function for Frederick-Armstrong, and Frederick-Armstrong-Modified models    */
void MatUpd_vMFA (double J2_pr, tensor_t dev_pr, double psi, double Su, tensor_t eta_n, tensor_t e_n1, double mu, double Lambda, double Sy,
		tensor_t *epl, tensor_t ep, double *ep_bar, double ep_barn, tensor_t *eta, tensor_t *sigma, tensor_t stresses,
		double *fs,  double *psi_n, double *loadunl_n, double *Tao_n, double *Tao_max ) {

	double C1, C2, C3, C4, C5, S_ss, S_aa, S_sa, dl, H_kin, H_nlin, G1;
	int i;

	// update psi value
	if ( ( theMaterialModel == VONMISES_FAM ) && *psi_n == 0.0 )
		*psi_n = psi ;

    S_ss    = 2.0 * J2_pr;
    S_aa    = 2.0 * tensor_J2 ( eta_n ); /* eta_n is already deviatoric */
    S_sa    = 2.0 * combtensor_J2(eta_n, dev_pr);

    if ( theMaterialModel == VONMISES_FAM ) {
    	H_kin  = (*psi_n) * mu;
    } else {
    	H_kin  = psi * mu;
    }

     H_nlin = H_kin/( sqrt(8.0/3.0) * Su - Sy );  /*Sy was already scaled by srt(2) before calling MatUpd_vMKH() */
     G1     = mu + H_kin/2.0;

    /* coefficients of the quartic function  */
    // Based on Auricchio and Tylor (1995) formulation
    C1 = pow( 2.0 * mu * H_nlin, 2.0 );
    C2 = (4.0 * Sy * mu * H_nlin + 8.0 * mu * G1 ) * H_nlin;
    C3 = ( H_nlin * H_nlin ) * ( Sy * Sy - S_ss ) + 4.0 * G1 * G1 + 4.0 * H_nlin * Sy * ( mu + G1 );
    C4 = 2.0 * H_nlin * ( Sy * Sy + S_sa - S_ss ) + 4.0 * Sy * G1;
    C5 = Sy*Sy - S_aa + 2.0 * S_sa - S_ss;

    /* roots finder */
    double coeff[5] = { C5, C4, C3, C2, C1 };
    double z[8];

    gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc (5);
    gsl_poly_complex_solve (coeff, 5, w, z);
    gsl_poly_complex_workspace_free (w);

    /* find the minimum positive root */
    dl  = FLT_MAX;

    for (i = 0; i < 4; i++) {
    	if ( ( z[2*i] >= 0.0 ) && ( fabs(z[2*i+1]) <= 1E-10 ) && ( z[2*i] < dl ) )
    		dl = z[2*i];
    }

    /* Sanity check. Should not get here !!!  */
    if ( dl == FLT_MAX ) {
        fprintf(stderr,"Material update error: "
                "could not find a positive root for von Mises with kinematic hardening\n");
        MPI_Abort(MPI_COMM_WORLD, ERROR);
        exit(1);
    }

    double T_lambda = 1.0 / ( 1.0 + H_nlin * dl );
    tensor_t Zo = scaled_tensor( eta_n, T_lambda );
    tensor_t Z  = subtrac_tensors( dev_pr, Zo );

    double magZ = sqrt( 2.0 * tensor_J2 ( Z ) );
    tensor_t   n = scaled_tensor( Z, 1.0 / magZ );

    /* updated variables */
    *epl    = add_tensors(ep, scaled_tensor( n, dl ) );              /* Updated plastic strains  */
    *ep_bar = ep_barn + dl;                                          /* Updated plastic equivalent strain  */
    *eta    = add_tensors( scaled_tensor( eta_n, T_lambda ), scaled_tensor( n, H_kin * T_lambda * dl ) ); /* Updated backstress tensor */

    *sigma  = subtrac_tensors( stresses, scaled_tensor( n , 2.0 * mu * dl ) );  /* Updated stresses tensor */

    /*  updated yield function. Must be ZERO */
    tensor_t Sdev = subtrac_tensors( dev_pr, scaled_tensor( n , 2.0 * mu * dl ) );
    *fs           = sqrt( 2.0 * tensor_J2( subtrac_tensors( Sdev,*eta ) ) ) - Sy;


	/* ================================================================ */
	/*  check for unloading. Only for vonMises_Modified (VONMISES_FAM)  */
	if ( theMaterialModel == VONMISES_FAM ) {

	    double loadunl = 2.0 * combtensor_J2(n, *eta);
	    double Tao_v   = sqrt(  tensor_J2( Sdev )  );


	    if (  ((*loadunl_n) * loadunl < 0.0) &&  ((*Tao_n) >= (*Tao_max))  ){

            *Tao_max  =  *Tao_n;

            // get elastic J2 at t-1
        	tensor_t stresses_n1    = point_stress ( e_n1, mu, Lambda );    /* elastic stress at t-1   */

        	/* compute the invariants of the elastic stress at t-1*/
        	double   I1_n1   = tensor_I1 ( stresses_n1 );
        	double   oct_n1  = tensor_octahedral ( I1_n1 );
        	tensor_t dev_n1  = tensor_deviator ( stresses_n1, oct_n1 );

        	double   Tao_e   = sqrt( tensor_J2 ( dev_n1 ) );

        	Su           = Su * 2.0/sqrt(3.0);
            *psi_n       = log( ( Su + *Tao_max ) / ( Su - *Tao_max ) ) * Su / ( Tao_e - *Tao_max );

            H_kin  = (*psi_n) * mu;
            H_nlin = H_kin/( sqrt(2.0) * Su - Sy );  /*Sy was already scaled by srt(2) before calling MatUpd_vMKH() */
            G1     = mu + H_kin/2.0;

            /* coefficients of the quartic function  */
            // Based on Auricchio and Tylor (1995) formulation
            C1 = pow( 2.0 * mu * H_nlin, 2.0 );
            C2 = (4.0 * Sy * mu * H_nlin + 8.0 * mu * G1 ) * H_nlin;
            C3 = ( H_nlin * H_nlin ) * ( Sy * Sy - S_ss ) + 4.0 * G1 * G1 + 4.0 * H_nlin * Sy * ( mu + G1 );
            C4 = 2.0 * H_nlin * ( Sy * Sy + S_sa - S_ss ) + 4.0 * Sy * G1;
            C5 = Sy*Sy - S_aa + 2.0 * S_sa - S_ss;

            /* roots finder */
            double coeff2[5] = { C5, C4, C3, C2, C1 };
            double z2[8];

            gsl_poly_complex_workspace * w2 = gsl_poly_complex_workspace_alloc (5);
            gsl_poly_complex_solve (coeff2, 5, w2, z2);
            gsl_poly_complex_workspace_free (w2);

            /* find the minimum positive root */
            dl  = FLT_MAX;

            for (i = 0; i < 4; i++) {
            	if ( ( z2[2*i] >= 0.0 ) && ( fabs(z2[2*i+1]) <= 1E-10 ) && ( z2[2*i] < dl ) )
            		dl = z2[2*i];
            }

            /* Sanity check. Should not get here !!!  */
            if ( dl == FLT_MAX ) {
                fprintf(stdout, "Tao_max=%f, Tao_e=%f, psi_n=%f,  Su=%f, C1=%f, C2=%f, C3=%f, C4=%f, C5=%f", *Tao_max, Tao_e, *psi_n, Su, C1, C2, C3, C4, C5);
                fprintf(stderr,"Material update error: "
                        "could not find a positive root for von Mises with kinematic hardening\n");
                MPI_Abort(MPI_COMM_WORLD, ERROR);
                exit(1);
            }

            double T_lambda = 1.0 / ( 1.0 + H_nlin * dl );
            tensor_t Zo = scaled_tensor( eta_n, T_lambda );
            tensor_t Z  = subtrac_tensors( dev_pr, Zo );

            double magZ = sqrt( 2.0 * tensor_J2 ( Z ) );
            tensor_t   n = scaled_tensor( Z, 1.0 / magZ );

            /* updated variables */
            *epl    = add_tensors(ep, scaled_tensor( n, dl ) );              /* Updated plastic strains  */
            *ep_bar = ep_barn + dl;                                          /* Updated plastic equivalent strain  */
            *eta    = add_tensors( scaled_tensor( eta_n, T_lambda ), scaled_tensor( n, H_kin * T_lambda * dl ) ); /* Updated backstress tensor */

            *sigma  = subtrac_tensors( stresses, scaled_tensor( n , 2.0 * mu * dl ) );  /* Updated stresses tensor */

            /*  updated yield function. Must be ZERO */
            tensor_t Sdev = subtrac_tensors( dev_pr, scaled_tensor( n , 2.0 * mu * dl ) );
            *fs           = sqrt( 2.0 * tensor_J2( subtrac_tensors( Sdev,*eta ) ) ) - Sy;
            loadunl       = 2.0 * combtensor_J2(n, *eta);
            Tao_v         = sqrt(  tensor_J2( Sdev )  );
            *Tao_max      = MAX(Tao_v, (*Tao_n) );

	    } else {

	        *Tao_max = MAX(Tao_v, (*Tao_max) );

	    }


	    *loadunl_n = loadunl; // update load-unload condition
	    *Tao_n     = Tao_v; // update tao_n


	}

}


void material_update ( nlconstants_t constants, tensor_t e_n, tensor_t e_n1, tensor_t ep, tensor_t eta_n,  double ep_barn, tensor_t sigma0, double dt,
		tensor_t *epl, tensor_t *eta, tensor_t *sigma, double *ep_bar, double *fs, double *psi_n, double *loadunl_n, double *Tao_n, double *Tao_max, double *kp, tensor_t *sigma_ref,
		int *flagTolSubSteps, int *flagNoSubSteps, double *ErrBA) {
	/* INPUTS:
	 * constants: Material constants
	 * e_n      : Total strain tensor
	 * e_n1     : Total strain tensor at t-1
	 * ep       : Plastic strain tensor at t-1. Used to compute the predictor state
	 * ep_barn  : equivalent plastic strain at t-1. Used to compute the predictor hardening function
	 * sigma0   : Approximated self-weight tensor.
	 * dt       : Time step
	 * eta_n    : backstress tensor at t-a. Used to compute the predictor state in vonMises kinematics
	 *
	 * OUTPUTS:
	 * fs       : Updated yield function value
	 * epl      : Updated plastic strain
	 * sigma    : Updated stress tensor
	 * eta      : Updated backstress tensor
	 * ep_bar   : Updated equivalent hardening variable
	 */

	double c, h, kappa, mu, Sy, beta, alpha, gamma, phi, dil, Fs_pr, Lambda, dLambda=0.0,
		   Tol_sigma = 1e-05, cond1, cond2, psi0, m;

	h      = constants.h;
	c      = constants.c;

	phi    = constants.phi;
	dil    = constants.dil_angle;

	mu     = constants.mu;
	Lambda = constants.lambda;
	kappa  = constants.lambda + 2.0 * mu / 3.0;

	beta   = constants.beta;
	alpha  = constants.alpha;
	gamma  = constants.gamma;
	Sy     = constants.Sstrain0*mu;

	psi0   = constants.psi0;
	m      = constants.m;

	//phi_pt = gamma / (3.0*beta);

	/* ---------      get the stress predictor tensor      ---------*/
	tensor_t estrain     = subtrac_tensors ( e_n, ep );             /* strain predictor   */
	tensor_t stresses    = point_stress ( estrain, mu, Lambda );    /* stress predictor   */
	tensor_t sigma_trial = add_tensors(stresses,sigma0);            /* stress predictor TOTAL  */

	/* compute the invariants of the stress predictor*/
	double   I1_pr   = tensor_I1 ( sigma_trial );
	double   oct_pr  = tensor_octahedral ( I1_pr );
	tensor_t dev_pr  = tensor_deviator ( sigma_trial, oct_pr );
    dev_pr           = subtrac_tensors ( dev_pr, eta_n );        /* Subtract backstress tensor  */

    double   J2_pr   = tensor_J2 ( dev_pr );
	double   J3_pr   = tensor_J3 (dev_pr);
	tensor_t dfds_pr = compute_dfds ( dev_pr, J2_pr, beta );
	vect1_t  sigma_ppal;

	/*  check predictor state */
	Fs_pr = compute_yield_surface_stateII ( J3_pr, J2_pr, I1_pr, alpha, phi, sigma_trial) - compute_hardening(gamma,c,Sy, h,ep_barn,phi, psi0); /* Fs predictor */

	if ( Fs_pr < Tol_sigma ) {
		*epl    = copy_tensor(ep);
		*sigma  = copy_tensor(stresses); /* return stresses without self-weight  */
		*ep_bar = ep_barn;
		*fs     = Fs_pr;
		return;
	}

	// corrector stage
	if ( ( theMaterialModel == VONMISES_EP ) || ( theMaterialModel == DRUCKERPRAGER ) ){
		Sy   = 0.0;
		psi0 = 0.0;

		dLambda = Fs_pr / ( mu + 9.0 * kappa * alpha * beta + h * gamma * gamma ); // ep_bar=gamma*dLambda as in Souza ec., (8.106)

		/* Updated plastic strains */
		epl->xx = ep.xx +  dLambda * dfds_pr.xx;
		epl->yy = ep.yy +  dLambda * dfds_pr.yy;
		epl->zz = ep.zz +  dLambda * dfds_pr.zz;
		epl->xy = ep.xy +  dLambda * dfds_pr.xy;
		epl->yz = ep.yz +  dLambda * dfds_pr.yz;
		epl->xz = ep.xz +  dLambda * dfds_pr.xz;

		/* Updated stresses */
		estrain     = subtrac_tensors ( e_n, *epl );
		stresses    = point_stress ( estrain, mu, Lambda );
		*sigma      = stresses;
		stresses    = add_tensors ( stresses, sigma0 );

		/* updated invariants*/
		double I1     = tensor_I1 ( stresses );
		double oct    = tensor_octahedral ( I1 );
		tensor_t dev  = tensor_deviator ( stresses, oct );
		double J2     = tensor_J2 ( dev );
		double J3     = tensor_J3 ( dev );

		/* Updated equivalent plastic strain */
		*ep_bar = ep_barn + dLambda * gamma;

		/* Updated yield function  */
		*fs = compute_yield_surface_stateII ( J3, J2, I1, alpha, phi, stresses) - compute_hardening(gamma,c,Sy, h,*ep_bar,phi, psi0);

		/* check for apex zone in DP model */
		if (  theMaterialModel == DRUCKERPRAGER  ){

			double Imin = I1_pr - 9.0 * kappa * beta / mu * sqrt(J2_pr);
			double Imax = I1_pr + sqrt(J2_pr)/alpha;

			if ( (I1 < Imin) || (I1 > Imax) || (sqrt(J2_pr) - mu * dLambda < 0.0) ) { /*return to the apex  */

				dLambda = (  compute_yield_surface_stateII ( 0.0, 0.0, I1_pr, alpha, phi, sigma_trial ) - compute_hardening(gamma,c,Sy,h,ep_barn,phi, psi0) ) / ( 9.0 * kappa * alpha * beta  + h * gamma * gamma );

				/* Updated equivalent plastic strain */
				*ep_bar = ep_barn + dLambda * gamma;

				double Skk = I1_pr - 9.0 * kappa * beta * dLambda;

				/* Updated stresses:
				 * It must be isotropic at the apex*/
				stresses.xx    = Skk/3.0;
				stresses.yy    = Skk/3.0;
				stresses.zz    = Skk/3.0;
				stresses.xy    = 0.0;
				stresses.xz    = 0.0;
				stresses.yz    = 0.0;
				*sigma      = subtrac_tensors ( stresses, sigma0 ); /* Sigma is still isotropic since sigma0 is isotropic */

				double Skk_rlt = tensor_I1 ( *sigma );  /* Relative trace. */

				/* Updated strains */
				estrain.xx = Skk_rlt / (3.0 * kappa);
				estrain.yy = Skk_rlt / (3.0 * kappa);
				estrain.zz = Skk_rlt / (3.0 * kappa);
				estrain.xy = 0.0;
				estrain.xz = 0.0;
				estrain.yz = 0.0;

				*epl  = subtrac_tensors ( e_n, estrain );

				*fs = alpha * Skk - compute_hardening(gamma,c,Sy,h,*ep_bar,phi,psi0);
			}
		}
	} else if ( theMaterialModel == VONMISES_FAM || theMaterialModel == VONMISES_FA ) {

		/* compute coefficients of the quartic function */
		dev_pr = add_tensors ( dev_pr, eta_n );       /* restore deviatoric predictor    */
		J2_pr   = tensor_J2 ( dev_pr );
		Sy      = sqrt(2.0)*Sy;                      // scale Sy to comply with the formulation for vonMises kinematic

		MatUpd_vMFA ( J2_pr,  dev_pr,  psi0,  c,  eta_n,  e_n1, mu,  Lambda,  Sy, epl,  ep,  ep_bar,  ep_barn,  eta,  sigma,  stresses, fs,  psi_n,  loadunl_n,  Tao_n,  Tao_max );
		return;

	}  else if ( theMaterialModel == VONMISES_BAE || theMaterialModel == VONMISES_BAH || theMaterialModel == VONMISES_GQH ) {

		MatUpd_vMGeneral ( constants,  kp,  e_n,  e_n1, sigma_ref, sigma, flagTolSubSteps, flagNoSubSteps, ErrBA );
		return;

	} else { /* Must be MohrCoulomb soil */

		/* Spectral decomposition of the sigma_trial tensor*/
		vect1_t   n1, n2, n3, sigma_ppal_trial;
		tensor_t  stressRecomp;
		int edge;

		if (theTensionCutoff == YES) {
			specDecomp(sigma_trial, &n1, &n2, &n3, &sigma_ppal_trial); /* eig_values.x > eig_values.y > eig_values.z   */
			double ST_fnc = get_ShearTensionLimits (phi, c, sigma_ppal_trial.x , sigma_ppal_trial.z); /* shear-tension function */

			if ( (ST_fnc > 0) && ( sigma_ppal_trial.x > 0 ) ) {
				/* Perform tension-cutoff  */
				TensionCutoff_Return( kappa, mu, phi, c, sigma_ppal_trial, &sigma_ppal );

				/* get updated stress tensor "sigma" */
				stressRecomp = specRecomp(sigma_ppal, n1, n2, n3);
				*sigma  = subtrac_tensors(stressRecomp,sigma0);

				estrain = elastic_strains (*sigma, mu, kappa);
				*epl  = subtrac_tensors ( e_n, estrain );

				/* updated invariants*/
				*fs = sigma_ppal.x;

				*ep_bar = 0.0; /* Todo: should think in a correct way to compute it when the tension cutoff option is on  */
				return;

			} else if ( (ST_fnc > 0) && ( sigma_ppal_trial.x <= 0 ) ) { /* This is an elastic state in the tension zone  */
				*epl    = copy_tensor(ep);
				*sigma  = copy_tensor(stresses); /* return stresses without self-weight  */
				*ep_bar = ep_barn;
				*fs     = sigma_ppal_trial.x;
				return;
			}

			/* set the spectral decomposition flag to 1 to avoid double computation  */
			// flagSpecDec = 1;
		}

		Fs_pr = compute_yield_surface_stateII ( J3_pr, J2_pr, I1_pr, alpha, phi, sigma_trial) - compute_hardening(gamma,c,Sy,h,ep_barn,phi,psi0);

		if ( Fs_pr < 0.0 ) {
			*epl    = copy_tensor(ep);
			*sigma  = copy_tensor(stresses); /* return stresses without self-weight  */
			*ep_bar = ep_barn;
			*fs     = Fs_pr;
			return;
		}

		/* Spectral decomposition */
		if ( theTensionCutoff == NO )
			specDecomp(sigma_trial, &n1, &n2, &n3, &sigma_ppal_trial); /* eig_values.x > eig_values.y > eig_values.z   */

		/* Return to the main plane */
		BOX85_l(ep_barn, sigma_ppal_trial, phi, dil, h, c, kappa, mu, &sigma_ppal, ep_bar); /* Return to the main plan */

		/* Check assumption of returning to the main plan */
		if ( ( sigma_ppal.x <= sigma_ppal.y ) || ( sigma_ppal.y <= sigma_ppal.z ) ) {

			if ( ( 1. - sin(dil) ) * sigma_ppal_trial.x - 2. * sigma_ppal_trial.y + ( 1. + sin(dil) ) * sigma_ppal_trial.z > 0 )
				edge = 1; /*return to the right edge*/
			else
				edge = -1; /*return to the left edge*/

			BOX86_l(ep_barn, sigma_ppal_trial, phi, dil, h, c, kappa, mu, edge, &sigma_ppal, ep_bar);

			cond1 = sigma_ppal.x - sigma_ppal.y;
			cond2 = sigma_ppal.y - sigma_ppal.z;
			double p_trial = ( sigma_ppal_trial.x + sigma_ppal_trial.y + sigma_ppal_trial.z )/3.0;

			if ( (cond1 <= 0.0 ) && ( fabs(cond1) >= Tol_sigma)  ) { /* return to the apex */
				if (theTensionCutoff == YES) {
					sigma_ppal.x = 0.0;
					sigma_ppal.y = 0.0;
					sigma_ppal.z = 0.0;
					*ep_bar      = 0.0; /* Todo: should think in a correct way to compute it when the tension cutoff option is on  */
				} else
					BOX87_l(ep_barn, p_trial, phi, dil, h, c, kappa, &sigma_ppal, ep_bar);
			} else if ( (cond2 <= 0.0) && (fabs(cond2) >= Tol_sigma) ){
				if (theTensionCutoff == YES) {
					sigma_ppal.x = 0.0;
					sigma_ppal.y = 0.0;
					sigma_ppal.z = 0.0;
					*ep_bar      = 0.0; /* Todo: should think in a correct way to compute it when the tension cutoff option is on  */
				} else
					BOX87_l(ep_barn, p_trial, phi, dil, h, c, kappa, &sigma_ppal, ep_bar);
			}
		}

		/*		 Tension cutoff check
		if ( (theTensionCutoff == YES) && ( sigma_ppal.x > 0 ) )  { Todo: Please check this one more time. Doriam

			 Perform tension-cutoff
			TensionCutoff_Return( kappa, mu, phi, c, sigma_ppal_trial, &sigma_ppal );

			 get updated stress tensor "sigma"
			stressRecomp = specRecomp(sigma_ppal, n1, n2, n3);
		 *sigma  = subtrac_tensors(stressRecomp,sigma0);

			estrain = elastic_strains (*sigma, mu, kappa);
		 *epl  = subtrac_tensors ( e_n, estrain );

			 updated invariants
		 *fs = sigma_ppal.x;
			return;
		}
		 Done tension cutoff check  */

		/* get updated stress tensor "sigma" */
		stressRecomp = specRecomp(sigma_ppal, n1, n2, n3);
		*sigma  = subtrac_tensors(stressRecomp,sigma0);

		estrain = elastic_strains (*sigma, mu, kappa);
		*epl  = subtrac_tensors ( e_n, estrain );

		/* updated invariants*/
		double I1     = tensor_I1 ( stressRecomp );
		double oct    = tensor_octahedral ( I1 );
		tensor_t dev  = tensor_deviator ( stressRecomp, oct );
		double J2     = tensor_J2 ( dev );
		double J3     = tensor_J3 ( dev );

		if ( (theTensionCutoff == YES) && ( sigma_ppal.x >= 0 ) )
			*fs = 0;
		else
			*fs = compute_yield_surface_stateII ( J3, J2, I1, alpha, phi, stressRecomp) - compute_hardening(gamma,c,Sy,h,*ep_bar,phi, psi0);

	}



	//	if ( thePlasticityModel == RATE_DEPENDANT ) { /*TODO: Add implementation for rate dependant material model  */
	//
	//		/* Rate dependant material is considered as a Drucker-Prager material   */
	//
	//		double factor      = fs / constants.k;
	//		double strainRate  = constants.strainrate;
	//		double sensitivity = constants.sensitivity;
	//		double oneOverM    = 1.0 / sensitivity;
	//
	//		dLambda =	strainRate * pow(factor, oneOverM);
	//
	//		epl->xx = ep.xx + dt * dLambda * dfds.xx;
	//		epl->yy = ep.yy + dt * dLambda * dfds.yy;
	//		epl->zz = ep.zz + dt * dLambda * dfds.zz;
	//		epl->xy = ep.xy + dt * dLambda * dfds.xy;
	//		epl->yz = ep.yz + dt * dLambda * dfds.yz;
	//		epl->xz = ep.xz + dt * dLambda * dfds.xz;
	//
	//		estrain   = subtrac_tensors ( e_n, *epl );
	//		*sigma    = point_stress ( estrain, constants.mu, constants.lambda );
	//
	//		*ep_bar   = ep_barn + gamma * dLambda;
	//
	//
	//	}


}

/* computes the derivatives of the flow potential for a Drucker-Prager material */
tensor_t compute_dfds (tensor_t dev, double J2, double beta) {

    tensor_t dfds;

    dfds.xx = dev.xx / ( 2.0 * sqrt(J2) ) + beta;
    dfds.yy = dev.yy / ( 2.0 * sqrt(J2) ) + beta;
    dfds.zz = dev.zz / ( 2.0 * sqrt(J2) ) + beta;
    dfds.xy = dev.xy / ( 2.0 * sqrt(J2) );
    dfds.yz = dev.yz / ( 2.0 * sqrt(J2) );
    dfds.xz = dev.xz / ( 2.0 * sqrt(J2) );

    return dfds;

}

/*tensor_t compute_pstrain2 ( nlconstants_t constants, tensor_t pstrain1, tensor_t tstrain,
							tensor_t dfds, double dLambda, double dt, double J2, double I1,
							double J2_st, double I1_st, double po ) {

    tensor_t pstrain2;
    double kappa;

	kappa  = ( constants.lambda + 2.0 * constants.mu / 3.0 );

	if ( dLambda == 0 )
		return pstrain1;

    if ( thePlasticityModel == RATE_DEPENDANT ) {
    	pstrain2.xx = pstrain1.xx + dt * dLambda * dfds.xx;
    	pstrain2.yy = pstrain1.yy + dt * dLambda * dfds.yy;
    	pstrain2.zz = pstrain1.zz + dt * dLambda * dfds.zz;
    	pstrain2.xy = pstrain1.xy + dt * dLambda * dfds.xy;
    	pstrain2.yz = pstrain1.yz + dt * dLambda * dfds.yz;
    	pstrain2.xz = pstrain1.xz + dt * dLambda * dfds.xz;

    } else if ( po >  0.0  ) {

        pstrain2.xx = tstrain.xx -  ( po - I1_st / 3. ) / ( 3.0 * kappa );
        pstrain2.yy = tstrain.yy -  ( po - I1_st / 3. ) / ( 3.0 * kappa );
        pstrain2.zz = tstrain.zz -  ( po - I1_st / 3. ) / ( 3.0 * kappa );
        pstrain2.xy = tstrain.xy;
        pstrain2.yz = tstrain.yz;
        pstrain2.xz = tstrain.xz;

    } else {
    	pstrain2.xx = pstrain1.xx +  dLambda * dfds.xx;
    	pstrain2.yy = pstrain1.yy +  dLambda * dfds.yy;
    	pstrain2.zz = pstrain1.zz +  dLambda * dfds.zz;
        pstrain2.xy = pstrain1.xy +  dLambda * dfds.xy;
        pstrain2.yz = pstrain1.yz +  dLambda * dfds.yz;
        pstrain2.xz = pstrain1.xz +  dLambda * dfds.xz;
    }

    return pstrain2;
}*/

int get_displacements(mysolver_t *solver, elem_t *elemp, fvector_t *u) {

    int i;
    int res = 0;

    /* Capture displacements for each node */
    for (i = 0; i < 8; i++) {

        int32_t    lnid;
        fvector_t *dis;

        lnid = elemp->lnid[i];
        dis  = solver->tm1 + lnid;

        res += vector_is_all_zero( dis );

        u[i].f[0] = dis->f[0];
        u[i].f[1] = dis->f[1];
        u[i].f[2] = dis->f[2];

    }

    return res;
}

/* -------------------------------------------------------------------------- */
/*                              Stability methods                             */
/* -------------------------------------------------------------------------- */

void check_yield_limit(mesh_t *myMesh, int32_t eindex, double vs, double fs,
                       double k, int qp)
{
    if ( fs > 1.5 * k ) {

        if ( superflag > 0 ) {

            double  north_m, east_m, depth_m;
            int32_t lnid0;

            lnid0 = myMesh->elemTable[eindex].lnid[0];

            north_m = (myMesh->ticksize)*(double)myMesh->nodeTable[lnid0].x;
            east_m  = (myMesh->ticksize)*(double)myMesh->nodeTable[lnid0].y;
            depth_m = (myMesh->ticksize)*(double)myMesh->nodeTable[lnid0].z;

            fprintf(stderr, "\n\n\tcompute_nonlinear_entities:"
                    "\n\tAn element exceeded the yield surface."
                    "\n\tThe element origin is at:\n"
                    "\n\tnorth (m) = %f"
                    "\n\teast  (m) = %f"
                    "\n\tdepth (m) = %f"
                    "\n\tVs  (m/s) = %f"
                    "\n\tFs        = %f"
                    "\n\tk         = %f"
                    "\n\tqp        = %d\n"
                    "\n\tA smaller dt or coarser mesh is required.\n",
                    north_m, east_m, depth_m, vs, fs, k, qp);

            MPI_Abort(MPI_COMM_WORLD, ERROR);
            exit(1);

        } else {

            superflag++;
        }
    }

    return;
}

void check_strain_stability ( double dLambda, double dt, mesh_t *myMesh,
                              int32_t eindex, double vs, double fs,
                              double k, int qp)
{
    double STRAINRATELIMIT = 0.002;

    if ( dLambda > STRAINRATELIMIT / dt ) {

        if ( superflag > 0 ) {

            double  north_m, east_m, depth_m;
            int32_t lnid0;

            lnid0 = myMesh->elemTable[eindex].lnid[0];

            north_m = (myMesh->ticksize)*(double)myMesh->nodeTable[lnid0].x;
            east_m  = (myMesh->ticksize)*(double)myMesh->nodeTable[lnid0].y;
            depth_m = (myMesh->ticksize)*(double)myMesh->nodeTable[lnid0].z;

            fprintf(stderr, "\n\n\tcompute_nonlinear_entities:"
                    "\n\tAn element violated dlambda condition"
                    "\n\tThe element origin is at:\n"
                    "\n\tnorth (m) = %f"
                    "\n\teast  (m) = %f"
                    "\n\tdepth (m) = %f"
                    "\n\tVs  (m/s) = %f"
                    "\n\tFs        = %f"
                    "\n\tk         = %f"
                    "\n\tqp        = %d\n"
                    "\n\tdLambda   = % 8e\n",
                    north_m, east_m, depth_m, vs, fs, k, qp, dLambda);

            MPI_Abort(MPI_COMM_WORLD, ERROR);
            exit(1);

        } else {

            superflag++;
        }
    }

    return;
}


/* -------------------------------------------------------------------------- */
/*                   Nonlinear core computational methods                     */
/* -------------------------------------------------------------------------- */
void BOX85_l(double ep_bar_n,vect1_t sigma_ppal_trial,double Phi, double Psi, double H, double c0, double K, double G, vect1_t *sigma_ppal, double *ep_bar_n1) {

/* Return mapping algorithm copied from:
   EA de Souza Neto, D Peric, D.O (2008). Computational methods for plasticity. Wiley */

	double a, dGamma;

	sigma_ppal->x = 0.0;
	sigma_ppal->y = 0.0;
	sigma_ppal->z = 0.0;

	a = ( 4. * G * ( 1. + 1./3. * sin(Psi) * sin(Phi) ) + 4. * K * sin(Phi) * sin(Psi) );

	dGamma = ( sigma_ppal_trial.x - sigma_ppal_trial.z + ( sigma_ppal_trial.x + sigma_ppal_trial.z ) * sin(Phi) - 2. * ( c0 + H * ep_bar_n ) * cos(Phi) ) /
			 ( 4. * H * cos(Phi) * cos(Phi) +  a );

	sigma_ppal->x = sigma_ppal_trial.x - ( 2. * G * ( 1. + 1./3. * sin(Psi) ) + 2. * K * sin(Psi) ) * dGamma;

	sigma_ppal->y = sigma_ppal_trial.y + ( 4./3. * G - 2. * K ) * sin(Psi) * dGamma;

	sigma_ppal->z = sigma_ppal_trial.z + ( 2. * G * ( 1. - 1./3. * sin(Psi) ) - 2. * K * sin(Psi) ) * dGamma;

	*ep_bar_n1 = ep_bar_n + 2.0 * cos(Phi) * dGamma;

}

void BOX86_l(double ep_bar_n,vect1_t sigma_ppal_trial,double Phi, double Psi, double H, double c0, double K, double G, double id, vect1_t *sigma_ppal, double *ep_bar_n1) {

	/* Return mapping algorithm copied from:
	   EA de Souza Neto, D Peric, D.O (2008). Computational methods for plasticity. Wiley */

	double aux, phi_a_bar, phi_b_bar, a, b, a11, b11, dGamma[2]={0}, sum_dGamma;

	sigma_ppal->x = 0.0;
	sigma_ppal->y = 0.0;
	sigma_ppal->z = 0.0;

    aux = 2. * cos(Phi) * ( c0 + H * ep_bar_n );

    phi_a_bar = sigma_ppal_trial.x - sigma_ppal_trial.z + ( sigma_ppal_trial.x + sigma_ppal_trial.z ) * sin(Phi) - aux;

    a = 4. * G * ( 1. + 1. / 3. * sin(Phi) * sin(Psi) ) + 4. * K * sin(Phi) * sin(Psi);

    if ( id == 1 ) {
        b = 2. * G * ( 1. + sin(Phi) + sin(Psi) - (1./3.) * sin(Phi) * sin(Psi) ) + 4. * K * sin(Phi) * sin(Psi);
        phi_b_bar = sigma_ppal_trial.x - sigma_ppal_trial.y + ( sigma_ppal_trial.x + sigma_ppal_trial.y ) * sin(Phi) - aux;
    } else {
        b = 2. * G * ( 1. - sin(Phi) - sin(Psi) - (1./3.) * sin(Phi) * sin(Psi) ) + 4. * K * sin(Phi) * sin(Psi);
        phi_b_bar = sigma_ppal_trial.y - sigma_ppal_trial.z + ( sigma_ppal_trial.y + sigma_ppal_trial.z ) * sin(Phi) - aux;
    }

    a11 = ( a + 4. * H * ( cos(Phi) * cos(Phi) ) );
    b11 = ( b + 4. * H * ( cos(Phi) * cos(Phi) ) );

    dGamma[0] = 1. / ( a11 * a11 - b11 * b11 ) * ( a11 * phi_a_bar - b11 * phi_b_bar );
    dGamma[1] = 1. / ( a11 * a11 - b11 * b11 ) * ( a11 * phi_b_bar - b11 * phi_a_bar );
    sum_dGamma    = dGamma[0] + dGamma[1];

   *ep_bar_n1 = ep_bar_n + 2. * cos(Phi) * sum_dGamma;

   if ( id == 1 ) {
    sigma_ppal->x = sigma_ppal_trial.x - ( 2. * G * ( 1. + (1./3.) * sin(Psi) ) + 2. * K * sin(Psi) ) * (sum_dGamma);
    sigma_ppal->y = sigma_ppal_trial.y + ( 4. * G / 3. - 2. * K ) * sin(Psi) * dGamma[0] + ( 2. * G * ( 1. - (1./3.) * sin(Psi) ) - 2.* K * sin(Psi) ) * dGamma[1];
    sigma_ppal->z = sigma_ppal_trial.z + ( 2. * G * ( 1. - (1./3.) * sin(Psi) ) - 2. * K * sin(Psi) ) * dGamma[0] + ( ( 4. * G / 3. ) - 2. * K ) * sin(Psi) * dGamma[1];
   } else {
    sigma_ppal->x = sigma_ppal_trial.x - ( 2. * G * ( 1. + (1./3.) * sin(Psi) ) + 2. * K * sin(Psi) ) * dGamma[0] + ( ( 4. * G / 3.) - 2. * K ) * sin(Psi) * dGamma[1];
    sigma_ppal->y = sigma_ppal_trial.y + ( 4. * G / 3. - 2. * K ) * sin(Psi) * dGamma[0] - ( 2. * G * ( 1. + (1./3.) * sin(Psi) ) + 2. * K * sin(Psi) ) * dGamma[1];
    sigma_ppal->z = sigma_ppal_trial.z + ( 2. * G * ( 1.- (1./3.) * sin(Psi) ) - 2. * K * sin(Psi) ) * sum_dGamma;
   }

}

void BOX87_l(double ep_bar_n,double p_trial,double Phi, double Psi, double H, double c0, double K, vect1_t *sigma_ppal, double *ep_bar_n1) {

	/* Return mapping algorithm copied from:
	   EA de Souza Neto, D Peric, D.O (2008). Computational methods for plasticity. Wiley */
	double omega, dep_bar, cot_phi, p;

	sigma_ppal->x = 0.0;
	sigma_ppal->y = 0.0;
	sigma_ppal->z = 0.0;

	omega = sin(Psi)/cos(Phi);
	cot_phi = cos(Phi) / sin (Phi);

	dep_bar = ( p_trial - ( c0 + H * ep_bar_n ) * cot_phi ) / ( H * cot_phi + omega * K); /*Here we deviate from the original formulation in
	                                                                                        order be able to use zero dilatancy  */
	*ep_bar_n1 = ep_bar_n + dep_bar;
	p = p_trial - K * omega * dep_bar;

	sigma_ppal->x = p;
	sigma_ppal->y = p;
	sigma_ppal->z = p;

}

tensor_t specRecomp(vect1_t eig_val, vect1_t n1, vect1_t n2, vect1_t n3) {
	tensor_t stress = zero_tensor();

	/* from first eigen_vector */
	stress.xx += eig_val.x * ( n1.x * n1.x);
	stress.yy += eig_val.x * ( n1.y * n1.y);
	stress.zz += eig_val.x * ( n1.z * n1.z);
	stress.xy += eig_val.x * ( n1.x * n1.y);
	stress.xz += eig_val.x * ( n1.x * n1.z);
	stress.yz += eig_val.x * ( n1.y * n1.z);

	/* from 2nd eigen_vector */
	stress.xx += eig_val.y * ( n2.x * n2.x);
	stress.yy += eig_val.y * ( n2.y * n2.y);
	stress.zz += eig_val.y * ( n2.z * n2.z);
	stress.xy += eig_val.y * ( n2.x * n2.y);
	stress.xz += eig_val.y * ( n2.x * n2.z);
	stress.yz += eig_val.y * ( n2.y * n2.z);

	/* from 3rd eigen_vector */
	stress.xx += eig_val.z * ( n3.x * n3.x);
	stress.yy += eig_val.z * ( n3.y * n3.y);
	stress.zz += eig_val.z * ( n3.z * n3.z);
	stress.xy += eig_val.z * ( n3.x * n3.y);
	stress.xz += eig_val.z * ( n3.x * n3.z);
	stress.yz += eig_val.z * ( n3.y * n3.z);

	return stress;

}



/*      Computes the spectral decomposition of a 3x3 symmetric matrix         */
/* -------------------------------------------------------------------------- */
/* Eigen decomposition code for symmetric 3x3 matrices, copied from the public
 domain Java Matrix library JAMA. */
// #define MAX(a, b) ((a)>(b)?(a):(b))

double hypot2(double x, double y) {
    return sqrt(x*x+y*y);
}

// Symmetric Householder reduction to tridiagonal form.

void tred2(double V[3][3], double *d, double *e) {

    int n=3, j, i, k;
    //  This is derived from the Algol procedures tred2 by
    //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutine in EISPACK.

    for (j = 0; j < n; j++) {
        d[j] = V[n-1][j];
    }

    // Householder reduction to tridiagonal form.

    for (i = n-1; i > 0; i--) {

        // Scale to avoid under/overflow.

        double scale = 0.0;
        double h = 0.0;
        for (k = 0; k < i; k++) {
            scale = scale + fabs(d[k]);
        }
        if (scale == 0.0) {
            e[i] = d[i-1];
            for ( j = 0; j < i; j++) {
                d[j] = V[i-1][j];
                V[i][j] = 0.0;
                V[j][i] = 0.0;
            }
        } else {

            // Generate Householder vector.

            for (k = 0; k < i; k++) {
                d[k] /= scale;
                h += d[k] * d[k];
            }
            double f = d[i-1];
            double g = sqrt(h);
            if (f > 0) {
                g = -g;
            }
            e[i] = scale * g;
            h = h - f * g;
            d[i-1] = f - g;
            for ( j = 0; j < i; j++) {
                e[j] = 0.0;
            }

            // Apply similarity transformation to remaining columns.

            for ( j = 0; j < i; j++) {
                f = d[j];
                V[j][i] = f;
                g = e[j] + V[j][j] * f;
                for ( k = j+1; k <= i-1; k++) {
                    g += V[k][j] * d[k];
                    e[k] += V[k][j] * f;
                }
                e[j] = g;
            }
            f = 0.0;
            for ( j = 0; j < i; j++) {
                e[j] /= h;
                f += e[j] * d[j];
            }
            double hh = f / (h + h);
            for ( j = 0; j < i; j++) {
                e[j] -= hh * d[j];
            }
            for ( j = 0; j < i; j++) {
                f = d[j];
                g = e[j];
                for ( k = j; k <= i-1; k++) {
                    V[k][j] -= (f * e[k] + g * d[k]);
                }
                d[j] = V[i-1][j];
                V[i][j] = 0.0;
            }
        }
        d[i] = h;
    }

    // Accumulate transformations.

    for ( i = 0; i < n-1; i++) {
        V[n-1][i] = V[i][i];
        V[i][i] = 1.0;
        double h = d[i+1];
        if (h != 0.0) {
            for ( k = 0; k <= i; k++) {
                d[k] = V[k][i+1] / h;
            }
            for ( j = 0; j <= i; j++) {
                double g = 0.0;
                for ( k = 0; k <= i; k++) {
                    g += V[k][i+1] * V[k][j];
                }
                for ( k = 0; k <= i; k++) {
                    V[k][j] -= g * d[k];
                }
            }
        }
        for ( k = 0; k <= i; k++) {
            V[k][i+1] = 0.0;
        }
    }
    for ( j = 0; j < n; j++) {
        d[j] = V[n-1][j];
        V[n-1][j] = 0.0;
    }
    V[n-1][n-1] = 1.0;
    e[0] = 0.0;

}



void tql2(double V[3][3], double* d, double* e) {

    //  This is derived from the Algol procedures tql2, by
    //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutine in EISPACK.
    int i, l, k, j,  n=3;


    for ( i = 1; i < n; i++) {
        e[i-1] = e[i];
    }
    e[n-1] = 0.0;

    double f = 0.0;
    double tst1 = 0.0;
    double eps = pow(2.0,-52.0);
    for ( l = 0; l < n; l++) {

        // Find small subdiagonal element

        tst1 = MAX(tst1,fabs(d[l]) + fabs(e[l]));
        int m = l;
        while (m < n) {
            if (fabs(e[m]) <= eps*tst1) {
                break;
            }
            m++;
        }

        // If m == l, d[l] is an eigenvalue,
        // otherwise, iterate.

        if (m > l) {
            int iter = 0;
            do {
                iter = iter + 1;  // (Could check iteration count here.)

                // Compute implicit shift

                double g = d[l];
                double p = (d[l+1] - g) / (2.0 * e[l]);
                double r = hypot2(p,1.0);
                if (p < 0) {
                    r = -r;
                }
                d[l] = e[l] / (p + r);
                d[l+1] = e[l] * (p + r);
                double dl1 = d[l+1];
                double h = g - d[l];
                for ( i = l+2; i < n; i++) {
                    d[i] -= h;
                }
                f = f + h;

                // Implicit QL transformation.

                p = d[m];
                double c = 1.0;
                double c2 = c;
                double c3 = c;
                double el1 = e[l+1];
                double s = 0.0;
                double s2 = 0.0;
                for ( i = m-1; i >= l; i--) {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = hypot2(p,e[i]);
                    e[i+1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i+1] = h + s * (c * g + s * d[i]);

                    // Accumulate transformation.

                    for ( k = 0; k < n; k++) {
                        h = V[k][i+1];
                        V[k][i+1] = s * V[k][i] + c * h;
                        V[k][i] = c * V[k][i] - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;

                // Check for convergence.

            } while (fabs(e[l]) > eps*tst1);
        }
        d[l] = d[l] + f;
        e[l] = 0.0;
    }

    // Sort eigenvalues and corresponding vectors.

    for ( i = 0; i < n-1; i++) {
         k = i;
        double p = d[i];
        for ( j = i+1; j < n; j++) {
            if (d[j] < p) {
                k = j;
                p = d[j];
            }
        }
        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for ( j = 0; j < n; j++) {
                p = V[j][i];
                V[j][i] = V[j][k];
                V[j][k] = p;
            }
        }
    }
}


void  specDecomp(tensor_t sigma, vect1_t *n1, vect1_t *n2, vect1_t *n3, vect1_t *eig_values){
  double V[3][3], d[3], e[3];

    V[0][0] = sigma.xx;
    V[0][1] = sigma.xy;
    V[0][2] = sigma.xz;
    V[1][0] = sigma.xy;
    V[1][1] = sigma.yy;
    V[1][2] = sigma.yz;
    V[2][0] = sigma.xz;
    V[2][1] = sigma.yz;
    V[2][2] = sigma.zz;

    tred2(V, d, e);
    tql2(V, d, e);

    eig_values->x = d[2];
    eig_values->y = d[1];
    eig_values->z = d[0];

    n3->x = V[0][0];
    n3->y = V[1][0];
    n3->z = V[2][0];

    n2->x = V[0][1];
    n2->y = V[1][1];
    n2->z = V[2][1];

    n1->x = V[0][2];
    n1->y = V[1][2];
    n1->z = V[2][2];

}

double get_ShearTensionLimits (double phi, double coh, double S1, double S3) {
/* Shear-tension limits defined as in FLAC3D User's manual version 5.01 pp 1-31   */

	double Nphi    = ( 1 + sin(phi) )/( 1 - sin(phi) );
	double alpha_p = sqrt( 1 + Nphi * Nphi )+ Nphi;
	double Sp      = -2 * coh * sqrt(Nphi);

	double h = S1 + alpha_p * ( S3 - Sp );

	return h;

}

void TensionCutoff_Return( double k, double mu, double phi, double coh, vect1_t sigma_ppal_pr, vect1_t *SigmaUP ) {

double N_phi, S_p, alpha_1, alpha_2, Phi_pr[5], Pl1, Pl2, aa, bb, det,
       DLam1, DLam2;
int    i, pos,  zone;

/* Compute maximum compression stress at sigma1 = 0 */
N_phi = (1+sin(phi))/(1-sin(phi));
S_p   = - 2. * coh * sqrt(N_phi);

/* Material constants */
alpha_1 = k + 4. * mu / 3.;
alpha_2 = k - 2. * mu / 3.;

/* check active zones in the predictor state */
vect1_t  sigma_TC1, sigma_TC2;
zone =  CornerZones( sigma_ppal_pr, S_p, &sigma_TC1, Phi_pr);

if ( zone == 0 ) {

    // get active planes. Note that the first plane is already an active plane
    Pl1 = Phi_pr[0];
    Pl2 = 0.0;
    for (i = 1; i < 5; i++) {
    	if ( Phi_pr[i] > 0 ) {
    		Pl2 = Phi_pr[i];
    		pos=i;
    		break;
    	}
    }

    if ( Pl2 != 0.0 ) {
    	aa = alpha_1;
    	bb = alpha_2;

    	if ( ( pos == 2 ) || ( pos == 5 ) )
    		bb=-bb;

    	det = ( aa * aa - bb * bb );

    	DLam1 = (  aa * Pl1 - bb * Pl2 ) / det;
    	DLam2 = ( -bb * Pl1 + aa * Pl2 ) / det;

    	switch ( pos ) {
    	case ( 2 ):
			sigma_TC2.x  = sigma_ppal_pr.x - ( alpha_1 * DLam1 - alpha_2 * DLam2  );
    		sigma_TC2.y  = sigma_ppal_pr.y - ( alpha_2 * DLam1 - alpha_2 * DLam2  );
    		sigma_TC2.z  = sigma_ppal_pr.z - ( alpha_2 * DLam1 - alpha_1 * DLam2  );
    		break;

    	case ( 3 ):
			sigma_TC2.x = sigma_ppal_pr.x - ( alpha_1 * DLam1 + alpha_2 * DLam2  );
    		sigma_TC2.y = sigma_ppal_pr.y - ( alpha_2 * DLam1 + alpha_1 * DLam2  );
    		sigma_TC2.z = sigma_ppal_pr.z - ( alpha_2 * DLam1 + alpha_2 * DLam2  );
    		break;

    	case ( 4 ):
			sigma_TC2.x = sigma_ppal_pr.x - ( alpha_1 * DLam1 + alpha_2 * DLam2  );
    		sigma_TC2.y = sigma_ppal_pr.y - ( alpha_2 * DLam1 + alpha_2 * DLam2  );
    		sigma_TC2.z = sigma_ppal_pr.z - ( alpha_2 * DLam1 + alpha_1 * DLam2  );
    		break;

    	case ( 5 ):
			sigma_TC2.x = sigma_ppal_pr.x - ( alpha_1 * DLam1 - alpha_2 * DLam2  );
    		sigma_TC2.y = sigma_ppal_pr.y - ( alpha_2 * DLam1 - alpha_1 * DLam2  );
    		sigma_TC2.z = sigma_ppal_pr.z - ( alpha_2 * DLam1 - alpha_2 * DLam2  );
    		break;
    	}

    } else {
    	sigma_TC2.x = sigma_ppal_pr.x - Pl1;
    	sigma_TC2.y = sigma_ppal_pr.y - ( alpha_2 * Pl1 / alpha_1 );
    	sigma_TC2.z = sigma_ppal_pr.z - ( alpha_2 * Pl1 / alpha_1 );
    }

    /* Check active zones for sigma tension cut 2 "sigma_TC2"  */
    zone =  CornerZones( sigma_TC2, S_p, SigmaUP, Phi_pr);

    if ( zone == 0 ) {
        SigmaUP->x = sigma_TC2.x;
        SigmaUP->y = sigma_TC2.y;
        SigmaUP->z = sigma_TC2.z;

        if( SigmaUP->z < S_p )
        	SigmaUP->z = S_p;

        if(SigmaUP->y < S_p)
        	SigmaUP->y = S_p;
    }

    return;

}

	SigmaUP->x = sigma_TC1.x;
	SigmaUP->y = sigma_TC1.y;
	SigmaUP->z = sigma_TC1.z;
}


int CornerZones( vect1_t Sigma, double S_p, vect1_t* SigmaUP, double Phi_pr[5]) {

	double Tol;
	int i, cnt=0, zone;

Phi_pr[0] = Sigma.x;
Phi_pr[1] = S_p - Sigma.z;
Phi_pr[2] = Sigma.y;
Phi_pr[3] = Sigma.z;
Phi_pr[4] = S_p - Sigma.y;

Tol=-1.0e-10;

/* sanity check. Cannot exist more that 3 active surfaces */
for (i = 0; i < 5; i++) {
	if ( Phi_pr[i] > 0 )
		++cnt;
}
if (cnt > 3) {
	fprintf(stderr, "Error: %d: No more that 3 active surfaces can coexist ",cnt);
	MPI_Abort(MPI_COMM_WORLD,ERROR);
	exit(1);
}

if ( ( Phi_pr[2] > Tol ) && ( Phi_pr[3] > Tol ) ) {
	SigmaUP->x = 0.0;
	SigmaUP->y = 0.0;
	SigmaUP->z = 0.0;
	zone=1;
} else if ( ( Phi_pr[3] > Tol ) && ( Phi_pr[4] > Tol ) ) {
	SigmaUP->x = 0.0;
	SigmaUP->y = S_p;
	SigmaUP->z = 0.0;
	zone=3;
} else if ( ( Phi_pr[1] > Tol ) && ( Phi_pr[4] > Tol ) ) {
	SigmaUP->x = 0.0;
	SigmaUP->y = S_p;
	SigmaUP->z = S_p;
	zone=5;
} else if ( ( Phi_pr[1] > Tol ) && ( Phi_pr[2] > Tol) ) {
	SigmaUP->x = 0.0;
	SigmaUP->y = 0.0;
	SigmaUP->z = S_p;
	zone=7;
} else {
	SigmaUP->x = 0.0;
	SigmaUP->y = 0.0;
	SigmaUP->z = 0.0;
	zone=0;
}

return zone;

}




/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */


double smooth_rise_factor(int32_t step, double dt) {

    /* TODO: Consider to eliminate t1, i.e. t1 = 0 */

    static noyesflag_t preloaded = NO;
    static double n1, n2, n3, n22, n31;
    static double C1, C2, B1, B2;
    static int    N;

    /* Pre-load constants, executed only once */
    if ( preloaded == NO ) {

        N = (int)(theGeostaticLoadingT / dt);

        n1 = (int)( 0.1 * N );
        n2 = (int)( 0.5 * N );
        n3 = (int)( 0.9 * N );

        n31 = n3 - n1;

        C1 = 2.0 / ( n31 * (n2 - n1) ) ;
        C2 = 2.0 / ( n31 * (n2 - n3) ) ;

        B1 = 0.5 * n1 * n1;
        B2 = 0.5 * ( n31*(n2-n3) + n3*n3 );

        preloaded = YES;
    }

    /* Started dynamic simulation (the most common case) */
    if ( step > n3 ) {
        return 1.0;
    }

    /* Has not even really started (a zeros-buffer zone) */
    if ( step <= n1 ) {
        return 0.0;
    }

    /* The actual geostatic loading cases */

    n22 = 0.5 * step * step;

    if ( (step > n1) && (step <= n2) ) {

        return C1 * (n22 - step*n1 + B1);

    } else if ( (step > n2) && (step <= n3) ) {

        return C2 * (n22 - step*n3 + B2);
    }

    /* The code should never get here */
    fprintf(stderr, "Smooth Rise Error %d %d\n", N, step);
    MPI_Abort(MPI_COMM_WORLD, ERROR);
    exit(1);
}

void add_force_reactions ( mesh_t     *myMesh,
                           mysolver_t *mySolver )
{

    int       i;
    int32_t   beindex, eindex;

    /* Loop on the number of elements */
    for (beindex = 0; beindex < myBottomElementsCount; beindex++) {

        elem_t    *elemp;

        eindex = myBottomElements[beindex].element_id;
        elemp  = &myMesh->elemTable[eindex];

        for (i = 4; i < 8; i++) {

            int32_t    lnid;
            fvector_t *nodalForce;

            lnid = elemp->lnid[i];
            nodalForce = mySolver->force + lnid;

            nodalForce->f[2] += myBottomElements[beindex].nodal_force[i-4];

        } /* element nodes */
    }

    return;
}

void check_balance( int32_t myID ) {

    double totalReaction = 0;
    int i;
    int32_t   beindex;

    for (beindex = 0; beindex < myBottomElementsCount; beindex++) {
        for ( i = 0; i < 4; i++ ) {
            totalReaction += myBottomElements[beindex].nodal_force[i];
        }
    }

    double theTotalWeight;
    double theTotalReaction;

    MPI_Reduce( &totalWeight, &theTotalWeight, 1,
                MPI_DOUBLE, MPI_SUM, 0, comm_solver );
    MPI_Reduce( &totalReaction, &theTotalReaction, 1,
                MPI_DOUBLE, MPI_SUM, 0, comm_solver );

    if ( myID == 0 ) {
        fprintf( stdout,
                 "\tTotal Weight   = %20.6f\n"
                 "\tTotal Reaction = %20.6f\n"
                 "\tDifference     = %20.6f\n\n",
                 theTotalWeight, theTotalReaction,
                 theTotalWeight+theTotalReaction );
    }

    return;
}

void compute_addforce_gravity( mesh_t     *myMesh,
                               mysolver_t *mySolver,
                               int         step,
                               double      dt )
{
    /*
     * Some gravity values:
     *
     * Los Angeles    9.796 m/s2      Athens    9.800 m/s2
     * Mexico City    9.779 m/s2      Auckland  9.799 m/s2
     * San Francisco  9.800 m/s2      Rome      9.803 m/s2
     * Istanbul       9.808 m/s2      Tokyo     9.798 m/s2
     * Vancouver      9.809 m/s2      Ottawa    9.806 m/s2
     *
     */

    #define G 9.8

    int32_t   eindex;

    /* Loop on the number of elements */
    for (eindex = 0; eindex < myMesh->lenum; eindex++) {

        int      i;
        elem_t  *elemp;
        edata_t *edata;
        double   h, h3;
        double   rho, W;

        /* Capture element data structure */
        elemp = &myMesh->elemTable[eindex];
        edata = (edata_t *)elemp->data;

        /* capture element data */
        h   = (double)edata->edgesize;
        rho = (double)edata->rho;

        /* Compute nodal total weight contribution */
        h3 = h * h * h;
        W  = h3 * rho * G * 0.125; /* volume x density x gravity / 8 */

        if ( step == theGeostaticFinalStep+10 ) {
            totalWeight += W*dt*dt*8;
        }

        /* Loop over the 8 element nodes:
         * Add the gravitational force contribution calculated with respect
         * to the current time-step to the nodal force vector.
         */
        for (i = 0; i < 8; i++) {

            int32_t    lnid;
            fvector_t *nodalForce;

            lnid       = elemp->lnid[i];
            nodalForce = mySolver->force + lnid;

            /* Force due to gravity is positive in Z-axis */
            nodalForce->f[2] += W * smooth_rise_factor(step, dt) * dt * dt;

        } /* element nodes */

    }

    if ( step > theGeostaticFinalStep ) {
        add_force_reactions(myMesh, mySolver);
    }

    return;
}

void compute_bottom_reactions ( mesh_t     *myMesh,
                                mysolver_t *mySolver,
                                fmatrix_t (*theK1)[8],
                                fmatrix_t (*theK2)[8],
                                int         step,
                                double      dt )
{
    if ( step != theGeostaticFinalStep ) {
        return;
    }

    int32_t   beindex;
    int32_t   eindex;
    fvector_t localForce[8];

    for ( beindex = 0; beindex < myBottomElementsCount; beindex++ ) {

        int        i, j;
        elem_t    *elemp;
        e_t*       ep;

        eindex = myBottomElements[beindex].element_id;
        elemp  = &myMesh->elemTable[eindex];
        ep     = &mySolver->eTable[eindex];

        /* -------------------------------
         * Ku DONE IN THE CONVENTIONAL WAY
         * ------------------------------- */

        /* step 1: calculate the force due to the element stiffness */
        memset( localForce, 0, 8 * sizeof(fvector_t) );

        /* contribution by node j to node i (only for nodes at the bottom) */
        for (i = 4; i < 8; i++) {

            fvector_t* toForce = &localForce[i];

            for (j = 0; j < 8; j++) {

                int32_t    nodeJ  = elemp->lnid[j];
                fvector_t *myDisp = mySolver->tm1 + nodeJ;

                MultAddMatVec( &theK1[i][j], myDisp, ep->c1, toForce );
                MultAddMatVec( &theK2[i][j], myDisp, ep->c2, toForce );
            }
        }

        /* step 2: store the forces */
        for (i = 4; i < 8; i++) {
            myBottomElements[beindex].nodal_force[i-4] = localForce[i].f[2];
        }

        edata_t *edata;
        double   h, h3;
        double   rho, W;

        edata = (edata_t *)elemp->data;
        h     = (double)edata->edgesize;
        rho   = (double)edata->rho;
        h3    = h * h * h;
        W     = h3 * rho * G * 0.125;

        for (i = 0; i < 4; i++) {
            myBottomElements[beindex].nodal_force[i] -= W * dt * dt;
        }
    }

    return;
}

void geostatic_displacements_fix( mesh_t     *myMesh,
                                  mysolver_t *mySolver,
                                  double      totalDomainDepth,
                                  double      dt,
                                  int         step )
{

    if ( step > theGeostaticFinalStep ) {
        return;
    }

    int32_t nindex;

    for ( nindex = 0; nindex < myMesh->nharbored; nindex++ ) {

        double z_m = (myMesh->ticksize)*(double)myMesh->nodeTable[nindex].z;

        if ( z_m == totalDomainDepth ) {
            fvector_t *tm2Disp;
            tm2Disp = mySolver->tm2 + nindex;
            tm2Disp->f[2] = 0;
        }
    }

    return;
}

/*
 * compute_addforce_nl: Adds a fictitious force that accounts for
 *                      nonlinearity in the soil --- a correction
 *                      to the 'trial' stress within K-elastic
 *
 * The computation done here corresponds to calculate the last expression
 * in the manuscript prepared by Ricardo (part I) for the implementation of
 * nonlinear soil into the program.
 *
 *      Integral(grad(Phi) * Cijkl * PlasticStrain) * delta_t^2
 */
void compute_addforce_nl (mesh_t     *myMesh,
                          mysolver_t *mySolver,
                          double      theDeltaTSquared)
{

    int       i, j;
    int32_t   eindex;
    int32_t   nl_eindex;
    fvector_t localForce[8];

    /* Loop on the number of elements */
    for (nl_eindex = 0; nl_eindex < myNonlinElementsCount; nl_eindex++) {

        elem_t  *elemp;
        edata_t *edata;
        double   h, h3, WiJi;
        double   mu, lambda;

        eindex = myNonlinElementsMapping[nl_eindex];

        qptensors_t stresses;

        /* Capture the table of elements from the mesh and the size
         * This is what gives me the connectivity to nodes */
        elemp = &myMesh->elemTable[eindex];
        edata = (edata_t *)elemp->data;

        h    = (double)edata->edgesize;
        h3   = h * h * h;
        WiJi = h3 * 0.125; /* (h^3)/8 */

        nlconstants_t ec = myNonlinSolver->constants[nl_eindex];

        mu     = ec.mu;
        lambda = ec.lambda;

//        if ( theMaterialModel == LINEAR ) {

            stresses = myNonlinSolver->stresses[nl_eindex];

//        } else {
//
//            qptensors_t tstrains, pstrains, estrains;
//
//            tstrains = myNonlinSolver->strains   [nl_eindex];
//            pstrains = myNonlinSolver->pstrains1[nl_eindex];
//
//            /* compute total strain - plastic strain */
//            estrains = subtrac_qptensors(tstrains, pstrains);
//
//            /* compute the corresponding stress */
//            stresses = compute_qp_stresses(estrains, mu, lambda);
//        }

        /* Clean memory for the local force vector */
        memset(localForce, 0, 8 * sizeof(fvector_t));

        /* Loop over the 8 element nodes:
         * Calculates the forces on each vertex */
        for (i = 0; i < 8; i++) {

            fvector_t *toForce;

            /* Points the loop force to the element force */
            toForce = &localForce[i];

            /* Loop over the 8 quadrature points */
            for (j = 0; j < 8; j++) {

                double dx, dy, dz;

                compute_qp_dxi (&dx, &dy, &dz, i, j, h);

                /* Gauss integration: Sum(DeltaPsi * Sigma * wi * Ji) */

                toForce->f[0] += ( ( dx * stresses.qp[j].xx )
                                 + ( dy * stresses.qp[j].xy )
                                 + ( dz * stresses.qp[j].xz ) ) * WiJi;

                toForce->f[1] += ( ( dy * stresses.qp[j].yy )
                                 + ( dx * stresses.qp[j].xy )
                                 + ( dz * stresses.qp[j].yz ) ) * WiJi;

                toForce->f[2] += ( ( dz * stresses.qp[j].zz )
                                 + ( dy * stresses.qp[j].yz )
                                 + ( dx * stresses.qp[j].xz ) ) * WiJi;

            } /* quadrature points */

        } /* element nodes */

        /* Loop over the 8 element nodes:
         * Add the contribution calculated above to the node
         * forces carried from the source and stiffness.
         */
        for (i = 0; i < 8; i++) {

            int32_t    lnid;
            fvector_t *nodalForce;

            lnid = elemp->lnid[i];

            nodalForce = mySolver->force + lnid;

            nodalForce->f[0] -= localForce[i].f[0] * theDeltaTSquared;
            nodalForce->f[1] -= localForce[i].f[1] * theDeltaTSquared;
            nodalForce->f[2] -= localForce[i].f[2] * theDeltaTSquared;

        	if (  ( ( nodalForce->f[0] >=0 ) || ( nodalForce->f[0] < 0 ) ) &&
        		  ( ( nodalForce->f[1] >=0 ) || ( nodalForce->f[1] < 0 ) ) &&
        		  ( ( nodalForce->f[2] >=0 ) || ( nodalForce->f[2] < 0 ) ) ) {
        	} else {

        	}


        } /* element nodes */

    } /* all elements */

    return;
}

/*
 * compute_nonlinear_state: Compute the necessary quantities to determine
 *                             if an element has undergone plastic deformation.
 *
 * The steps described in this method are in reference to the manuscript
 * prepared by Ricardo about the implementation of nonlinear soil in the
 * program (Part III).
 *
 * The convention for the strain and stress tensors is:
 *
 *      T[6] = {Txx, Tyy, Tzz, Txy, Tyz, Txz}
 */
void compute_nonlinear_state ( mesh_t     *myMesh,
                               mysolver_t *mySolver,
                               int32_t     theNumberOfStations,
                               int32_t     myNumberOfStations,
                               station_t  *myStations,
                               double      theDeltaT,
                               int         step )
{
	/* In general, j-index refers to the quadrature point in a loop (0 to 7 for
	 * eight points), and i-index refers to the tensor component (0 to 5), with
	 * the following order xx[0], yy[1], zz[2], xy[3], yz[4], xz[5]. i-index is
	 * also some times used for the number of nodes (8, 0 to 7).
	 */

	int     i;
	int32_t eindex, nl_eindex;


	/* Loop over the number of local elements */
	for (nl_eindex = 0; nl_eindex < myNonlinElementsCount; nl_eindex++) {

		elem_t        *elemp;
		edata_t       *edata;
		nlconstants_t *enlcons;

		double         h;          /* Element edge-size in meters   */
		double         mu, lambda; /* Elasticity material constants */
		double         XI, QC;
		fvector_t      u[8];
		qptensors_t   *stresses, *tstrains, *tstrains1, *pstrains1, *pstrains2, *alphastress1, *alphastress2, *Sref;
		qpvectors_t   *epstr1, *epstr2,   *psi_n,   *lounlo_n,   *Sv_n,   *Sv_max, *kappa;

		/* Capture data from the element and mesh */
		eindex = myNonlinElementsMapping[nl_eindex];

		elemp = &myMesh->elemTable[eindex];
		edata = (edata_t *)elemp->data;
		h     = edata->edgesize;

		/* Capture data from the nonlinear element structure */

		enlcons = myNonlinSolver->constants + nl_eindex;

		mu     = enlcons->mu;
		lambda = enlcons->lambda;

		/* Capture the current state in the element */
		tstrains     = myNonlinSolver->strains      + nl_eindex;
		tstrains1     = myNonlinSolver->strains1    + nl_eindex;
		stresses     = myNonlinSolver->stresses     + nl_eindex;
		pstrains1    = myNonlinSolver->pstrains1    + nl_eindex;   /* Previous plastic tensor  */
		pstrains2    = myNonlinSolver->pstrains2    + nl_eindex;   /* Current  plastic tensor  */
		alphastress1 = myNonlinSolver->alphastress1 + nl_eindex;   /* Previous backstress tensor  */
		alphastress2 = myNonlinSolver->alphastress2 + nl_eindex;   /* Current  backstress tensor  */
		epstr1       = myNonlinSolver->ep1          + nl_eindex;
		epstr2       = myNonlinSolver->ep2          + nl_eindex;

		psi_n        = myNonlinSolver->psi_n        + nl_eindex;
		lounlo_n     = myNonlinSolver->LoUnlo_n     + nl_eindex;
		Sv_n         = myNonlinSolver->Sv_n         + nl_eindex;
		Sv_max       = myNonlinSolver->Sv_max       + nl_eindex;
		kappa        = myNonlinSolver->kappa        + nl_eindex;
		Sref         = myNonlinSolver->Sref         + nl_eindex;

		/* initialize kappa */
		if ( ( theMaterialModel == VONMISES_BAE  ||  theMaterialModel == VONMISES_BAH || theMaterialModel == VONMISES_GQH ) && ( step == 0 ) ){
			for (i = 0; i < 8; i++) {
				kappa->qv[i] = 1E+06;
			}
		}

		/* Capture displacements */
		if ( get_displacements(mySolver, elemp, u) == 0 ) {
			/* If all displacements are zero go for next element */
			continue;
		}

		/* Loop over the quadrature points */
		for (i = 0; i < 8; i++) {

			tensor_t  sigma0;

			/* Quadrature point local coordinates */
			double lx = xi[0][i] * qc ;
			double ly = xi[1][i] * qc ;
			double lz = xi[2][i] * qc ;

			tstrains1->qp[i]    = copy_tensor ( tstrains->qp[i] ); // get elastic strains at t-1
			/* Calculate total strains */
			tstrains->qp[i] = point_strain(u, lx, ly, lz, h);

			/* strain and backstress predictor  */
	        pstrains1->qp[i]    = copy_tensor ( pstrains2->qp[i] );     /* The strain predictor assumes that the current plastic strains equal those from the previous step   */
	        alphastress1->qp[i] = copy_tensor ( alphastress2->qp[i] );

			/* Calculate stresses */
			if ( ( theMaterialModel == LINEAR ) || ( step <= theGeostaticFinalStep ) ){
				stresses->qp[i]  = point_stress ( tstrains->qp[i], mu, lambda );
				continue;
			} else {

				if ( theApproxGeoState == YES )
					sigma0 = ApproxGravity_tensor(enlcons->sigmaZ_st, enlcons->phi, h, lz, edata->rho);
				else
					sigma0 = zero_tensor();

				int flagTolSubSteps=0, flagNoSubSteps=0;
				double ErrBA=0;

				/*				double po=90;
				if (i==5 && eindex == 111412 && ( step == 240 ) ) {
					po=89;
				}*/

				material_update ( *enlcons,           tstrains->qp[i],      tstrains1->qp[i],   pstrains1->qp[i],  alphastress1->qp[i], epstr1->qv[i],   sigma0,        theDeltaT,
						          &pstrains2->qp[i],  &alphastress2->qp[i], &stresses->qp[i],   &epstr2->qv[i],    &enlcons->fs[i],     &psi_n->qv[i],
						          &lounlo_n->qv[i], &Sv_n->qv[i], &Sv_max->qv[i], &kappa->qv[i], &Sref->qp[i], &flagTolSubSteps, &flagNoSubSteps, &ErrBA);

				if ( ( theMaterialModel == VONMISES_BAE || theMaterialModel == VONMISES_BAH || theMaterialModel == VONMISES_GQH ) ) {
					enlcons->fs[i] = ErrBA;
					//if (flagTolSubSteps==1)
					//	fprintf(stdout,"Exceeded Error Tolerance:%f at GP:%d, eindex: %d, step: %d \n", ErrBA, i, eindex, step);

					//if (flagNoSubSteps==1)
					//	fprintf(stdout,"Exceeded number of sub-steps at GP:%d, eindex: %d -- INCREASE SUBSTEPS NUMBER \n", i, eindex);

					//if (step>110)
						//fprintf(stdout,"stuck at GP:%d, eindex: %d, step: %d \n", i, eindex, step);
				}


			}
		} /* for all quadrature points */
	} /* for all nonlinear elements */
}

/* -------------------------------------------------------------------------- */
/*                        Nonlinear Finalize and Stats                        */
/* -------------------------------------------------------------------------- */

void nonlinear_yield_stats(mesh_t *myMesh, int32_t myID, int32_t theTotalSteps, int32_t theGroupSize) {

    static double VSMIN = 0;
    static double VSMAX = 10000;

    int32_t  nl_eindex, eindex;
    int      r;
    int      ranges = thePropertiesCount+1;
    double  *myFsMaxs;
    double  *myFsAvgs;
    double   vs, vs0, vs1;
    int32_t *myFsAvgCount;

    myFsMaxs     =  (double *)calloc(ranges, sizeof(double));
    myFsAvgs     =  (double *)calloc(ranges, sizeof(double));
    myFsAvgCount = (int32_t *)calloc(ranges, sizeof(int32_t));

    for ( nl_eindex = 0; nl_eindex < myNonlinElementsCount; nl_eindex++ ) {

        elem_t        *elemp;
        edata_t       *edata;
        nlconstants_t  ec;

        eindex = myNonlinElementsMapping[nl_eindex];
        elemp  = &myMesh->elemTable[eindex];
        edata  = (edata_t *)elemp->data;
        ec     = myNonlinSolver->constants[nl_eindex];
        vs     = edata->Vs;

        for ( r = 0; r < ranges; r++ ) {

            /* set bottom vs limit */
            if ( r == 0 ) {
                vs0 = VSMIN;
            } else {
                vs0 = theVsLimits[r-1];
            }

            /* set top vs limit */
            if ( r == ranges-1 ) {
                vs1 = VSMAX;
            } else {
                vs1 = theVsLimits[r];
            }

            if ( (vs > vs0) && (vs <= vs1) ) {
                myFsAvgs[r] += ec.avgFs;
                myFsAvgCount[r]++;
                if ( ec.maxFs > myFsMaxs[r] ) {
                    myFsMaxs[r] = ec.maxFs;
                }
            }
        }
    }

    if ( myNonlinElementsCount > 0 ) {
        for ( r = 0; r < ranges; r++ ) {
            myFsAvgs[r] /= theTotalSteps;
        }
    }


    double  *theFsMaxs     = NULL;
    double  *theFsAvgs     = NULL;
    int32_t *theFsAvgCount = NULL;

    theFsMaxs     =  (double *)calloc(ranges, sizeof(double));
    theFsAvgs     =  (double *)calloc(ranges, sizeof(double));
    theFsAvgCount = (int32_t *)calloc(ranges, sizeof(int32_t));

    MPI_Reduce(myFsMaxs, theFsMaxs, ranges, MPI_DOUBLE, MPI_MAX, 0, comm_solver );
    MPI_Reduce(myFsAvgs, theFsAvgs, ranges, MPI_DOUBLE, MPI_SUM, 0, comm_solver );
    MPI_Reduce(myFsAvgCount, theFsAvgCount, ranges, MPI_INT, MPI_SUM, 0, comm_solver );

    if ( myID == 0 ) {

        for ( r = 0; r < ranges; r++ ) {
            if ( theFsAvgCount[r] > 0 ) {
                theFsAvgs[r] /= theFsAvgCount[r];
            }
        }

        FILE *fp = hu_fopen( "stat-fs-yield.txt", "w" );

        fputs( "\n"
               "# ------------------------------------------- \n"
               "# Nonlinear Fs maximum and average values:    \n"
               "# ------------------------------------------- \n"
               "#   Vs >    Vs <=           Avg           Max \n"
               "# ------------------------------------------- \n", fp );

        for ( r = 0; r < ranges; r++ ) {

            /* set bottom vs limit */
            if ( r == 0 ) {
                vs0 = VSMIN;
            } else {
                vs0 = theVsLimits[r-1];
            }

            /* set top vs limit */
            if ( r == ranges-1 ) {
                vs1 = VSMAX;
            } else {
                vs1 = theVsLimits[r];
            }

            fprintf( fp, "%8.0f %8.0f % 10e % 10e\n",
                     vs0, vs1, theFsAvgs[r], theFsMaxs[r]);

        }

        fprintf( fp, "# ------------------------------------------- \n\n");

        hu_fclosep( &fp );

    }
}

/* -------------------------------------------------------------------------- */
/*                        Nonlinear Output to Stations                        */
/* -------------------------------------------------------------------------- */

void nonlinear_stations_init(mesh_t    *myMesh,
                             station_t *myStations,
                             int32_t    myNumberOfStations)
{

    if ( myNumberOfStations == 0 ) {
        return;
    }

    int32_t     eindex, nl_eindex;
    int32_t     iStation=0;
    vector3D_t  point;
    octant_t   *octant;
    int32_t     lnid0;

    myNumberOfNonlinStations = 0;
    for (iStation = 0; iStation < myNumberOfStations; iStation++) {

        for ( nl_eindex = 0; nl_eindex < myNonlinElementsCount; nl_eindex++ ) {

            /* capture the stations coordinates */
            point = myStations[iStation].coords;

            /* search the octant */
            if ( search_point(point, &octant) != 1 ) {
                fprintf(stderr,
                        "nonlinear_stations_init: "
                        "No octant with station coords\n");
                MPI_Abort(MPI_COMM_WORLD, ERROR);
                exit(1);
            }

            eindex = myNonlinElementsMapping[nl_eindex];

            lnid0 = myMesh->elemTable[eindex].lnid[0];

            if ( (myMesh->nodeTable[lnid0].x == octant->lx) &&
                 (myMesh->nodeTable[lnid0].y == octant->ly) &&
                 (myMesh->nodeTable[lnid0].z == octant->lz) ) {

                /* I have a match for the element's origin */

                /* Now, perform level sanity check */
                if (myMesh->elemTable[eindex].level != octant->level) {
                    fprintf(stderr,
                            "nonlinear_stations_init: First pass: "
                            "Wrong level of octant\n");
                    MPI_Abort(MPI_COMM_WORLD, ERROR);
                    exit(1);
                }

                myNumberOfNonlinStations++;

                break;
            }
        }
    }

    XMALLOC_VAR_N( myStationsElementIndices, int32_t, myNumberOfNonlinStations);
    XMALLOC_VAR_N( myNonlinStationsMapping, int32_t, myNumberOfNonlinStations);
 //   XMALLOC_VAR_N( myNonlinStations, nlstation_t, myNumberOfNonlinStations);

    int32_t nlStationsCount = 0;
    for (iStation = 0; iStation < myNumberOfStations; iStation++) {

        for ( nl_eindex = 0; nl_eindex < myNonlinElementsCount; nl_eindex++ ) {

            /* capture the stations coordinates */
            point = myStations[iStation].coords;

            /* search the octant */
            if ( search_point(point, &octant) != 1 ) {
                fprintf(stderr,
                        "nonlinear_stations_init: "
                        "No octant with station coords\n");
                MPI_Abort(MPI_COMM_WORLD, ERROR);
                exit(1);
            }

            eindex = myNonlinElementsMapping[nl_eindex];

            lnid0 = myMesh->elemTable[eindex].lnid[0];

            if ( (myMesh->nodeTable[lnid0].x == octant->lx) &&
                 (myMesh->nodeTable[lnid0].y == octant->ly) &&
                 (myMesh->nodeTable[lnid0].z == octant->lz) ) {

                /* I have a match for the element's origin */

                /* Now, perform level sanity check */
                if (myMesh->elemTable[eindex].level != octant->level) {
                    fprintf(stderr,
                            "nonlinear_stations_init: Second pass: "
                            "Wrong level of octant\n");
                    MPI_Abort(MPI_COMM_WORLD, ERROR);
                    exit(1);
                }

                if ( nlStationsCount >= myNumberOfNonlinStations ) {
                    fprintf(stderr,
                            "nonlinear_stations_init: Second pass: "
                            "More stations than initially counted\n");
                    MPI_Abort(MPI_COMM_WORLD, ERROR);
                    exit(1);
                }

                /* Store the element index and mapping to stations */
                myStationsElementIndices[nlStationsCount] = nl_eindex;
                myNonlinStationsMapping[nlStationsCount] = iStation;

                nlStationsCount++;

                break;
            }

        } /* for all my elements */

    } /* for all my stations */

/*    for ( iStation = 0; iStation < myNumberOfNonlinStations; iStation++ ) {

        tensor_t *stress, *strain, *pstrain1, *pstrain2;
        double   *ep1;

        strain   = &(myNonlinStations[iStation].strain);
        stress   = &(myNonlinStations[iStation].stress);
        pstrain1 = &(myNonlinStations[iStation].pstrain1);
        pstrain2 = &(myNonlinStations[iStation].pstrain2);
        ep1      = &(myNonlinStations[iStation].ep );
        *ep1     = 0.;

        init_tensorptr(strain);
        init_tensorptr(stress);
        init_tensorptr(pstrain1);
        init_tensorptr(pstrain2);

    }*/

}

void print_nonlinear_stations(mesh_t     *myMesh,
                              mysolver_t *mySolver,
                              station_t  *myStations,
                              int32_t     myNumberOfStations,
                              double      dt,
                              int         step,
                              int         rate)
{

    int32_t eindex;
    int32_t nl_eindex;
    int32_t iStation;
    int32_t mappingIndex;

    for ( iStation = 0; iStation < myNumberOfNonlinStations; iStation++ ) {
    	tensor_t       *stress, *tstrain, tstress;
    	qptensors_t    *stressF, *tstrainF;
    	double         bStrain = 0., bStress = 0., Fy, h;
    	tensor_t       sigma0;

    	elem_t         *elemp;
		edata_t        *edata;
    	nlconstants_t  *enlcons;

    	nl_eindex    = myStationsElementIndices[iStation];
    	eindex       = myNonlinElementsMapping[nl_eindex];
    	mappingIndex = myNonlinStationsMapping[iStation];
    	enlcons      = myNonlinSolver->constants + nl_eindex;

		elemp = &myMesh->elemTable[eindex];
		edata = (edata_t *)elemp->data;
		h     = edata->edgesize;

		/* compute the self-weight stresses at the first Gauss point*/
		double lz = -0.577350269189;
		if ( theApproxGeoState == YES )
			sigma0 = ApproxGravity_tensor(enlcons->sigmaZ_st, enlcons->phi, h, lz, edata->rho);
		else
			sigma0 = zero_tensor();

    	/* Capture data from the nonlinear element structure
    	 * corresponding to the first Gauss point*/
    	tstrainF   = myNonlinSolver->strains   + nl_eindex;
    	stressF    = myNonlinSolver->stresses  + nl_eindex;

    	stress      = &(stressF->qp[0]);            /* relative stresses of the first Gauss point */
    	tstress     = add_tensors(*stress,sigma0); /* compute the total stress tensor */

    	tstrain    = &(tstrainF->qp[0]);

    	Fy         = (myNonlinSolver->constants   + nl_eindex)->fs[0];

    	bStrain = tstrain->xx + tstrain->yy + tstrain->zz;
    	bStress = tstress.xx + tstress.yy + tstress.zz;

    	if (step % rate == 0) {
    		fprintf( myStations[mappingIndex].fpoutputfile,

    				" % 8e % 8e"
    				" % 8e % 8e"
    				" % 8e % 8e"
    				" % 8e % 8e"
    				" % 8e % 8e"
    				" % 8e % 8e"
    				" % 8e % 8e"
    				" % 8e",

    				tstrain->xx, tstress.xx, // 11 12
    				tstrain->yy, tstress.yy, // 13 14
    				tstrain->zz, tstress.zz, // 15 16
    				bStrain,     bStress,    // 17 18
    				tstrain->xy, tstress.xy, // 19 20
    				tstrain->yz, tstress.yz, // 21 22
    				tstrain->xz, tstress.xz,
    				Fy); // 23 24
    	}
    } /* for all my stations */

}
