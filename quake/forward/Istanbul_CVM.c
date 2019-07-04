/*
 * Istanbul_CVM.c
 *
 *  Created on: Oct 22, 2018
 *      Author: eafit
 */

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <stdint.h>
# include <string.h>
#include <inttypes.h>

# include "Istanbul_CVM.h"

#include "psolve.h"
#include "octor.h"
#include "util.h"
#include "cvm.h"
#include "geometrics.h"

static int      *n_layer2;
static int      *Layer_start_ID;

static double   *Zcoord_2;
static double   *Soil_Vs_data;
static double   *Soil_rho_data;
static double   *Soil_depth_data;

/*
#define ZCOORD2 Zcoord_2[4558] =

#define NLAYER2 n_layer2[4558] =

#define LAYERSTART_ID Layer_start_ID[4558] =

#define SOILVS_DATA Soil_Vs_data[28651] =

#define SOILRHO_DATA Soil_rho_data[28651] =

#define SOILDEPTH_DATA Soil_depth_data[28651] = {
*/


void Istanbul_init ( int32_t myID ) {

    /* Capturing data from file --- only done by PE0 */
    if (myID == 0) {
        if ( Istanbul_initparameters( ) != 0 ) {
            fprintf(stderr,"Thread 0: Istanbul_local_init: "
                    "Istanbul_initparameters error\n");
            MPI_Abort(MPI_COMM_WORLD, ERROR);
            exit(1);
        }
    }

    if (myID != 0) {
    	Zcoord_2        = (double*)malloc( sizeof(double) * 4558 );
    	n_layer2        = (int*)malloc( sizeof(int) * 4558 );
    	Layer_start_ID  = (int*)malloc( sizeof(int) * 4558 );

    	Soil_Vs_data    = (double*)malloc( sizeof(double) * 28651 );
    	Soil_rho_data   = (double*)malloc( sizeof(double) * 28651 );
    	Soil_depth_data = (double*)malloc( sizeof(double) * 28651 );
    }

    /* Broadcast table of properties */
    MPI_Bcast(Zcoord_2,         4558, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(n_layer2,         4558, MPI_INT, 0, comm_solver);
    MPI_Bcast(Layer_start_ID,   4558, MPI_INT, 0, comm_solver);

    MPI_Bcast(Soil_Vs_data,    28651, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(Soil_rho_data,   28651, MPI_DOUBLE, 0, comm_solver);
    MPI_Bcast(Soil_depth_data, 28651, MPI_DOUBLE, 0, comm_solver);

    return;

}

int32_t
Istanbul_initparameters ( ) {

	int i;
    FILE               *fp_Zcoord, *fp_nlay, *fp_Layerid, *fp_soilVs, *fp_soilDepth, *fp_soilrho;
    char                zcoord_file[256], nlayer_file[256], layerID_file[256],
    					soilVs_file[256], soilDepth_file[256], soilRho_file[256];

    /* read material info */
	sprintf( zcoord_file,    "%s/Zcoord.in",          "inputfiles" );
	sprintf( nlayer_file,    "%s/n_layer.in",         "inputfiles" );
	sprintf( layerID_file,   "%s/Layer_start_ID.in",  "inputfiles" );
	sprintf( soilVs_file,    "%s/Soil_Vs_data.in",    "inputfiles" );
	sprintf( soilDepth_file, "%s/Soil_depth_data.in", "inputfiles" );
	sprintf( soilRho_file,   "%s/Soil_rho_data.in",   "inputfiles" );

	if ( ( ( fp_Zcoord    = fopen ( zcoord_file ,    "r") ) == NULL ) ||
		 ( ( fp_nlay      = fopen ( nlayer_file ,    "r") ) == NULL ) ||
		 ( ( fp_Layerid   = fopen ( layerID_file ,   "r") ) == NULL ) ||
		 ( ( fp_soilVs    = fopen ( soilVs_file ,    "r") ) == NULL ) ||
		 ( ( fp_soilDepth = fopen ( soilDepth_file , "r") ) == NULL ) ||
		 ( ( fp_soilrho   = fopen ( soilRho_file ,   "r") ) == NULL ) ) {
	    fprintf(stderr, "Istanbul material data files not found \n" );
	    return -1;
	}

	Zcoord_2        = (double*)malloc( sizeof(double) * 4558 );
	n_layer2        = (int*)malloc( sizeof(int) * 4558 );
	Layer_start_ID  = (int*)malloc( sizeof(int) * 4558 );

	Soil_Vs_data    = (double*)malloc( sizeof(double) * 28651 );
	Soil_rho_data   = (double*)malloc( sizeof(double) * 28651 );
	Soil_depth_data = (double*)malloc( sizeof(double) * 28651 );

	if ( ( Zcoord_2        == NULL ) ||
		 ( n_layer2        == NULL ) ||
		 ( Layer_start_ID  == NULL ) ||
		 ( Soil_Vs_data    == NULL ) ||
		 ( Soil_rho_data   == NULL ) ||
		 ( Soil_depth_data == NULL ) ) {
		fprintf( stderr, "Error allocating transient arrays for Istanbul material data"
				"in Istanbul_initparameters " );
		return -1;
	}

	for ( i = 0; i < 4558; ++i) {
	    fscanf(fp_Zcoord,   " %lf ", &(Zcoord_2[i]) );
	    fscanf(fp_nlay,     " %d ", &(n_layer2[i]) );
	    fscanf(fp_Layerid,  " %d ", &(Layer_start_ID[i]) );
	}

	for ( i = 0; i < 28651; ++i) {
	    fscanf(fp_soilVs,    " %lf ", &(Soil_Vs_data[i]) );
	    fscanf(fp_soilDepth, " %lf ", &(Soil_depth_data[i]) );
	    fscanf(fp_soilrho,   " %lf ", &(Soil_rho_data[i]) );
	}

    fclose(fp_Zcoord);
    fclose(fp_nlay);
    fclose(fp_Layerid);
    fclose(fp_soilVs);
    fclose(fp_soilDepth);
    fclose(fp_soilrho);

    return 0;
}


void material_property_relative_V10_local(double x_input, double y_input, double z_input, double output[3] ) {

	int i;
	int k;

	//double ZCOORD2, SOILVS_DATA, SOILRHO_DATA, SOILDEPTH_DATA ;
	//int    NLAYER2, LAYERSTART_ID;

	double vs_layer, rho_layer, depth_layer;

	double Soil_depth;

	double zd_vs[4];
	double zd_vp[4];
	double zd_rho[4];

	int row1;

	double zi_elevation, zi_vs, zi_vp, zi_rho;

	double xmin = 388500.0;
	double xmax = 414750.0;
	double ymin = 4536400.0;
	double ymax = 4546900.0;

	int x_grid_ID, y_grid_ID, Rec_node_ID[4];

	double x_coord[4], y_coord[4], z_coord[4], xi, eta, N1, N2, N3, N4;

	if (z_input > 0.0) {

		output[0] = 0.0;
		output[1] = 0.0;
		output[2] = 0.0;

		return ;
	}

	if (!(x_input >= xmin && x_input <= xmax && y_input >= ymin && y_input <= ymax && z_input >= -300.0)) {

		output[0] = 3200.0;
		output[1] = 5500.0;
		output[2] = 2800.0;

		return ;
	}

	if (x_input == xmax) {
		x_grid_ID = 105;
	}
	else {
		x_grid_ID = floor((x_input - xmin) / 250.0)+1;
	}

	if (y_input == ymax) {
		y_grid_ID = 42;
	}
	else {
		y_grid_ID = floor((y_input - ymin) / 250.0)+1;
	}

	    Rec_node_ID[0] = 106*(y_grid_ID-1) +  x_grid_ID - 1;
	    Rec_node_ID[1] = 106*(y_grid_ID-1) +  x_grid_ID;
	    Rec_node_ID[2] = 106*y_grid_ID     +  x_grid_ID;
	    Rec_node_ID[3] = 106*y_grid_ID     +  x_grid_ID - 1;

	    x_coord[0] = xmin + (x_grid_ID-1)*250.0;
	    x_coord[1] = x_coord[0]+250.0;
	    x_coord[2] = x_coord[1];
	    x_coord[3] = x_coord[0];

	    y_coord[0] = ymin + (y_grid_ID-1)*250.0;
	    y_coord[1] = y_coord[0]+250.0;
	    y_coord[2] = y_coord[1];
	    y_coord[3] = y_coord[0];

	    z_coord[0] = Zcoord_2[Rec_node_ID[0]];
	    z_coord[1] = Zcoord_2[Rec_node_ID[1]];
	    z_coord[2] = Zcoord_2[Rec_node_ID[2]];
	    z_coord[3] = Zcoord_2[Rec_node_ID[3]];

	    xi  = (x_input - x_coord[0])/250.0*2.0 - 1.0;
	    eta = (y_input - y_coord[0])/250.0*2.0 - 1.0;

	    N1 = (1.0/4.0)*(1.0-xi)*(1.0-eta);
	    N2 = (1.0/4.0)*(1.0+xi)*(1.0-eta);
	    N3 = (1.0/4.0)*(1.0+xi)*(1.0+eta);
	    N4 = (1.0/4.0)*(1.0-xi)*(1.0+eta);

	    zi_elevation = N1*z_coord[0] + N2*z_coord[1] + N3*z_coord[2] + N4*z_coord[3];

		z_input = z_input + zi_elevation;

		row1 = 0;

		for (i = 0; i < 4; i++)
		{

			Soil_depth = 0.0;

			row1 = Layer_start_ID[Rec_node_ID[i]]-1;

			for (k = 0; k < n_layer2[Rec_node_ID[i]]; k++)
			{
				//printf("%f\n", Vs_Data[row1+k][0]);

				vs_layer     = Soil_Vs_data[row1+k];
				rho_layer    = Soil_rho_data[row1+k];
				depth_layer  = Soil_depth_data[row1+k];

				//printf("%lf %lf %lf %lf %lf\n", Soil_Property[0], Soil_Property[1], Soil_Property[2], Soil_Property[3], Soil_Property[4]);

				Soil_depth = Soil_depth + depth_layer;

				if (z_input > z_coord[i] - Soil_depth)
				{
					zd_vs[i] = vs_layer;
					zd_vp[i] = zd_vs[i] * 1.8;
					zd_rho[i] = rho_layer;
					break;
				}


				zd_vs[i] = vs_layer;
				zd_vp[i] = zd_vs[i] * 2.0;
				zd_rho[i] = rho_layer;

			}
			//printf("%d %d\n", row1,n_Layer[i][0]);
		}

		/*
		  Evaluate the interpolant.
		*/

		zi_vs = N1*zd_vs[0] + N2*zd_vs[1] + N3*zd_vs[2] + N4*zd_vs[3];
		zi_vp = N1*zd_vp[0] + N2*zd_vp[1] + N3*zd_vp[2] + N4*zd_vp[3];
		zi_rho = N1*zd_rho[0] + N2*zd_rho[1] + N3*zd_rho[2] + N4*zd_rho[3];

		if ( zi_vs <= 0.0 || zi_vp <= 0.0 || zi_rho <= 0.0  ) {
			 fprintf(stdout,"zero properties found at xm =%f, ym=%f, zm=%f\n", x_input, y_input, z_input);
			 zi_vs = 1000.00;
			 zi_vp = 2000.00;
			 zi_rho = 2.2;
		}

		output[0] = zi_vs;
		output[1] = zi_vp;
		output[2] = zi_rho*1000;

}
