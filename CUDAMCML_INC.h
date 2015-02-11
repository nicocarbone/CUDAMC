/*	This file is part of CUDAMCML_INC.

    CUDAMCML_INC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CUDAMCML_INC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CUDAMCML_INC.  If not, see <http://www.gnu.org/licenses/>.*/

// DEFINES
#define NUM_THREADS_PER_BLOCK 96//192 //Keep above 192 to eliminate global memory access overhead However, keep low to allow enough registers per thread
#define NUM_BLOCKS 15//14//56  //Keep numblocks a multiple of the #MP's of the GPU (8800GT=14MP)
#define NUM_THREADS 1440//1344//1792//896//1344//2688//896//10752//17920

#define NUMSTEPS_GPU 2000
#define PI 3.141592654f
#define RPI 0.318309886f
#define MAX_LAYERS 100
#define STR_LEN 200
#define MAX_STEP 14000
#define MAX_INCS 100
#define TAM_GRILLA 25
#define RAD_FIB_BAN 0.2
#define NUM_CAN_TEMP 1000
#define MAX_TEMP 50e-9
#define GANANCIA_TEMP 5 
#define TEMP_CAN (MAX_TEMP/(NUM_CAN_TEMP*GANANCIA_TEMP))
#define	LIGHTSPEED 2.997925E10 
#define GLASS_STEP 100.0f

//#define WEIGHT 0.0001f
#define WEIGHTI 429497u //0xFFFFFFFFu*WEIGHT
#define CHANCE 0.1f


// TYPEDEFS
typedef struct __align__(16)
{
	float z_min;		// Layer z_min [cm]
	float z_max;		// Layer z_max [cm]
	float mutr;			// Reciprocal mu_total [cm]
	float mua;			// Absorption coefficient [1/cm]
	float g;			// Anisotropy factor [-]
	float n;			// Refractive index [-]
}LayerStruct;

typedef struct __align__(16) 
{

	float x;		// Global x coordinate [cm]
	float y;		// Global y coordinate [cm]
	float z;		// Global z coordinate [cm]
	double t;		// Tiempo recorrido
	float dx;		// (Global, normalized) x-direction
	float dy;		// (Global, normalized) y-direction
	float dz;		// (Global, normalized) z-direction
	unsigned int weight;			// Photon weight
	int layer;				// Current layer
	//unsigned int* xx;		// Array de posicion en cada step, coordenada x
	//unsigned int* zz;		// Array de posicion en cada step, coordenada z
	unsigned int step;			// Step actual
	int inc;				// Current inclusion
}PhotonStruct;

typedef struct __align__(16)
{
	float dx;		// Detection grid resolution, x-direction [cm]
	float dy;		// Detection grid resolution, y-direction [cm]
	//float dz;		// Detection grid resolution, z-direction [cm]
	
	int nx;			// Number of grid elements in x-direction
	int ny;			// Number of grid elements in y-direction
	//int nz;			// Number of grid elements in z-direction
	//float sep;		// Separacion fibra de detaccion - eje optico
	float fix;		// Posición en x de la fibra de detección
	float fiy;		// Posición en y de la fibra de detección
	unsigned int face; 	// Face where the detection fiber is locarted
}DetStruct;

typedef struct //__align__(16)
{
	float x; 			// Inclusion's x coordinate
	float y; 			// Inclusion's y coordinate
	float z; 			// Inclusion's z coordinate
	float r; 			// Inclusion's radius
	float mutr; 			// Mu_total reciproco de la inclusion
	float mua; 			// Absorption coefficient of the inclusion
	float g;			// Anisotropy coefficient
	float n;			// Refractive index
	unsigned int layer;		// Inclusion layer
	unsigned int type;		// Inclusion type		
}IncStruct;

typedef struct 
{
	unsigned long number_of_photons;
	int ignoreAdetection;
	unsigned int n_layers;
	unsigned int start_weight;
	char outp_filename[STR_LEN];
	char inp_filename[STR_LEN];
	unsigned int n_inclusions;
	long begin,end;
	char AorB;
	DetStruct det;
	LayerStruct* layers;
	IncStruct* inclusion;
	float esp;
	float fx; 			// Posición en x de la fuente
	float fy;			// Posición en y de la fuente

}SimulationStruct;


typedef struct 
{
	PhotonStruct* p;					// Pointer to structure array containing all the photon data
	unsigned long long* x;				// Pointer to the array containing all the WMC x's
	unsigned int* a;					// Pointer to the array containing all the WMC a's
	unsigned int* thread_active;		// Pointer to the array containing the thread active status
	unsigned int* num_terminated_photons;	//Pointer to a scalar keeping track of the number of terminated photons

	unsigned long long* Rd_xy;			// Matriz 2D for reflexion
	//unsigned long long* A_xyz;			// Matriz 3D for absorcion
	unsigned long long* Tt_xy;			// Matriz 2D for transmission
	unsigned int* xx;		// Array de posicion en cada step, coordenada x
	unsigned int* zz;		// Array de posicion en cada step, coordenada z
	unsigned long long* banana;			// Matriz 2D para la banana
	unsigned long long* histo_temp;			// Array para el histograma temporal
}MemStruct;

