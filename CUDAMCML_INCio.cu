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

#define NFLOATS 9
#define NINTS 5

#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>

int interpret_arg(int argc, char* argv[], unsigned long long* seed, int* ignoreAdetection)
{

	int unknown_argument;
	for(int i=2;i<argc;i++)
	{
		unknown_argument=1;
		if(!strcmp(argv[i],"-A")) 
		{
			unknown_argument=0;
			*ignoreAdetection=1; //This option is not yet implemented. Therefore, this option has no effect.
			printf("Ignoring A-detection!\n");
		}
		if(!strncmp(argv[i],"-S",2) && sscanf(argv[i],"%*2c %llu",seed))
		{
		unknown_argument=0;
		printf("Seed=%llu\n",*seed);
		}
		if(unknown_argument)
		{
			printf("Unknown argument %s!\n",argv[i]);
			return 1;
		}
	}
	return 0;
}

int Write_Simulation_Results(MemStruct* HostMem, SimulationStruct* sim, clock_t simulation_time)
{
	FILE* pFile_inp;
	char mystring[STR_LEN];

	// Copy stuff from sim->det to make things more readable:
	double dx=(double)sim->det.dx;		// Detection grid resolution, x-direction [cm]
	double dy=(double)sim->det.dy;		// Detection grid resolution, y-direction [cm]
	
	int nx=sim->det.nx;			// Number of grid elements in x-direction
	int ny=sim->det.ny;			// Number of grid elements in y-direction
	

	int x,y;//,z;
	//unsigned int l;
	int i;

	//unsigned long long temp=0;
	double scale1 = (double)0xFFFFFFFF*(double)sim->number_of_photons; // Number of photons (used to normalize)
	double scale2;

	


		
	// Open the input and output files
	pFile_inp = fopen (sim->inp_filename , "r");
	if (pFile_inp == NULL){perror ("Error opening input file");return 0;}

	//pFile_outp = fopen (sim->outp_filename , "w");
	//if (pFile_outp == NULL){perror ("Error opening output file");return 0;}

	// Write other stuff here first!

	char salida_simdata[256];
	strcpy(salida_simdata, sim->outp_filename);
	strcat(salida_simdata,"_simdata.dat");

	char salida_trans[256];
	strcpy(salida_trans, sim->outp_filename);
	strcat(salida_trans,"_trans.dat");

	char salida_ref[256];
	strcpy(salida_ref, sim->outp_filename);
	strcat(salida_ref,"_ref.dat");

	
	char salida_banana[256];
	strcpy(salida_banana, sim->outp_filename);
	strcat(salida_banana,"_banana.dat");

	char salida_t[256];
	strcpy(salida_t, sim->outp_filename);
	strcat(salida_t,"_temp.dat");
	
	FILE* dFile_out;
	dFile_out = fopen (salida_simdata,"w");

	fprintf(dFile_out,"A1 	# Version number of the file format.\n\n");
	fprintf(dFile_out,"####\n");
	fprintf(dFile_out,"# Data categories include: \n");
	fprintf(dFile_out,"# InParm, RAT, \n");
	fprintf(dFile_out,"# A_l, A_z, Rd_r, Rd_a, Tt_r, Tt_a, \n");
	fprintf(dFile_out,"# A_rz, Rd_ra, Tt_ra \n");
	fprintf(dFile_out,"####\n\n");

	// Write simulation time
	fprintf(dFile_out,"# User time: %.2f sec\n\n",(double)simulation_time/CLOCKS_PER_SEC);


	fprintf(dFile_out,"InParam\t\t# Input parameters:\n");
	// Copy the input data from inp_filename
	fseek(pFile_inp, sim->begin, SEEK_SET);
	while(sim->end>ftell(pFile_inp))
	{
		
		fgets(mystring , STR_LEN , pFile_inp);
		fputs(mystring , dFile_out );
	}

	
	fclose(pFile_inp);
	fclose(dFile_out);
	
	FILE* trFile_out;
	trFile_out = fopen (salida_trans,"w");
	
	i=0;
	//fprintf(pFile_outp,"\n\n# T[x][y]. [1/(cm2)].\n# T[0][0], [0][1],..[0][ny-1]\n# T[1][0], [1][1],..[1][ny-1]\n# ...\n# T[nx-1][0], [nx-1][1],..[nx-1][ny-1]\nT\n");
	for(y=0;y<ny;y++)
	{
		for(x=0;x<nx;x++)
		{
			scale2=scale1*dx*dy; // Normalization Constant
			fprintf(trFile_out," %E ",(double)HostMem->Tt_xy[y*nx+x]/scale2);
		}
		fprintf(trFile_out," \n ");
	}

	fclose(trFile_out);

	FILE* reFile_out;
	reFile_out = fopen (salida_ref,"w");

	i=0;
	//fprintf(pFile_outp,"\n\n# R[x][y]. [1/(cm2)].\n# R[0][0], [0][1],..[0][ny-1]\n# R[1][0], [1][1],..[1][ny-1]\n# ...\n# R[nx-1][0], [nx-1][1],..[nx-1][ny-1]\nR\n");
	for(y=0;y<ny;y++)
	{
		for(x=0;x<nx;x++)
		{
			scale2=scale1*dx*dy; // Normalization Constant
			fprintf(reFile_out," %E ",(double)HostMem->Rd_xy[y*nx+x]/scale2);
		}
		fprintf(reFile_out," \n ");
	}

	fclose(reFile_out);

		
	
	FILE* bFile_out;
	bFile_out = fopen (salida_banana,"w");
	int max_z=(int)(sim->esp)*TAM_GRILLA;
	int max_x=(int)2*(sim->esp)*TAM_GRILLA;	
	for(int ix=0;ix<max_x;ix++){                                            /*Reconvierte i y j a x y z*/
        	for(int jz=0;jz<max_z;jz++){
            		fprintf(bFile_out, "%E\t", (double)HostMem->banana[ix*max_z+jz]/scale1);
        	}
        	fprintf(bFile_out,"\n");
    	}
	fclose (bFile_out);
	
	FILE* tFile_out;
	tFile_out = fopen (salida_t,"w");
	for(int itemp=0;itemp<NUM_CAN_TEMP;itemp++){                                            /*Reconvierte i y j a x y z*/
        		fprintf(tFile_out, " %E\t", (double)HostMem->histo_temp[itemp]/scale1);
        	fprintf(tFile_out,"\n");
    	}
	fclose (tFile_out);
	
	return 0;

}


int isnumeric(char a)
{
	if(a>=(char)48 && a<=(char)57) return 1;
	else return 0;
}

int readfloats(int n_floats, float* temp, FILE* pFile)
{
	int ii=0;
	char mystring [200];

	if(n_floats>NFLOATS) return 0; //cannot read more than NFLOATS floats

	while(ii<=0)
	{
		if(feof(pFile)) return 0; //if we reach EOF here something is wrong with the file!
		fgets(mystring , 200 , pFile);
		memset(temp,0,NFLOATS*sizeof(float));
		ii=sscanf(mystring,"%f %f %f %f %f %f %f %f %f",&temp[0],&temp[1],&temp[2],&temp[3],&temp[4],&temp[5],&temp[6],&temp[7],&temp[8]);
		if(ii>n_floats) return 0; //if we read more number than defined something is wrong with the file!
	}
	return 1; // Everyting appears to be ok!
}

int readints(int n_ints, int* temp, FILE* pFile) //replace with template?
{
	int ii=0;
	char mystring[STR_LEN];

	if(n_ints>NINTS) return 0; //cannot read more than NFLOATS floats

	while(ii<=0)
	{
		if(feof(pFile)) return 0; //if we reach EOF here something is wrong with the file!
		fgets(mystring , STR_LEN , pFile);
		memset(temp,0,NINTS*sizeof(int));
		ii=sscanf(mystring,"%d %d %d %d %d",&temp[0],&temp[1],&temp[2],&temp[3],&temp[4]);
		if(ii>n_ints) return 0; //if we read more number than defined something is wrong with the file!
	}
	return 1; // Everyting appears to be ok!
}

int ischar(char a)
{
	if((a>=(char)65 && a<=(char)90)||(a>=(char)97 && a<=(char)122)) return 1;
	else return 0;
}

int read_simulation_data(char* filename, SimulationStruct** simulations, int ignoreAdetection)
{
	int i=0;
	int ii=0;
	int iii=0;
	unsigned long number_of_photons;
	int n_inclusions = 0;
	unsigned int start_weight;
	int n_simulations = 0;
	int n_layers = 0;
	FILE * pFile;
	char mystring [STR_LEN];
	char str[STR_LEN];
	char AorB;
	float dtot=0;


	float ftemp[NFLOATS];//Find a more elegant way to do this...
	int itemp[NINTS];


	pFile = fopen(filename , "r");
	if (pFile == NULL){perror ("Error opening file");return 0;}
	
	// First read the first data line (file version) and ignore
	if(!readfloats(1, ftemp, pFile)){perror ("Error reading file version");return 0;}
	//printf("File version: %f\n",ftemp[0]);

	// Second, read the number of runs
	if(!readints(1, itemp, pFile)){perror ("Error reading number of runs");return 0;}
	n_simulations = itemp[0];
	printf("Number of runs: %d\n",n_simulations);
	
	// Allocate memory for the SimulationStruct array
	*simulations = (SimulationStruct*) malloc(sizeof(SimulationStruct)*n_simulations);
	if(*simulations == NULL){perror("Failed to malloc simulations.\n");return 0;}//{printf("Failed to malloc simulations.\n");return 0;}

	for(i=0;i<n_simulations;i++)
	{
		// Store the input filename
		strcpy((*simulations)[i].inp_filename,filename);
		// Echo the Filename
		//printf("Input filename: %s\n",filename);

		// Store ignoreAdetection data
		(*simulations)[i].ignoreAdetection=ignoreAdetection;

		// Read the output filename and determine ASCII or Binary output
		ii=0;
		while(ii<=0)
		{
			(*simulations)[i].begin=ftell(pFile);
			fgets (mystring , STR_LEN , pFile);
			ii=sscanf(mystring,"%s %c",str,&AorB);
			if(feof(pFile)|| ii>2){perror("Error reading output filename");return 0;}
			if(ii>0)ii=ischar(str[0]);
		}
		// Echo the Filename and AorB
		//printf("Output filename: %s, AorB=%c\n",str,AorB);
		strcpy((*simulations)[i].outp_filename,str);
		(*simulations)[i].AorB=AorB;

		//printf("begin=%d\n",(*simulations)[i].begin);



		// Read the number of photons
		ii=0;
		while(ii<=0)
		{
			fgets(mystring , STR_LEN , pFile);
			number_of_photons=0;
			ii=sscanf(mystring,"%lu",&number_of_photons);
			if(feof(pFile) || ii>1){perror("Error reading number of photons");return 0;} //if we reach EOF or read more number than defined something is wrong with the file!
			//printf("ii=%d temp=%f %f %f %f %f\n",ii,temp[0],temp[1],temp[2],temp[3],temp[4]);
		}
		//printf("Number of photons: %lu\n",number_of_photons);
		(*simulations)[i].number_of_photons=number_of_photons;

		// Read dr and dz (3x float)
		if(!readfloats(2, ftemp, pFile)){perror ("Error reading dr and dz");return 0;}
		//printf("dz=%f, dx=%f, dy=%f\n",ftemp[0],ftemp[1],ftemp[2]);
		//(*simulations)[i].det.dz=ftemp[0];
		(*simulations)[i].det.dx=ftemp[0];
		(*simulations)[i].det.dy=ftemp[1];
		
		// Read No. of dz, dr and da  (3x int)
		if(!readints(2, itemp, pFile)){perror ("Error reading No. of dz, dr and da");return 0;}

		//printf("No. of dz=%d, dx=%d, dy=%d\n",itemp[0],itemp[1],itemp[2]);
		//(*simulations)[i].det.nz=itemp[0];
		(*simulations)[i].det.nx=itemp[0];
		(*simulations)[i].det.ny=itemp[1];

		// Leer separacion fuente-detector	
		//if(!readfloats(1, ftemp, pFile)){perror ("Error leyendo separacion fuente-detector");return 0;}
		//printf("Useparacion fuente-detector=%f\n",ftemp[0]);
		//(*simulations)[i].det.sep=ftemp[0];
		
		// Leer posicion de la fibra (3x float)
		if(!readfloats(3, ftemp, pFile)){perror ("Error leyendo fix and fiy");return 0;}
		(*simulations)[i].det.face=(int)ftemp[0];
		(*simulations)[i].det.fix=ftemp[1];
		(*simulations)[i].det.fiy=ftemp[2];
		
		printf("Detector: cara=%u, x=%f, y=%f\n",(*simulations)[i].det.face,(*simulations)[i].det.fix,(*simulations)[i].det.fiy);

		
		
		// Leer posicion de la fuente (2x float)
		if(!readfloats(2, ftemp, pFile)){perror ("Error leyendo fx and fy");return 0;}
		(*simulations)[i].fx=ftemp[0];
		(*simulations)[i].fy=ftemp[1];
		printf("Fuente: x=%f, y=%f\n",(*simulations)[i].fx,(*simulations)[i].fy);


		// Read No. of layers (1xint)
		if(!readints(1, itemp, pFile)){perror ("Error reading No. of layers");return 0;}
		printf("No. of layers=%d\n",itemp[0]);
		n_layers = itemp[0];
		(*simulations)[i].n_layers = itemp[0];
		printf("No. of layers of %i = %d\n",i,(*simulations)[i].n_layers);

		// Allocate memory for the layers (including one for the upper and one for the lower)
		(*simulations)[i].layers = (LayerStruct*) malloc(sizeof(LayerStruct)*(n_layers+2));
		if((*simulations)[i].layers == NULL){perror("Failed to malloc layers.\n");return 0;}//{printf("Failed to malloc simulations.\n");return 0;}


		// Read upper refractive index (1xfloat)
		if(!readfloats(1, ftemp, pFile)){perror ("Error reading upper refractive index");return 0;}
		printf("Upper refractive index=%f\n",ftemp[0]);
		(*simulations)[i].layers[0].n=ftemp[0];

		dtot=0;
		for(ii=1;ii<=n_layers;ii++)
		{
			// Read Layer data (5x float)
			if(!readfloats(5, ftemp, pFile)){perror ("Error reading layer data");return 0;}
			printf("n=%f, mua=%f, mus=%f, g=%f, d=%f\n",ftemp[0],ftemp[1],ftemp[2],ftemp[3],ftemp[4]);
			(*simulations)[i].layers[ii].n=ftemp[0];
			(*simulations)[i].layers[ii].mua=ftemp[1];
			(*simulations)[i].layers[ii].g=ftemp[3];
			(*simulations)[i].layers[ii].z_min=dtot;
			dtot+=ftemp[4];
			(*simulations)[i].layers[ii].z_max=dtot;
			if(ftemp[2]==0.0f)(*simulations)[i].layers[ii].mutr=FLT_MAX; //Glass layer
			else(*simulations)[i].layers[ii].mutr=1.0f/(ftemp[1]+ftemp[2]);
		}//end ii<n_layers
		
		//Calcular espesor
		printf("Espesor=%f\n",dtot);
		(*simulations)[i].esp=dtot;
		
		// Read lower refractive index (1xfloat)
		if(!readfloats(1, ftemp, pFile)){perror ("Error reading lower refractive index");return 0;}
		printf("Lower refractive index=%f\n",ftemp[0]);
		(*simulations)[i].layers[n_layers+1].n=ftemp[0];            
		
		// Read number of inclusions (1xint)
		if(!readints(1, itemp, pFile)){perror ("Error leyendo el numero de inclusiones");return 0;}
		printf("Number of inclusions=%d\n",itemp[0]);
		n_inclusions = itemp[0];
		(*simulations)[i].n_inclusions = itemp[0];
		
		// Allocate memory for the inclusions 
		(*simulations)[i].inclusion = (IncStruct*) malloc(sizeof(IncStruct)*(n_inclusions));
		if((*simulations)[i].inclusion == NULL){perror("Failed to malloc inclusions.\n");return 0;}                                                                             

		// Read inclusion data (9xfloat)
		for(ii=0;ii<n_inclusions;ii++)	{
			if(!readfloats(9, ftemp, pFile)){perror ("Error leyendo datos de inclusion");return 0;}
			printf("type= %f, x=%f, y=%f, z=%f, r=%f, n=%f, mua=%f, mus=%f, g=%f\n",ftemp[0],ftemp[1],ftemp[2],ftemp[3],ftemp[4],ftemp[5],ftemp[6],ftemp[7],ftemp[8]);
			(*simulations)[i].inclusion[ii].type=(int)ftemp[0];
			(*simulations)[i].inclusion[ii].x=ftemp[1];
			(*simulations)[i].inclusion[ii].y=ftemp[2];
			(*simulations)[i].inclusion[ii].z=ftemp[3];
			(*simulations)[i].inclusion[ii].r=ftemp[4];
			(*simulations)[i].inclusion[ii].n=ftemp[5];
			(*simulations)[i].inclusion[ii].mua=ftemp[6];
			(*simulations)[i].inclusion[ii].g=ftemp[8];
			if(ftemp[7]==0.0f)(*simulations)[i].inclusion[ii].mutr=FLT_MAX; //Inclusion with mus=0
			//else(*simulations)[i].inclusion.mutr=1.0f/(ftemp[4]+ftemp[5]);
			else //Calculates the corrected mus and mutr
			{
				for (iii=1; iii<=(*simulations)[i].n_layers; iii++){
					float z_min_temp = (*simulations)[i].layers[iii].z_min;
					printf ("iii %i , z min %f\n",iii,z_min_temp);
					float z_max_temp = (*simulations)[i].layers[iii].z_max;
					printf ("iii %i , z max %f\n",iii,z_max_temp);
					float z_temp = ftemp[3];
					printf ("iii %i , z %f\n",iii,z_temp);
					if (z_temp>z_min_temp && z_temp<=z_max_temp){
						(*simulations)[i].inclusion[ii].layer = iii;
						float corr=((*simulations)[i].layers[iii].n/ftemp[5])*((*simulations)[i].layers[iii].n/ftemp[5]);
						printf ("layer= %u, nmed: %f, ninc: %f, corr: %f\n",(*simulations)[i].inclusion[ii].layer,(*simulations)[i].layers[iii].n,ftemp[4],corr);
						(*simulations)[i].inclusion[ii].mutr=1.0f/(ftemp[6]+corr*ftemp[7]);
					}
				}
			}	
		
			printf("inclusion %i: type:%u, x=%f, y=%f, z=%f, r=%f, n=%f, mua=%f, mutr=%f, g=%f\n",ii,(*simulations)[i].inclusion[ii].type,(*simulations)[i].inclusion[ii].x,(*simulations)[i].inclusion[ii].y,(*simulations)[i].inclusion[ii].z,(*simulations)[i].inclusion[ii].r,(*simulations)[i].inclusion[ii].n,(*simulations)[i].inclusion[ii].mua,(*simulations)[i].inclusion[ii].mutr,(*simulations)[i].inclusion[ii].g);
		
		}	
		

		(*simulations)[i].end=ftell(pFile);
		//printf("end=%d\n",(*simulations)[i].end);
		
		
		
		//calculate start_weight
		double n1=(*simulations)[i].layers[0].n;
		double n2=(*simulations)[i].layers[1].n;
		double r = (n1-n2)/(n1+n2);
		r = r*r;
		start_weight = (unsigned int)((double)0xffffffff*(1-r));
		//printf("Start weight=%u\n",start_weight);
		(*simulations)[i].start_weight=start_weight;
		

	}//end for i<n_simulations
	return n_simulations;
}
