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

int CopyDeviceToHostMem(MemStruct* HostMem, MemStruct* DeviceMem, SimulationStruct* sim)
{ //Copy data from Device to Host memory

	int xy_size = sim->det.nx*sim->det.ny;
	int max_x=(int)(2*(sim->esp)*TAM_GRILLA);
	int max_z=(int)((sim->esp)*TAM_GRILLA);
	int banana_size = max_x*max_z+max_z;
	//int xyz_size = sim->det.nx*sim->det.ny*sim->det.nz;

	//Copy Rd_xy, Tt_xy and A_xyz
	CUDA_SAFE_CALL( cudaMemcpy(HostMem->Rd_xy,DeviceMem->Rd_xy,xy_size*sizeof(unsigned long long),cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaMemcpy(HostMem->Tt_xy,DeviceMem->Tt_xy,xy_size*sizeof(unsigned long long),cudaMemcpyDeviceToHost) );
	

	CUDA_SAFE_CALL( cudaMemcpy(HostMem->banana,DeviceMem->banana,banana_size*sizeof(unsigned long long),cudaMemcpyDeviceToHost) );	
	CUDA_SAFE_CALL( cudaMemcpy(HostMem->histo_temp,DeviceMem->histo_temp,NUM_CAN_TEMP*sizeof(unsigned long long),cudaMemcpyDeviceToHost) );	
	
	//Also copy the state of the RNG's
	CUDA_SAFE_CALL( cudaMemcpy(HostMem->x,DeviceMem->x,NUM_THREADS*sizeof(unsigned long long),cudaMemcpyDeviceToHost) );

	return 0;
}
int InitDCMem(SimulationStruct* sim)
{
	unsigned int temp=0xFFFFFFFF;
	// Copy ignoreAdetection to constant device memory
	if(sim->ignoreAdetection) temp=0;
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(ignoreAdetection_dc,&temp,sizeof(unsigned int)) );

	// Copy det-data to constant device memory
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(det_dc,&(sim->det),sizeof(DetStruct)) );
	
	// Copy inclusion data to constant device memory
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(inclusion_dc,sim->inclusion,(sim->n_inclusions)*sizeof(IncStruct)) );
	
	// Copy n_inclusions_dc to constant device memory
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(n_inclusions_dc,&(sim->n_inclusions),sizeof(unsigned int)));
	
	// Copy num_photons_dc to constant device memory
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(n_layers_dc,&(sim->n_layers),sizeof(unsigned int)));

	// Copy start_weight_dc to constant device memory
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(start_weight_dc,&(sim->start_weight),sizeof(unsigned int)));

	// Copy layer data to constant device memory
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(layers_dc,sim->layers,(sim->n_layers+2)*sizeof(LayerStruct)) );

	// Copy num_photons_dc to constant device memory
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(num_photons_dc,&(sim->number_of_photons),sizeof(unsigned int)));
	
	// Copy posicion de la fuente en x a la memoria constante del dispositivo
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(fx_dc,&(sim->fx),sizeof(float)));
	
	// Copy posicion de la fuente en y a la memoria constante del dispositivo
	CUDA_SAFE_CALL( cudaMemcpyToSymbol(fy_dc,&(sim->fy),sizeof(float)));

	return 0;
	
}

int InitMemStructs(MemStruct* HostMem, MemStruct* DeviceMem, SimulationStruct* sim)
{
	int xy_size;//,xyz_size;

	xy_size = sim->det.nx*sim->det.ny;

	int max_x=(int)(2*(sim->esp)*TAM_GRILLA);
	int max_z=(int)((sim->esp)*TAM_GRILLA);
	int banana_size = max_x*max_z+max_z;
	int xxzz_size = (MAX_STEP*NUM_THREADS+NUM_THREADS)*sizeof(unsigned int);
	//xyz_size = sim->det.nx*sim->det.ny*sim->det.nz;

	
	// Allocate p on the device!!
	CUDA_SAFE_CALL( cudaMalloc((void**)&DeviceMem->p,NUM_THREADS*sizeof(PhotonStruct)) );
	
	// Reservar memoria para xx y zz en GPU	
		
	CUDA_SAFE_CALL( cudaMalloc((void**)&DeviceMem->xx,xxzz_size ));
	CUDA_SAFE_CALL( cudaMemset(DeviceMem->xx,0,xxzz_size ));
	CUDA_SAFE_CALL( cudaMalloc((void**)&DeviceMem->zz, xxzz_size ));
	CUDA_SAFE_CALL( cudaMemset(DeviceMem->zz,0, xxzz_size ));


	// Allocate Rd_xy on CPU and GPU
	HostMem->Rd_xy = (unsigned long long*) malloc(xy_size*sizeof(unsigned long long));
	if(HostMem->Rd_xy==NULL){printf("Error allocating HostMem->Rd_xy"); exit (1);}
	CUDA_SAFE_CALL( cudaMalloc((void**)&DeviceMem->Rd_xy,xy_size*sizeof(unsigned long long)) );
	CUDA_SAFE_CALL( cudaMemset(DeviceMem->Rd_xy,0,xy_size*sizeof(unsigned long long)) );

	// Allocate Tt_xy on CPU and GPU
	HostMem->Tt_xy = (unsigned long long*) malloc(xy_size*sizeof(unsigned long long));
	if(HostMem->Tt_xy==NULL){printf("Error allocating HostMem->Tt_xy"); exit (1);}
	CUDA_SAFE_CALL( cudaMalloc((void**)&DeviceMem->Tt_xy,xy_size*sizeof(unsigned long long)) );
	CUDA_SAFE_CALL( cudaMemset(DeviceMem->Tt_xy,0,xy_size*sizeof(unsigned long long)) );

	// Allocate banana on CPU and GPU
	HostMem->banana = (unsigned long long*) malloc(banana_size*sizeof(unsigned long long));
	if(HostMem->banana==NULL){printf("Error allocating HostMem->banana"); exit (1);}
	CUDA_SAFE_CALL( cudaMalloc((void**)&DeviceMem->banana,banana_size*sizeof(unsigned long long)) );
	CUDA_SAFE_CALL( cudaMemset(DeviceMem->banana,0,banana_size*sizeof(unsigned long long)) );
	
	// Allocate histo_temp on CPU and GPU
	HostMem->histo_temp = (unsigned long long*) malloc(NUM_CAN_TEMP*sizeof(unsigned long long));
	if(HostMem->histo_temp==NULL){printf("Error allocating HostMem->histo_temp"); exit (1);}
	CUDA_SAFE_CALL( cudaMalloc((void**)&DeviceMem->histo_temp,NUM_CAN_TEMP*sizeof(unsigned long long)) );
	CUDA_SAFE_CALL( cudaMemset(DeviceMem->histo_temp,0,NUM_CAN_TEMP*sizeof(unsigned long long)) );

	// Allocate x and a on the device (For MWC RNG)
    	CUDA_SAFE_CALL(cudaMalloc((void**)&DeviceMem->x,NUM_THREADS*sizeof(unsigned long long)));
    	CUDA_SAFE_CALL(cudaMemcpy(DeviceMem->x,HostMem->x,NUM_THREADS*sizeof(unsigned long long),cudaMemcpyHostToDevice));
	
    	CUDA_SAFE_CALL(cudaMalloc((void**)&DeviceMem->a,NUM_THREADS*sizeof(unsigned int)));
    	CUDA_SAFE_CALL(cudaMemcpy(DeviceMem->a,HostMem->a,NUM_THREADS*sizeof(unsigned int),cudaMemcpyHostToDevice));


	// Allocate thread_active on the device and host
	HostMem->thread_active = (unsigned int*) malloc(NUM_THREADS*sizeof(unsigned int));
	if(HostMem->thread_active==NULL){printf("Error allocating HostMem->thread_active"); exit (1);}
	for(int i=0;i<NUM_THREADS;i++)HostMem->thread_active[i]=1u;

	CUDA_SAFE_CALL( cudaMalloc((void**)&DeviceMem->thread_active,NUM_THREADS*sizeof(unsigned int)) );
	CUDA_SAFE_CALL( cudaMemcpy(DeviceMem->thread_active,HostMem->thread_active,NUM_THREADS*sizeof(unsigned int),cudaMemcpyHostToDevice));


	//Allocate num_launched_photons on the device and host
	HostMem->num_terminated_photons = (unsigned int*) malloc(sizeof(unsigned int));
	if(HostMem->num_terminated_photons==NULL){printf("Error allocating HostMem->num_terminated_photons"); exit (1);}
	*HostMem->num_terminated_photons=0;

	CUDA_SAFE_CALL( cudaMalloc((void**)&DeviceMem->num_terminated_photons,sizeof(unsigned int)) );
	CUDA_SAFE_CALL( cudaMemcpy(DeviceMem->num_terminated_photons,HostMem->num_terminated_photons,sizeof(unsigned int),cudaMemcpyHostToDevice));

	return 1;
}

void FreeMemStructs(MemStruct* HostMem, MemStruct* DeviceMem)
{
	//free(HostMem->A_xyz);
	free(HostMem->Rd_xy);
	free(HostMem->Tt_xy);
	free(HostMem->banana);
	free(HostMem->histo_temp);
	free(HostMem->thread_active);
	free(HostMem->num_terminated_photons);
	
	//cudaFree(DeviceMem->A_xyz);
	cudaFree(DeviceMem->Rd_xy);
	cudaFree(DeviceMem->Tt_xy);
	cudaFree(DeviceMem->xx);
	cudaFree(DeviceMem->zz);
	cudaFree(DeviceMem->banana);    	
	cudaFree(DeviceMem->histo_temp);    		
	cudaFree(DeviceMem->x);
    	cudaFree(DeviceMem->a);
	cudaFree(DeviceMem->thread_active);
	cudaFree(DeviceMem->num_terminated_photons);

}

void FreeSimulationStruct(SimulationStruct* sim, int n_simulations)
{
	for(int i=0;i<n_simulations;i++)free(sim[i].layers);
	free(sim);
}

