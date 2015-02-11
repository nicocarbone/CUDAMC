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

// forward declaration of the device code
__global__ void MCd(MemStruct);
__device__ float rand_MWC_oc(unsigned long long*,unsigned int*);
__device__ float rand_MWC_co(unsigned long long*,unsigned int*);
__device__ void LaunchPhoton(PhotonStruct*,unsigned long long*, unsigned int*, MemStruct);
__global__ void LaunchPhoton_Global(MemStruct);
__device__ void Spin(PhotonStruct*, float,unsigned long long*,unsigned int*);
__device__ unsigned int Reflect(PhotonStruct*, int, unsigned long long*, unsigned int*);
__device__ unsigned int PhotonSurvive(PhotonStruct*, unsigned long long*, unsigned int*);
__device__ void AtomicAddULL(unsigned long long* address, unsigned int add);

__global__ void MCd(MemStruct DeviceMem)
{
    //Block index
    int bx=blockIdx.x;

    //Thread index
    int tx=threadIdx.x;

	int num_hilo=bx*NUM_THREADS_PER_BLOCK+tx;	
	//printf("%i\n",num_hilo);
	

    //First element processed by the block
    int begin=NUM_THREADS_PER_BLOCK*bx;
	
	float esp=layers_dc[(*n_layers_dc)].z_max;	
	unsigned int max_z = __float2uint_rz(esp*TAM_GRILLA);
	
	unsigned long long int x=DeviceMem.x[begin+tx];//coherent
	unsigned int a=DeviceMem.a[begin+tx];//coherent

	float s;	//step length
	//int index_x,index_y;
	int index;
	int b_index;
	double temp;
	int temp_index;
	unsigned int w;
	
	PhotonStruct p = DeviceMem.p[begin+tx];


	int new_layer;

	//float __fdividef(det_dc[0].dx*det_dc[0].nx,2)=__fdividef(det_dc[0].dx*det_dc[0].nx,2);
	//float (-(__fdividef(det_dc[0].dx*det_dc[0].nx,2)))=-(__fdividef(det_dc[0].dx*det_dc[0].nx,2));
	//float __fdividef(det_dc[0].dy*det_dc[0].ny,2)=__fdividef(det_dc[0].dy*det_dc[0].ny,2);
	//float (-(__fdividef(det_dc[0].dy*det_dc[0].ny,2)))=-(__fdividef(det_dc[0].dy*det_dc[0].ny,2));
	//float max_det_z=__fdividef(det_dc[0].dz*det_dc[0].nz,2);
	//float min_det_z=-(__fdividef(det_dc[0].dz*det_dc[0].nz,2));

	//First, make sure the thread (photon) is active
	unsigned int ii = 0;
	if(!DeviceMem.thread_active[begin+tx]) ii = NUMSTEPS_GPU;

	for(;ii<NUMSTEPS_GPU;ii++) //this is the main while loop
	{
		
		if(p.step<MAX_STEP){
			if(fabsf(p.x)<esp) DeviceMem.xx[p.step*NUM_THREADS+num_hilo]= __float2uint_rz((p.x+esp)*TAM_GRILLA);
				else DeviceMem.xx[p.step*NUM_THREADS+num_hilo] = __float2uint_rz((2*esp)*TAM_GRILLA);
		
			if(p.z<esp) DeviceMem.zz[p.step*NUM_THREADS+num_hilo]= __float2uint_rz((p.z)*TAM_GRILLA);
				else DeviceMem.zz[p.step*NUM_THREADS+num_hilo]= __float2uint_rz(esp*TAM_GRILLA);
		}
	
		if(p.inc<0)
		{// Outside of inclusions		
			if(layers_dc[p.layer].mutr!=FLT_MAX)
				s = -__logf(rand_MWC_oc(&x,&a))*layers_dc[p.layer].mutr;//sample step length [cm] //HERE AN OPEN_OPEN FUNCTION WOULD BE APPRECIATED
			else
				s = GLASS_STEP;//temporary, say the step in glass is 100 cm.
		}
		else
		{// Inside inclusion number inc
			if(inclusion_dc[p.inc].mutr!=FLT_MAX)
				s = -__logf(rand_MWC_oc(&x,&a))*(inclusion_dc[p.inc].mutr);//sample step length [cm] //HERE AN OPEN_OPEN FUNCTION WOULD BE APPRECIATED
			else
				s = GLASS_STEP;//temporary, say the step in glass is 100 cm.
		}
		p.inc = -1; // Set inc to default value, outside inclusions
		
		//Check for layer transitions and in case, calculate s
		new_layer = p.layer;
		if(p.z+s*p.dz<layers_dc[p.layer].z_min){new_layer--; s = __fdividef(layers_dc[p.layer].z_min-p.z,p.dz);} //Check for upwards reflection/transmission & calculate new s
		if(p.z+s*p.dz>layers_dc[p.layer].z_max){new_layer++; s = __fdividef(layers_dc[p.layer].z_max-p.z,p.dz);} //Check for downward reflection/transmission
		
		

		p.t += fabsf(s)/(LIGHTSPEED/layers_dc[p.layer].n);
		p.x += p.dx*s;
		p.y += p.dy*s;
		p.z += p.dz*s;
		if (p.step<(MAX_STEP-1)) p.step ++;
		
		
		
		if(new_layer!=p.layer)
		{
			// set the remaining step length to 0
			s = 0.0f;  
 
			if(Reflect(&p,new_layer,&x,&a)==0u)//Check for reflection
			{ // Photon is transmitted
				if(new_layer == 0)
				{ //Diffuse reflectance
										
					//index = __float2int_rz(acosf(-p.dz)*2.0f*RPI*det_dc[0].na)*det_dc[0].nr+min(__float2int_rz(__fdividef(sqrtf(p.x*p.x+p.y*p.y),det_dc[0].dr)),(int)det_dc[0].nr-1);
					
					if(det_dc[0].face==1){//Detector in entry face
						if(((p.x-det_dc[0].fix)*(p.x-det_dc[0].fix)+(p.y-det_dc[0].fiy)*(p.y-det_dc[0].fiy))<=(RAD_FIB_BAN*RAD_FIB_BAN)){
						
						// Temporal Histogram
						
							temp_index =  __float2int_rz(p.t/TEMP_CAN);
   							if (temp_index >= NUM_CAN_TEMP) temp_index=NUM_CAN_TEMP-1;
							if (temp_index < 0) temp_index=0;
							AtomicAddULL(&DeviceMem.histo_temp[temp_index],p.weight);
						
						// Banana
						
							for (int i=0; i<p.step; i++){
								if((fabsf(DeviceMem.xx[i])<esp)||(DeviceMem.zz[i]<esp)){							
									b_index = DeviceMem.xx[i*NUM_THREADS+num_hilo]*max_z+DeviceMem.zz[i*NUM_THREADS+num_hilo];
									AtomicAddULL(&DeviceMem.banana[b_index], p.weight);
								//printf ("%i, %i, %lld\n",i , b_index, DeviceMem.banana[b_index]);
								}				 
							}
						}
					}
					
					if(p.x<__fdividef(det_dc[0].dx*det_dc[0].nx,2) && p.x>(-(__fdividef(det_dc[0].dx*det_dc[0].nx,2))) && p.y<__fdividef(det_dc[0].dy*det_dc[0].ny,2) && p.y>(-(__fdividef(det_dc[0].dy*det_dc[0].ny,2))))
					{//Inside the detector
						// Calculates the position in the exit matrix			
						//index=(__float2int_rz(fabsf(__fdividef(p.x+__fdividef(det_dc[0].dx*det_dc[0].nx,2),det_dc[0].dx))))*det_dc[0].ny+__float2int_rz(fabsf(__fdividef(p.y+__fdividef(det_dc[0].dy*det_dc[0].ny,2),det_dc[0].dy)));
						index=(__float2int_rn(fabsf(__fdividef(p.y+__fdividef(det_dc[0].dy*det_dc[0].ny,2),det_dc[0].dy))))*det_dc[0].nx+__float2int_rn(fabsf(__fdividef(p.x+__fdividef(det_dc[0].dx*det_dc[0].nx,2),det_dc[0].dx)));
						AtomicAddULL(&DeviceMem.Rd_xy[index], p.weight);
					}
				p.weight = 0; // Set the remaining weight to 0, effectively killing the photon
				}
				if(new_layer > *n_layers_dc)
				{	//Transmitted
					//index = __float2int_rz(acosf(p.dz)*2.0f*RPI*det_dc[0].na)*det_dc[0].nr+min(__float2int_rz(__fdividef(sqrtf(p.x*p.x+p.y*p.y),det_dc[0].dr)),(int)det_dc[0].nr-1);
					if(det_dc[0].face==0){//Detector in exit face
						if(((p.x-det_dc[0].fix)*(p.x-det_dc[0].fix)+(p.y-det_dc[0].fiy)*(p.y-det_dc[0].fiy))<=(RAD_FIB_BAN*RAD_FIB_BAN)){
						
						// Temporal Histogram
						
							temp_index =  __float2int_rz(p.t/TEMP_CAN);
   							if (temp_index >= NUM_CAN_TEMP) temp_index=NUM_CAN_TEMP-1;
							if (temp_index < 0) temp_index=0;
							AtomicAddULL(&DeviceMem.histo_temp[temp_index],p.weight);
						
						// Banana
						
							for (int i=0; i<p.step; i++){
								if((fabsf(DeviceMem.xx[i])<esp)||(DeviceMem.zz[i]<esp)){							
									b_index = DeviceMem.xx[i*NUM_THREADS+num_hilo]*max_z+DeviceMem.zz[i*NUM_THREADS+num_hilo];
									AtomicAddULL(&DeviceMem.banana[b_index], p.weight);
								//printf ("%i, %i, %lld\n",i , b_index, DeviceMem.banana[b_index]);
								}				 
							}
						}
					}
					if(fabsf(p.x)<__fdividef(det_dc[0].dx*det_dc[0].nx,2) && fabsf(p.y)<__fdividef(det_dc[0].dy*det_dc[0].ny,2))
					{//Estoy dentro del detector
						// Calculates the position in the exit matrix				
						//index=(__float2int_rz(fabsf(__fdividef(p.x+__fdividef(det_dc[0].dx*det_dc[0].nx,2),det_dc[0].dx))+0.5))*det_dc[0].ny+__float2int_rz(fabsf(__fdividef(p.y+__fdividef(det_dc[0].dy*det_dc[0].ny,2),det_dc[0].dy))+0.5);
						index=(__float2int_rn(fabsf(__fdividef(p.y+__fdividef(det_dc[0].dy*det_dc[0].ny,2),det_dc[0].dy))))*det_dc[0].nx+__float2int_rn(fabsf(__fdividef(p.x+__fdividef(det_dc[0].dx*det_dc[0].nx,2),det_dc[0].dx)));
						AtomicAddULL(&DeviceMem.Tt_xy[index], p.weight);
					}
					//AtomicAddULL(&DeviceMem.Tt_xy[index], p.weight);
					p.weight = 0; // Set the remaining weight to 0, effectively killing the photon
				}
			}
		}

		w=0;
		
		if(s > 0.0f)
		{
			// Drop weight (apparently only when the photon is scattered)
			//if((p.layer == inclusion_dc[0].layer)&&((((p.x-inclusion_dc[0].x)*(p.x-inclusion_dc[0].x))+((p.y-inclusion_dc[0].y)*(p.y-inclusion_dc[0].y))+((p.z-inclusion_dc[0].z)*(p.z-inclusion_dc[0].z)))<=(inclusion_dc[0].r*inclusion_dc[0].r)))
			
			int inc=0;
			
			while (inc<*n_inclusions_dc && p.inc<0)			
			{
				if( ( (inclusion_dc[inc].type==1) && (p.layer == inclusion_dc[inc].layer) && ((((p.x-inclusion_dc[inc].x)*(p.x-inclusion_dc[inc].x))+((p.z-inclusion_dc[inc].z)*(p.z-inclusion_dc[inc].z)))<=(inclusion_dc[inc].r*inclusion_dc[inc].r)) ) || ( (inclusion_dc[inc].type==0) && (p.layer == inclusion_dc[inc].layer) && ((((p.x-inclusion_dc[inc].x)*(p.x-inclusion_dc[inc].x))+((p.y-inclusion_dc[inc].y)*(p.y-inclusion_dc[inc].y))+((p.z-inclusion_dc[inc].z)*(p.z-inclusion_dc[inc].z)))<=(inclusion_dc[inc].r*inclusion_dc[inc].r)) ) )
				{// Estamos dentro de la inclusion
				//printf("entre \n");
				w = __float2uint_rn(inclusion_dc[inc].mua*inclusion_dc[inc].mutr*__uint2float_rn(p.weight));
				p.weight -= w;//__int_as_float(data.w);
				p.inc = inc;	
				Spin(&p,inclusion_dc[inc].g,&x,&a);
				}
				else{	
					inc++;
				}
			
			
			}	
			
			if (p.inc<0) {
				w = __float2uint_rn(layers_dc[p.layer].mua*layers_dc[p.layer].mutr*__uint2float_rn(p.weight));
				p.weight -= w;//__int_as_float(data.w);
				Spin(&p,layers_dc[p.layer].g,&x,&a);
			}
		
		}
		
		w = w&(*ignoreAdetection_dc); //this will set w to 0 if user has specified to ignore detection of absorbed weight

		//if(w!=0u) AtomicAddULL(&DeviceMem.A_xyz[index], w);


		if(!PhotonSurvive(&p,&x,&a)) // Check if photons survives or not
		{
			if(atomicAdd(DeviceMem.num_terminated_photons,1u) < (*num_photons_dc-NUM_THREADS))
			{	// Ok to launch another photon
				LaunchPhoton(&p,&x,&a, DeviceMem);//Launch a new photon
			}
			else
			{	// No more photons should be launched. 
				DeviceMem.thread_active[begin+tx] = 0u; // Set thread to inactive
				ii = NUMSTEPS_GPU;				// Exit main loop
			}
			
		}
	}//end main for loop!
	__syncthreads();//necessary?

	//save the state of the MC simulation in global memory before exiting
	DeviceMem.p[begin+tx] = p;	//This one is incoherent!!!
	DeviceMem.x[begin+tx] = x; //this one also seems to be coherent
	
	

}//end MCd




__device__ void LaunchPhoton(PhotonStruct* p, unsigned long long* x, unsigned int* a, MemStruct DeviceMem)
{
	// We are currently not using the RNG but might do later
	//float input_fibre_radius = 0.03;//[cm]
	//p->x=input_fibre_radius*sqrtf(rand_MWC_co(x,a));

	p->x  = (float)*fx_dc;
	p->y  = (float)*fy_dc;
	p->z  = 0.0f;
	p->t  = 0.0f;
	p->dx = 0.0f;
	p->dy = 0.0f;
	p->dz = 1.0f;
	p->step= 0;
	/*for (int i =0; i<MAX_STEP; i++){
		DeviceMem.xx[i]=0;
		DeviceMem.zz[i]=0;
	}*/
	p->inc=-1;
	
	p->layer = 1;
	p->weight = *start_weight_dc; //specular reflection!

}

__global__ void LaunchPhoton_Global(MemStruct DeviceMem)//PhotonStruct* pd, unsigned long long* x, unsigned int* a)
{
	int bx=blockIdx.x;
    	int tx=threadIdx.x;	

    	//First element processed by the block
    	int begin=NUM_THREADS_PER_BLOCK*bx;

	PhotonStruct p;
	unsigned long long int x=DeviceMem.x[begin+tx];//coherent

	unsigned int a=DeviceMem.a[begin+tx];//coherent

	LaunchPhoton(&p,&x,&a, DeviceMem);

	//__syncthreads();//necessary?
	DeviceMem.p[begin+tx]=p;//incoherent!?
}


__device__ void Spin(PhotonStruct* p, float g, unsigned long long* x, unsigned int* a)
{
	float cost, sint;	// cosine and sine of the 
						// polar deflection angle theta. 
	float cosp, sinp;	// cosine and sine of the 
						// azimuthal angle psi. 
	float temp;

	float tempdir=p->dx;

	//This is more efficient for g!=0 but of course less efficient for g==0
	temp = __fdividef((1.0f-(g)*(g)),(1.0f-(g)+2.0f*(g)*rand_MWC_co(x,a)));//Should be close close????!!!!!
	cost = __fdividef((1.0f+(g)*(g) - temp*temp),(2.0f*(g)));
	if(g==0.0f)
		cost = 2.0f*rand_MWC_co(x,a) -1.0f;//Should be close close??!!!!!

	sint = sqrtf(1.0f - cost*cost);

	__sincosf(2.0f*PI*rand_MWC_co(x,a),&cosp,&sinp);// spin psi [0-2*PI)
	
	temp = sqrtf(1.0f - p->dz*p->dz);

	if(temp==0.0f) //normal incident.
	{
		p->dx = sint*cosp;
		p->dy = sint*sinp;
		p->dz = copysignf(cost,p->dz*cost);
	}
	else // regular incident.
	{
		p->dx = __fdividef(sint*(p->dx*p->dz*cosp - p->dy*sinp),temp) + p->dx*cost;
		p->dy = __fdividef(sint*(p->dy*p->dz*cosp + tempdir*sinp),temp) + p->dy*cost;
		p->dz = -sint*cosp*temp + p->dz*cost;
	}

	//normalisation seems to be required as we are using floats! Otherwise the small numerical error will accumulate
	temp=rsqrtf(p->dx*p->dx+p->dy*p->dy+p->dz*p->dz);
	p->dx = p->dx*temp;
	p->dy = p->dy*temp;
	p->dz = p->dz*temp;
}// end Spin

			

__device__ unsigned int Reflect(PhotonStruct* p, int new_layer, unsigned long long* x, unsigned int* a)
{
	//Calculates whether the photon is reflected (returns 1) or not (returns 0)
	// Reflect() will also update the current photon layer (after transmission) and photon direction (both transmission and reflection)


	float n1 = layers_dc[p->layer].n;
	float n2 = layers_dc[new_layer].n;
	float r;
	float cos_angle_i = fabsf(p->dz);

	if(n1==n2)//refraction index matching automatic transmission and no direction change
	{	
		p->layer = new_layer;
		return 0u;
	}

	if(n2*n2<n1*n1*(1-cos_angle_i*cos_angle_i))//total internal reflection, no layer change but z-direction mirroring
	{
		p->dz *= -1.0f;
		return 1u; 
	}

	if(cos_angle_i==1.0f)//normal incident
	{		
		r = __fdividef((n1-n2),(n1+n2));
		if(rand_MWC_co(x,a)<=r*r)
		{
			//reflection, no layer change but z-direction mirroring
			p->dz *= -1.0f;
			return 1u;
		}
		else
		{	//transmission, no direction change but layer change
			p->layer = new_layer;
			return 0u;
		}
	}
	
	//gives almost exactly the same results as the old MCML way of doing the calculation but does it slightly faster
	// save a few multiplications, calculate cos_angle_i^2;
	float e = __fdividef(n1*n1,n2*n2)*(1.0f-cos_angle_i*cos_angle_i); //e is the sin square of the transmission angle
	r=2*sqrtf((1.0f-cos_angle_i*cos_angle_i)*(1.0f-e)*e*cos_angle_i*cos_angle_i);//use r as a temporary variable
	e=e+(cos_angle_i*cos_angle_i)*(1.0f-2.0f*e);//Update the value of e
	r = e*__fdividef((1.0f-e-r),((1.0f-e+r)*(e+r)));//Calculate r	

	if(rand_MWC_co(x,a)<=r)
	{ 
		// Reflection, mirror z-direction!
		p->dz *= -1.0f;
		return 1u;
	}
	else
	{	
		// Transmission, update layer and direction
		r = __fdividef(n1,n2);
		e = r*r*(1.0f-cos_angle_i*cos_angle_i); //e is the sin square of the transmission angle
		p->dx *= r;
		p->dy *= r;
		p->dz = copysignf(sqrtf(1-e) ,p->dz);
		p->layer = new_layer;
		return 0u;
	}

}

__device__ unsigned int PhotonSurvive(PhotonStruct* p, unsigned long long* x, unsigned int* a)
{	//Calculate wether the photon survives (returns 1) or dies (returns 0)

	if(p->weight>WEIGHTI) return 1u; // No roulette needed
	if(p->weight==0u) return 0u;	// Photon has exited slab, i.e. kill the photon

	if(rand_MWC_co(x,a)<CHANCE)
	{
		p->weight = __float2uint_rn(__fdividef((float)p->weight,CHANCE));
		return 1u;
	}

	//else
	return 0u;
}

//Device function to add an unsigned integer to an unsigned long long using CUDA Compute Capability 1.1
__device__ void AtomicAddULL(unsigned long long* address, unsigned int add)
{
	if(atomicAdd((unsigned int*)address,add)+add<add)
		atomicAdd(((unsigned int*)address)+1,1u);
}
