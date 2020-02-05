/********************************************************************************
*
*   Copyright (C) 2015 Culham Centre for Fusion Energy,
*   United Kingdom Atomic Energy Authority, Oxfordshire OX14 3DB, UK
*
*   Licensed under the Apache License, Version 2.0 (the "License");
*   you may not use this file except in compliance with the License.
*   You may obtain a copy of the License at
*
*       http://www.apache.org/licenses/LICENSE-2.0
*
*   Unless required by applicable law or agreed to in writing, software
*   distributed under the License is distributed on an "AS IS" BASIS,
*   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*   See the License for the specific language governing permissions and
*   limitations under the License.
*
********************************************************************************
*
*   Program: SPILADY - A Spin-Lattice Dynamics Simulation Program
*   Version: 1.0
*   Date:    Aug 2015
*   Author:  Pui-Wai (Leo) MA
*   Contact: info@spilady.ccfe.ac.uk
*   Address: Culham Centre for Fusion Energy, OX14 3DB, United Kingdom
*
********************************************************************************/

#if defined CPU

#include "spilady.h"

#if defined extbrat

void external_brat_CPU(int current_step, int &n_cycle){

    static bool infile_extvelocity = 0;
    static bool reverse_flag = 0;
    static bool forward_flag = 0;
    
    if (current_step == -1) n_cycle = 0;

    if (current_step ==  -1){
        ifstream infile("extvelocity_rupt.in");

        if (infile) {
            cout << "Reading external velocities file!!!" << '\n';
            infile_extvelocity = 1;
            
            int temp;
            infile >> natom;
            for (int i = 0; i < natom; ++i){
                struct atom_struct* atom_ptr;
                atom_ptr = first_atom_ptr + i;
                infile >> temp >> atom_ptr->fext.x >> atom_ptr->fext.y >> atom_ptr->fext.z;
            }
        }
	}
        
    if ((current_step + 1)%(5*interval_of_print_out) == 0){
	    cout << "Checking minimum cross-section" << '\n';

		int ner = 0;
	    double alatt = 2.8655;
	    double radius = pow((0.6802*3.0/(8.0*Pi_num)),1.0/3.0)*alatt;
	    double v0 = 4.0/3.0*Pi_num*pow(radius,3);
	    double deltaz = 2.0*radius;
	    double dz = deltaz;
	    double zmaxbrat = -1000.0;
	    double zminbrat = - zmaxbrat;
	    double vlength = 0.0;
	    
	    const double cte = Pi_num/3.0;
	       
		for (int i = 0; i < natom; ++i){
		    struct atom_struct* atom_ptr;
		    atom_ptr = first_atom_ptr + i;
		    if (atom_ptr->r.z < zminbrat) zminbrat = atom_ptr->r.z;
		    if (atom_ptr->r.z > zmaxbrat) zmaxbrat = atom_ptr->r.z;
		}
		
		vlength = zmaxbrat - zminbrat;
		
		double nlayer = vlength/dz + 1.0;
		int nbrat =  static_cast<int>(nlayer);
		double volume[nbrat];
		double lowerz, upperz, ztemp, up, down, height, smin, satinlayer;
		
		lowerz = ztemp = 0.0;
		smin = 1000.0;
		
		for (int i = 0; i < nbrat; ++i){
			volume[i] = 0.0;
			upperz = lowerz + deltaz;
			for (int j = 0; j < natom; ++j){
				struct atom_struct* atom_ptr;
				atom_ptr = first_atom_ptr + j;
				ztemp = atom_ptr->r.z - zminbrat;
				if (ztemp >= lowerz & ztemp < upperz) {
                	up = upperz - ztemp;
                	down = ztemp - lowerz;
                	if (up >= radius & down >= radius) volume[i] = volume[i] + v0;
                	else if (up >= radius & down <= radius) {
	                	height = radius + down;
                    	volume[i] = volume[i] + cte*pow(height,2)*(3*radius - height);
                	}
                	else if (up <= radius & down >= radius) {
	                	height = radius + up;                                       
	                	volume[i] = volume[i] + cte*pow(height,2)*(3*radius - height);
                	}
                	else if (up <= radius & down <= radius) {
	                	ner = 1;
	                	break;
                	}
                	else ner = 2;
        		}
    		}
        	
        	satinlayer = volume[i]/v0;
        	if (satinlayer <= smin) smin = satinlayer;
        	lowerz = lowerz + dz;
    	}
    	cout << "The current minimum cross-section is " << smin << '\n';
    	
	    if (reverse_flag == 1 & forward_flag == 0) {
		    forward_flag = 1;
		    ifstream infile("extvelocity_cont.in");
		     
		   	if (infile) {
			   	cout << "Reading external velocities file!!!" << '\n';
	            	
	           	int temp; 
	            infile >> natom;
	            for (int i = 0; i < natom; ++i){
		            struct atom_struct* atom_ptr;
	               	atom_ptr = first_atom_ptr + i;                	
	               	infile >> temp >> atom_ptr->fext.x >> atom_ptr->fext.y >> atom_ptr->fext.z;		    			    
		     	}
	     	}
	 	}	
		    	
	    if (reverse_flag == 0 & forward_flag == 1) {
		    forward_flag = 0;
		    ifstream infile("extvelocity_rupt.in");
		     
		   	if (infile) {
			   	cout << "Reading external velocities file!!!" << '\n';
	            	
	           	int temp; 
	            infile >> natom;
	            for (int i = 0; i < natom; ++i){
		            struct atom_struct* atom_ptr;
	               	atom_ptr = first_atom_ptr + i;                	
	               	infile >> temp >> atom_ptr->fext.x >> atom_ptr->fext.y >> atom_ptr->fext.z;		    			    
		     	}
	     	}
	 	}			        	
    
    	if (smin < 0.00001 & reverse_flag == 0) {
	    	int n_tmp = n_cycle + 1;
	    	n_cycle = n_tmp;
	    	cout << "This is rupture cycle " << n_cycle << '\n';
	    	cout << "The contact rupture minimum cross-section is " << smin << '\n';   		    	
	   		reverse_flag = 1;
    	}

	    if (smin > 6.0 & reverse_flag == 1) {
		    cout << "The contact formation minimum cross-section is " << smin << '\n';   
	   		reverse_flag = 0;
		}
		    	
	}

    if (infile_extvelocity == 0){

        if (current_step == 0) cout <<  "User defined external velocities apply." << '\n';

        double pushing_force = 8e8;
        double time_const = 4e-11; // = total_time = 1 ps

        #pragma omp parallel for
        for (int i = 0; i < natom; ++i){
            struct atom_struct* atom_ptr;
            atom_ptr = first_atom_ptr + i;
            if (atom_ptr->r.z < 2e0*unit_cell_edge_z || atom_ptr->r.z > 8e0*unit_cell_edge_z){
                atom_ptr->fext.z = pushing_force*time_const *(atom_ptr->r.z - box_length_half.z)/box_length_half.z;
            } else {
                atom_ptr->fext.z = 0e0;
            }
            atom_ptr->fext.x = 0e0;
            atom_ptr->fext.y = 0e0;
        }
    }
}

void external_brat(int current_step, int &n_cycle){
    external_brat_CPU(current_step,n_cycle);
}

#endif
#endif
