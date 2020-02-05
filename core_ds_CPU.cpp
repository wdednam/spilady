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

#if (defined SDH || defined SDHL || defined SLDH || defined SLDHL || defined SLDNC) && defined CPU

#include "spilady.h"

void core_ds_CPU(double dt){	

    double dt_half = dt/2e0;
    #ifdef magmom
    #pragma omp parallel for
    for(int i = 0; i < natom ; ++i) {
        struct atom_struct *atom_ptr;
        atom_ptr = first_atom_ptr + i;
        atom_ptr->s  = vec_divide(atom_ptr->m,-el_g);
        atom_ptr->s0 = vec_length(atom_ptr->s);
    }
    #endif

    double dt_k = dt_half/5e0;
    for (int k = 0 ; k < 5 ; ++k){
        //Suzuki-Trotter decompsition. Forward and backward. See Eq. (12) of P-W Ma and C. Woo, Phys. Rev. E, 79, 046703, (2009)
        for (int i = 0 ; i < ngroups ; ++i){
            #pragma omp parallel for
            for (int j = 0 ; j < *(allocate_threads_ptr+i); ++j){
                struct atom_struct *atom_ptr;
                atom_ptr = (*(allocate_cell_ptr_ptr + i*max_no_of_members + j))->head_ptr;
                while(atom_ptr != NULL){  
                	calculate_spin(atom_ptr, dt_k);
                    atom_ptr = atom_ptr->next_atom_ptr;
                }
            }
        }
        
        for (int i = ngroups - 1 ; i >=0 ; --i){
            #pragma omp parallel for
            for (int j = 0 ; j < *(allocate_threads_ptr+i); ++j){
                struct atom_struct *atom_ptr;
                atom_ptr = (*(allocate_cell_ptr_ptr + i*max_no_of_members + j))->tail_ptr;
                while(atom_ptr != NULL){      
                	calculate_spin(atom_ptr, dt_k);                               
                    atom_ptr = atom_ptr->prev_atom_ptr;                                                                                                 
                }
            }
        }
 }   
        
        #ifdef magmom
        #pragma omp parallel for
        for(int i = 0; i < natom ; ++i) {
            struct atom_struct *atom_ptr;
            atom_ptr = first_atom_ptr + i;
            atom_ptr->m  = vec_times(-el_g, atom_ptr->s);
            atom_ptr->m0 = vec_length(atom_ptr->m);
        }
        #endif
    }

void core_ds(double dt){
    core_ds_CPU(dt);
}


#endif
