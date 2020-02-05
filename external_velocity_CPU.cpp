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

#if defined extvelocity

void external_velocity_CPU(int current_step){

    static bool infile_extvelocity = 0;

    if (current_step ==  -1){
        ifstream infile("extvelocity.in");

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

void external_velocity(int current_step){
    external_velocity_CPU(current_step);
}

#endif
#endif
