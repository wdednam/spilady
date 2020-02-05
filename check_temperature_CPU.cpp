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

#ifdef CPU

#include "spilady.h"

#if defined SDH || defined SDHL || defined SLDH || defined SLDHL
void inner_spin_temp_up(atom_struct *atom_ptr);
void inner_spin_temp_dn(atom_struct *atom_ptr);
#endif

void check_temperature_CPU(int current_step){

    #if defined MD || defined SLDH || defined SLDHL || defined SLDNC
    double sum_ke = 0e0;
    #pragma omp parallel for reduction(+:sum_ke)
    for (int i = 0; i < natom; ++i) sum_ke += (first_atom_ptr + i)->ke;
    double Tl = sum_ke*2e0/3e0/natom/boltz;
    #endif

    #if defined SDH || defined SDHL || defined SLDH || defined SLDHL
    double sum_R_up = 0e0;
    double sum_R_dn = 0e0;
      #ifdef eltemp
      //cell_ptr->sum_R_up and cell_ptr->sum_R_dn were calculated in calculate_temperature();
      #pragma omp parallel for reduction(+:sum_R_up,sum_R_dn)
      for (int i = 0; i < ncells; ++i){
          struct cell_struct *cell_ptr;
          cell_ptr = first_cell_ptr + i;
          sum_R_up += cell_ptr->sum_R_up;
          sum_R_dn += cell_ptr->sum_R_dn;
      }
      #else
      #pragma omp parallel for reduction(+:sum_R_up,sum_R_dn)
      for (int i = 0; i < natom; ++i){
          struct atom_struct *atom_ptr;
          atom_ptr = first_atom_ptr + i;
          #ifdef extfield
          atom_ptr->Heff_H = atom_ptr->Hext;
          #else
          atom_ptr->Heff_H = vec_zero();
          atom_ptr->s_cross_Heff = vec_zero();
          atom_ptr->Ts_dn = 0e0;
          #endif
          inner_spin_temp_up(atom_ptr); //calculate the effective field of current atom
          
          //inner_spin(atom_ptr);
                       
          sum_R_up += vec_sq(vec_cross(atom_ptr->s,atom_ptr->s_cross_Heff));
          inner_spin_temp_dn(atom_ptr);
          sum_R_dn += atom_ptr->Ts_dn;
          
          //sum_R_up += vec_sq(vec_cross(atom_ptr->s, atom_ptr->Heff_H));
          //sum_R_dn += vec_dot(atom_ptr->s,atom_ptr->Heff_H);
      }
      #endif

    double Ts_R = sum_R_up/sum_R_dn/boltz;
    #endif

    #if defined SDHL || defined SLDHL
    double sum_L_up  = 0e0;
    double sum_L_dn  = 0e0;
      #ifdef eltemp
      //cell_ptr->sum_L_up and cell_ptr->sum_L_dn were calculated in calculate_temperature();
      #pragma omp parallel for reduction(+:sum_L_up,sum_L_dn)
      for (int i = 0; i < ncells; ++i){
          struct cell_struct *cell_ptr;
          cell_ptr = first_cell_ptr + i;
          sum_L_up += cell_ptr->sum_L_up;
          sum_L_dn += cell_ptr->sum_L_dn;
      }
      #else
      #pragma omp parallel for reduction(+:sum_L_up,sum_L_dn)
      for (int i = 0; i < natom; ++i){
          struct atom_struct *atom_ptr;
          atom_ptr = first_atom_ptr + i;

          #ifdef SLDHL
          double A = LandauA(atom_ptr->rho);
          double B = LandauB(atom_ptr->rho);
          double C = LandauC(atom_ptr->rho);
          double D = LandauD(atom_ptr->rho);
          #endif
          #ifdef SDHL
          double A = LandauA(1);
          double B = LandauB(1);
          double C = LandauC(1);
          double D = LandauD(1);
          #endif


          #ifdef SLDHL
          atom_ptr->sum_Jij_sj = 0e0;
          inner_sum_Jij_sj(atom_ptr);
          sum_L_dn += 2e0*atom_ptr->sum_Jij_sj/vec_length(atom_ptr->s);
          #endif

          double s_sq = vec_sq(atom_ptr->s);

          #if defined SDHL || defined SLDHL
          atom_ptr->Heff_L = vec_zero();
          #endif
          atom_ptr->Heff_L = vec_times(-(2e0*A + 4e0*B*s_sq + 6e0*C*pow(s_sq,2) + 8e0*D*pow(s_sq,3)), atom_ptr->s);
          #ifdef SLDHL
          atom_ptr->Heff_HC = vec_times(-atom_ptr->sum_Jij_sj/vec_length(atom_ptr->s), atom_ptr->s);
          atom_ptr->Heff_L = vec_add(atom_ptr->Heff_L, atom_ptr->Heff_HC);
          #endif
             
          sum_L_up += vec_sq(vec_add(atom_ptr->Heff_H, atom_ptr->Heff_L));
          sum_L_dn += 6e0*A + 20e0*B*s_sq + 42e0*C*pow(s_sq,2) + 72e0*D*pow(s_sq,3);

      }
      #endif
    double Ts_L = sum_L_up/sum_L_dn/boltz;
    #endif

    #ifdef eltemp
    double Te = 0e0;
    #pragma omp parallel for reduction(+:Te)
    for (int i = 0; i < ncells; ++i) Te += (first_cell_ptr + i)->Te;
    Te /= ncells;
    Te /= boltz;
    #endif

    char out_tmp_front[] = "tmp-";
    char out_tmp[256];
    strcpy(out_tmp,out_tmp_front);
    strcat(out_tmp,out_body);
    strcat(out_tmp,".dat");

    ofstream out_file(out_tmp,ios::app);
    out_file << setiosflags(ios::scientific) << setprecision(15);

    out_file << current_step
             << " " << total_time
             #if defined MD || defined SLDH || defined SLDHL || defined SLDNC
             << " " << Tl
             #endif
             #if defined SDH || defined SDHL || defined SLDH || defined SLDHL
             << " " << Ts_R
             #endif
             #if defined SDHL || defined SLDHL
             << " " << Ts_L
             #endif
             #ifdef eltemp
             << " " << Te
             #endif
             << '\n';

     out_file.close();
}

#if defined SDH || defined SDHL || defined SLDH || defined SLDHL
void inner_spin_temp_up(atom_struct *atom_ptr){

    struct atom_struct *work_ptr;

    struct cell_struct *ccell_ptr;
    struct cell_struct *wcell_ptr;

    ccell_ptr = first_cell_ptr + atom_ptr->new_cell_index;
    
    for (int i = 0; i <= 26; ++i){
        if (i == 26)
            wcell_ptr = ccell_ptr;
        else
            wcell_ptr = first_cell_ptr + (ccell_ptr->neigh_cell[i]);

        work_ptr = wcell_ptr->head_ptr;
        while (work_ptr != NULL){

            vector rij = vec_sub(atom_ptr->r, work_ptr->r);

            //find image of j closest to i
            find_image(rij);

            double rsq  = vec_sq(rij);

            if (rsq < rcut_mag_sq && atom_ptr != work_ptr){

                double rij0 = sqrt(rsq);
                double Jij_rij = Jij(rij0);
                //vector si_cross_sj = vec_cross(atom_ptr->s,work_ptr->s);
                atom_ptr->s_cross_Heff = vec_add(atom_ptr->s_cross_Heff, vec_times(Jij_rij,work_ptr->s));
                
            }

            // Tally the contributions to the anisotropy corrections from atom i's neighbors
            if (rsq <= rcut_phi_sq && atom_ptr != work_ptr){
	            
	            double rij0 = sqrt(rsq);
	            double C1 = 0.2;
	            double C2 = 0.1;
	            
	            double d_ij = (dphi_ij(rij0)/rij0);
	            //vector si_cross_rij = vec_cross(atom_ptr->s,rij);
	            
	            //First anistropy correction
	            //vector dphi_rij = vec_times(d_ij,si_cross_rij);
	            
	            atom_ptr->s_cross_Heff = vec_add(atom_ptr->s_cross_Heff, vec_times(C1*d_ij,rij));
	            
	            // Second anisotropy correction
	            
	            double dd_ij = (ddphi_ij(rij0) - d_ij)/rsq;
	            
	            double si_dot_rij = vec_dot(atom_ptr->s,rij);
	         
	            //vector ddphi_rij_si = vec_times(dd_ij*si_dot_rij,rij); 
	            
	            atom_ptr->s_cross_Heff = vec_add(atom_ptr->s_cross_Heff, vec_times(2e0*C2*dd_ij*si_dot_rij,rij));
            
        	}    
        	
            work_ptr = work_ptr->next_atom_ptr;
        }
    }
}

void inner_spin_temp_dn(atom_struct *atom_ptr){

    struct atom_struct *work_ptr;

    struct cell_struct *ccell_ptr;
    struct cell_struct *wcell_ptr;

    ccell_ptr = first_cell_ptr + atom_ptr->new_cell_index;
    
    for (int i = 0; i <= 26; ++i){
        if (i == 26)
            wcell_ptr = ccell_ptr;
        else
            wcell_ptr = first_cell_ptr + (ccell_ptr->neigh_cell[i]);

        work_ptr = wcell_ptr->head_ptr;
        while (work_ptr != NULL){

            vector rij = vec_sub(atom_ptr->r, work_ptr->r);

            //find image of j closest to i
            find_image(rij);

            double rsq  = vec_sq(rij);

            if (rsq < rcut_mag_sq && atom_ptr != work_ptr){

                double rij0 = sqrt(rsq);
                double Jij_rij = Jij(rij0);
                double si_dot_sj = vec_dot(atom_ptr->s, work_ptr->s);
                atom_ptr->Ts_dn +=  2e0*Jij_rij*si_dot_sj;
                
            }

            // Tally the contributions to the anisotropy corrections from atom i's neighbors
            if (rsq <= rcut_phi_sq && atom_ptr != work_ptr){
	            
	            double rij0 = sqrt(rsq);
	            double C1 = 0.2;
	            double C2 = 0.1;
	            
                //First anistropy correction
                double d_ij = (dphi_ij(rij0)/rij0);	            
	            
	            vector dphi_rij = vec_times(d_ij,rij);
	            
	            double dphi_rij_dot_si = vec_dot(dphi_rij,atom_ptr->s);
	            
	            atom_ptr->Ts_dn += 2e0*C1*dphi_rij_dot_si;
	            
	            // Second anisotropy correction
	            
	            double si_dot_rij_sq = vec_dot(atom_ptr->s,rij)*vec_dot(atom_ptr->s,rij);            
	                                                             
	            double si_cross_rij_sq = vec_sq(vec_cross(atom_ptr->s,rij));	  
	            
	            double dd_ij = (ddphi_ij(rij0) - d_ij)/rsq; 
	            
	            atom_ptr->Ts_dn += dd_ij*C2*(2e0*si_dot_rij_sq - si_cross_rij_sq);
            
        	}    
        	
            work_ptr = work_ptr->next_atom_ptr;
        }
    }
}


#endif

void check_temperature(int current_step){
    check_temperature_CPU(current_step);
}

#endif
