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
	
vector spin_rotation(vector Heff, vector s, double dt);

void calculate_spin(atom_struct *atom_ptr, double dt){

    #ifdef extfield
    atom_ptr->Heff_H = atom_ptr->Hext;
    #else
    atom_ptr->Heff_H = vec_zero();
    #endif
    #if defined SDHL || defined SLDHL
    atom_ptr->Heff_L = vec_zero();
    #endif

    #ifndef spinlang
    //exact solution; no thermostat; PRL 86, 898 (2001) I. P. Omelyan
    
    #if (defined extvelocity || defined extbrat)
    if (atom_ptr->fext.x != 0.0 || atom_ptr->fext.y != 0.0 || atom_ptr->fext.z != 0.0){    
	    atom_ptr->Heff_H = vec_zero();
	    atom_ptr->s = spin_rotation(atom_ptr->Heff_H, atom_ptr->s, dt);                                                                                                                    
    }
    else {
        atom_ptr->s_curr = atom_ptr->s;
        atom_ptr->s_next = atom_ptr->s_curr;
	    for (int k = 0; k < 4; ++k) {
		    inner_spin(atom_ptr);
		    atom_ptr->s_next = spin_rotation(atom_ptr->Heff_H, atom_ptr->s_curr, dt);
		    atom_ptr->Heff_H = vec_zero();
		}
        atom_ptr->s = atom_ptr->s_next;                           
	}
	#endif
	
	#if not (defined extvelocity || defined extbrat)
    atom_ptr->s_curr = atom_ptr->s;
    atom_ptr->s_next = atom_ptr->s_curr;     
	for (int k = 0; k < 4; ++k) {
		inner_spin(atom_ptr);
		atom_ptr->s_next = spin_rotation(atom_ptr->Heff_H, atom_ptr->s_curr, dt);
		atom_ptr->Heff_H = vec_zero();
	}
	atom_ptr->s = atom_ptr->s_next;
	#endif
    
    
    #endif /*no spin langevin*/

      #ifdef spinlang
      
      #if (defined extvelocity || defined extbrat)                                           
      if (atom_ptr->fext.x != 0.0 || atom_ptr->fext.y != 0.0 || atom_ptr->fext.z != 0.0){    
	      atom_ptr->Heff_H = vec_zero();                                                           
      }                                                                                      
      else {                                                                                 
          inner_spin(atom_ptr); //calculate the effective field of current atom                    
      }	                                                                                     
      #endif                                                                                 
                                                                                             
      #if not (defined extvelocity || defined extbrat)                                       
      inner_spin(atom_ptr); //calculate the effective field of current atom                  
      #endif                                                                                         
      
      double dt_half = dt/2e0;

      #if defined SDH || defined SLDH
      //Pui-Wai Ma and S. L. Dudarev, PHYSICAL REVIEW B 83, 134418 (2011)
      //There are 3 parts. Deterministic -> Stochastic -> Deterministic

      vector s_temp = atom_ptr->s;

      //1st part
      
      vector s_cross_Heff = vec_cross(s_temp, atom_ptr->Heff_H);
      double Heff_H0 = vec_length(atom_ptr->Heff_H);            

      double cos_a  = cos(Heff_H0*(dt_half/hbar));
      double sin_a  = sin(Heff_H0*(dt_half/hbar));
      double exp_b  = exp(-Heff_H0*atom_ptr->s0*gamma_S_H*(dt_half/hbar));
      double exp_2b = exp_b*exp_b;
      double normalized_S_dot_H = 0e0;
      if (atom_ptr->s0 > 0e0 && Heff_H0 > 0e0)
          normalized_S_dot_H = vec_dot(s_temp, atom_ptr->Heff_H)/atom_ptr->s0/Heff_H0;

      double denominator = (1e0 + exp_2b + normalized_S_dot_H*(1e0 - exp_2b))*Heff_H0;
      double factor      = (1e0 - exp_2b + normalized_S_dot_H*(1e0 + exp_2b - 2e0*cos_a*exp_b))*atom_ptr->s0;

      atom_ptr->s = vec_divide( vec_add(vec_add(vec_times(2e0*cos_a*exp_b*Heff_H0, s_temp),
                                                vec_times(2e0*sin_a*exp_b, s_cross_Heff)),
                                        vec_times(factor, atom_ptr->Heff_H)), denominator);

      //2nd part
      #ifdef eltemp
      double random_h = sqrt(2e0*(first_cell_ptr+(atom_ptr->new_cell_index))->Te*gamma_S_H*hbar/dt);
      #else
      double random_h = sqrt(2e0*temperature*gamma_S_H*hbar/dt);
      #endif
      int thread_index = omp_get_thread_num();
      vector dh;
      dh.x = random_h*normal_rand(thread_index);
      dh.y = random_h*normal_rand(thread_index);
      dh.z = random_h*normal_rand(thread_index);

      atom_ptr->s = spin_rotation(dh, atom_ptr->s, dt);

      //3rd part
      
      s_temp = atom_ptr->s;
      s_cross_Heff = vec_cross(s_temp, atom_ptr->Heff_H);

      normalized_S_dot_H = 0e0;
      if (atom_ptr->s0 > 0e0 && Heff_H0 > 0e0)
          normalized_S_dot_H = vec_dot(s_temp,atom_ptr->Heff_H)/atom_ptr->s0/Heff_H0;

      denominator = (1e0 + exp_2b + normalized_S_dot_H*(1e0 - exp_2b))*Heff_H0;
      factor      = (1e0 - exp_2b + normalized_S_dot_H*(1e0 + exp_2b -2e0*cos_a*exp_b))*atom_ptr->s0;

      atom_ptr->s = vec_divide( vec_add(vec_add(vec_times(2e0*cos_a*exp_b*Heff_H0, s_temp),
                                                vec_times(2e0*sin_a*exp_b, s_cross_Heff)),
                                        vec_times(factor, atom_ptr->Heff_H)), denominator);

      #endif

      #if defined SDHL || defined SLDHL
      // In 5 parts.
      //part 1

      atom_ptr->s = spin_rotation(atom_ptr->Heff_H, atom_ptr->s, dt_half);

      //part 2
      double dt_quad = dt/4e0;
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

      double s_sq;
      double s0;

      //RK2
      s_sq = vec_sq(atom_ptr->s);
      atom_ptr->Heff_L = vec_times(-(2e0*A + 4e0*B*s_sq + 6e0*C*pow(s_sq,2) + 8e0*D*pow(s_sq,3)), atom_ptr->s);
      #ifdef SLDHL
        atom_ptr->sum_Jij_sj = 0e0;
        inner_sum_Jij_sj(atom_ptr);
        s0 = sqrt(s_sq);
        atom_ptr->Heff_HC = vec_zero();
        if (s0 > 0e0) atom_ptr->Heff_HC = vec_times(-atom_ptr->sum_Jij_sj/s0, atom_ptr->s);
        atom_ptr->Heff_L = vec_add(atom_ptr->Heff_L, atom_ptr->Heff_HC);
      #endif
      vector s_temp = vec_add(atom_ptr->s, vec_times(gamma_S_HL*dt_quad, vec_add(atom_ptr->Heff_H, atom_ptr->Heff_L)));

      s_sq = vec_sq(s_temp);
      atom_ptr->Heff_L = vec_times(-(2e0*A + 4e0*B*s_sq + 6e0*C*pow(s_sq,2) + 8e0*D*pow(s_sq,3)), s_temp);
      #ifdef SLDHL
        s0 = sqrt(s_sq);
        atom_ptr->Heff_HC = vec_zero();
        if (s0 > 0e0) atom_ptr->Heff_HC = vec_times(-atom_ptr->sum_Jij_sj/s0, s_temp);
        atom_ptr->Heff_L = vec_add(atom_ptr->Heff_L, atom_ptr->Heff_HC);
      #endif
      atom_ptr->s = vec_add(atom_ptr->s, vec_times(gamma_S_HL*dt_half, vec_add(atom_ptr->Heff_H, atom_ptr->Heff_L)));

      //part 3
      #ifdef eltemp
        double random_S = sqrt(2e0*(first_cell_ptr+(atom_ptr->new_cell_index))->Te*gamma_S_HL/dt);
      #else
        double random_S = sqrt(2e0*temperature*gamma_S_HL/dt);
      #endif
      int thread_index = omp_get_thread_num();
      vector dS;
      dS.x = random_S*normal_rand(thread_index);
      dS.y = random_S*normal_rand(thread_index);
      dS.z = random_S*normal_rand(thread_index);
      atom_ptr->s = vec_add(atom_ptr->s, vec_times(dt, dS));

      //part 4; RK2
      s_sq = vec_sq(atom_ptr->s);
      atom_ptr->Heff_L = vec_times(-(2e0*A + 4e0*B*s_sq + 6e0*C*pow(s_sq,2) + 8e0*D*pow(s_sq,3)), atom_ptr->s);
      #ifdef SLDHL
        s0 = sqrt(s_sq);
        atom_ptr->Heff_HC = vec_zero();
        if (s0 > 0e0) atom_ptr->Heff_HC = vec_times(-atom_ptr->sum_Jij_sj/s0, atom_ptr->s);
        atom_ptr->Heff_L = vec_add(atom_ptr->Heff_L, atom_ptr->Heff_HC);
      #endif
      s_temp = vec_add(atom_ptr->s, vec_times(gamma_S_HL*dt_quad, vec_add(atom_ptr->Heff_H, atom_ptr->Heff_L)));

      s_sq = vec_sq(s_temp);
      atom_ptr->Heff_L = vec_times(-(2e0*A + 4e0*B*s_sq + 6e0*C*pow(s_sq,2) + 8e0*D*pow(s_sq,3)), s_temp);
      #ifdef SLDHL
        s0 = sqrt(s_sq);
        atom_ptr->Heff_HC = vec_zero();
        if (s0 > 0e0) atom_ptr->Heff_HC = vec_times(-atom_ptr->sum_Jij_sj/s0, s_temp);
        atom_ptr->Heff_L = vec_add(atom_ptr->Heff_L, atom_ptr->Heff_HC);
      #endif
      atom_ptr->s = vec_add(atom_ptr->s, vec_times(gamma_S_HL*dt_half, vec_add(atom_ptr->Heff_H, atom_ptr->Heff_L)));

      //part 5
      atom_ptr->s = spin_rotation(atom_ptr->Heff_H, atom_ptr->s, dt_half);

      #endif
    #endif /*spinlang*/

    atom_ptr->s0 = vec_length(atom_ptr->s);

}

void inner_spin(atom_struct *atom_ptr){

    struct atom_struct *work_ptr;

    struct cell_struct *ccell_ptr;
    struct cell_struct *wcell_ptr;

    ccell_ptr = first_cell_ptr + atom_ptr->new_cell_index;
    
    vector sum_s = vec_add(atom_ptr->s_curr,atom_ptr->s_next);
    
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
                atom_ptr->Heff_H = vec_add(atom_ptr->Heff_H, vec_times(Jij_rij, work_ptr->s));
                
            }

            // Tally the contribution to the first anisotropy correction from atom i's neighbors
            if (rsq <= rcut_phi_sq && atom_ptr != work_ptr){
	            
	            double rij0 = sqrt(rsq);
	            double C1 = 0.2;
	            double C2 = 0.1;

                double d_ij = (dphi_ij(rij0)/rij0);        	            
                	            	            
	            //First anistropy correction

                atom_ptr->Heff_H = vec_add(atom_ptr->Heff_H, vec_times(C1*d_ij,rij));
                
                //Subtract second anisotropy contribution
                
                double dd_ij = (ddphi_ij(rij0) - d_ij)/rsq;
	            
	            vector dphi_rij_sum_s = vec_times(d_ij,sum_s);
	            
	            double sum_s_dot_rij = vec_dot(sum_s,rij);
	            
	            vector ddphi_rij_sum_s = vec_times(dd_ij*sum_s_dot_rij,rij); 
	            
	            atom_ptr->Heff_H = vec_add(atom_ptr->Heff_H, vec_times(C2,vec_add(dphi_rij_sum_s,ddphi_rij_sum_s)));                 
                      
        	}    
        	
            work_ptr = work_ptr->next_atom_ptr;
        }
    }
}

vector spin_rotation(vector Heff, vector s, double dt){

    vector omega  = vec_divide(Heff, -hbar);
    double omega0 = vec_length(omega);

    if (omega0 > 0e0){
        omega = vec_divide(omega, omega0);
    } else {
        omega = vec_zero();
    }

    double omega_12 = omega.x*omega.y;
    double omega_23 = omega.y*omega.z;
    double omega_13 = omega.x*omega.z;

    double omega1_sq = omega.x*omega.x;
    double omega2_sq = omega.y*omega.y;
    double omega3_sq = omega.z*omega.z;

    double A = sin(omega0*dt);
    double B = 1e0 - cos(omega0*dt);

    vector s_temp;
    s_temp.x = s.x
             + (s.x*B*(-omega2_sq - omega3_sq)
             +  s.y*(B*omega_12 - A*omega.z)
             +  s.z*(A*omega.y + B*omega_13));
    s_temp.y = s.y
             + (s.y*B*(-omega1_sq - omega3_sq)
             +  s.z*(B*omega_23 - A*omega.x)
             +  s.x*(A*omega.z + B*omega_12));
    s_temp.z = s.z
             + (s.z*B*(-omega1_sq - omega2_sq)
             +  s.x*(B*omega_13 - A*omega.y)
             +  s.y*(A*omega.x + B*omega_23));

    return s_temp;
       
}

#ifdef SLDHL
void inner_sum_Jij_sj(atom_struct *atom_ptr){

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

            double rsq = vec_sq(rij);

            if (rsq < rcut_mag_sq && atom_ptr != work_ptr){

                double rij0 = sqrt(rsq);
                double Jij_rij = Jij(rij0);
                double sj = vec_length(work_ptr->s);
                atom_ptr->sum_Jij_sj += Jij_rij*sj;
            }
      
            work_ptr = work_ptr->next_atom_ptr;

        }
    }
}
#endif


#endif
