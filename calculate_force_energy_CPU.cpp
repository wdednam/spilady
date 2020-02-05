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

void inner_loop(atom_struct *atom_ptr){

    #if defined SDH || defined SDHL || defined SLDH || defined SLDHL || defined SLDNC
    double si_sq = vec_sq(atom_ptr->s);
    double C1 = 0.2;
    double C2 = 0.1;
       
    #endif
    
    struct atom_struct *work_ptr;

    struct cell_struct *ccell_ptr;
    struct cell_struct *wcell_ptr;

    ccell_ptr = first_cell_ptr + atom_ptr->new_cell_index;               

    for (int i = 0; i <= 13; ++i){
        if (i == 13)
            wcell_ptr = ccell_ptr;
        else
            wcell_ptr = first_cell_ptr + (ccell_ptr->neigh_cell[i]);

        work_ptr = wcell_ptr->head_ptr;
        
        #if defined SLDH || defined SLDHL                                  
        vector dphi_rij = vec_zero();                                         
        vector siT_dot_Li = vec_zero();
        #endif
        
        
        while (work_ptr != NULL){

            if (work_ptr == atom_ptr && i == 13) break;

            vector rij = vec_sub(atom_ptr->r, work_ptr->r);
            
            vector sij = vec_sub(atom_ptr->s, work_ptr->s);

            find_image(rij);

            double rsq  = vec_sq(rij);                           

            #if defined MD || defined SLDH || defined SLDHL || defined SLDNC
            double pair_enr = 0e0;
            double dudr = 0e0;
            #endif

            #if defined SDH || defined SDHL || defined SLDH || defined SLDHL
            double dudr_spin = 0e0;
            #endif

            #if defined SDH || defined SDHL || defined SLDH || defined SLDHL || defined SLDNC
            double sj_sq = vec_sq(work_ptr->s);
            #endif
            
            #if defined SLDH || defined SLDHL                       
            vector dudr_spin_aniso = vec_zero();
            #endif                                                              
                       
            if (rsq < rcut_max_sq && atom_ptr != work_ptr){

                double rij0 = sqrt(rsq);

                #ifdef localvol
                if (rij0 < rcut_vol){
                    #pragma omp atomic
                    atom_ptr->sum_rij_m1 += 1e0/rij0;
                    #pragma omp atomic
                    atom_ptr->sum_rij_m2 += 1e0/rsq;
                    #pragma omp atomic
                    work_ptr->sum_rij_m1 += 1e0/rij0;
                    #pragma omp atomic
                    work_ptr->sum_rij_m2 += 1e0/rsq;
                }
                #endif

                #if defined MD || defined SLDH || defined SLDHL || defined SLDNC
                if (rij0 < rcut_pot){

                    double dsmallf_rij = dsmallf(rij0);

                    dudr = (dbigf(atom_ptr->rho) + dbigf(work_ptr->rho))*dsmallf_rij + dpair(rij0);

                    #if defined SLDHL
                      dudr += (dLandauA(atom_ptr->rho)*si_sq
                            + dLandauB(atom_ptr->rho)*pow(si_sq,2)
                            + dLandauC(atom_ptr->rho)*pow(si_sq,3)
                            + dLandauD(atom_ptr->rho)*pow(si_sq,4))*dsmallf_rij;
                      dudr += (dLandauA(work_ptr->rho)*sj_sq
                            + dLandauB(work_ptr->rho)*pow(sj_sq,2)
                            + dLandauC(work_ptr->rho)*pow(sj_sq,3)
                            + dLandauD(work_ptr->rho)*pow(sj_sq,4))*dsmallf_rij;
                    #endif

                    pair_enr = pairij(rij0);
                    #pragma omp atomic
                    atom_ptr->pe += 0.5e0*pair_enr;
                    #pragma omp atomic
                    work_ptr->pe += 0.5e0*pair_enr;
                }
                #endif

                #if defined SDH || defined SDHL || defined SLDH || defined SLDHL
                if (rij0 < rcut_mag){

                    double si_dot_sj   = vec_dot(atom_ptr->s, work_ptr->s);               //Si.Sj
                    double si_times_sj = sqrt(si_sq*sj_sq); //|Si|.|Sj|

                    #if defined SLDH || defined SLDHL
                    dudr_spin = -dJij(rij0)*(si_dot_sj - si_times_sj); // -dJdr_ij(Si dot Sj  - |Si||Sj|);
                    #endif

                    double Jij_half = Jij(rij0)/2e0;
                    double J_times =  Jij_half*si_times_sj;
                    double J_dot   = -Jij_half*si_dot_sj;

                    #pragma omp atomic
                    atom_ptr->me0 += J_times;
                    #pragma omp atomic
                    atom_ptr->me  += J_dot;
                    #pragma omp atomic
                    work_ptr->me0 += J_times;
                    #pragma omp atomic
                    work_ptr->me  += J_dot;
                }
                #endif
            
                #if defined SLDH || defined SLDHL
                if (rij0 < rcut_phi && atom_ptr != work_ptr){
	                
	                double si_dot_rij = vec_dot(atom_ptr->s,rij); 
	                double sj_dot_rij = vec_dot(work_ptr->s,rij); 
	                                                              
	                double sij_dot_rij = vec_dot(sij,rij);        
	                                                   
                    double d_ij = (dphi_ij(rij0)/rij0);
                    double dd_ij = (ddphi_ij(rij0) - d_ij)/rsq;
                    
                    // First anisotropy correction
                    dphi_rij = vec_add(dphi_rij,vec_times(d_ij,rij));
                    
                    vector sij_dot_dKi = vec_add(vec_times(d_ij,sij),vec_times(dd_ij*sij_dot_rij,rij));    
                           
                    // Second anisotropy correction - energy contribution
                                        
                    vector siT_dot_dLi = vec_add(vec_times(d_ij,atom_ptr->s),vec_times(dd_ij*si_dot_rij,rij));
                                                                        
                    siT_dot_Li = vec_add(siT_dot_Li,siT_dot_dLi);              
                    
                    // Second anisotropy correction - force contribution
                    
                    double ddd_ij = (dddphi_ij(rij0) - 3e0*ddphi_ij(rij0)/rij0 + 3e0*d_ij/rij0)/pow(rij0,3);
                    
                    vector siT_dot_dLi_dot_si_I = vec_times(dd_ij*(si_sq + sj_sq),rij);
                    
                    vector siT_dot_dLi_dot_si_II = vec_times(ddd_ij*(pow(si_dot_rij,2)+pow(sj_dot_rij,2)),rij);
                    
                    vector siT_dot_dLi_dot_si_III = vec_times(2e0*dd_ij,vec_add(vec_times(si_dot_rij,atom_ptr->s),vec_times(sj_dot_rij,work_ptr->s)));
                    
                    vector siT_dot_dLi_dot_si = vec_add(siT_dot_dLi_dot_si_I,vec_add(siT_dot_dLi_dot_si_II,siT_dot_dLi_dot_si_III));
                    
                    dudr_spin_aniso = vec_add(dudr_spin_aniso,vec_add(vec_times(C1,sij_dot_dKi),vec_times(C2,siT_dot_dLi_dot_si)));
                                                                         
                }
                #endif                              
                                
                #if defined SLDH || defined SLDHL
                dudr += dudr_spin;               
                #endif                           

                #if defined MD || defined SLDH || defined SLDHL || defined SLDNC
                double force = -dudr/rij0;
                
                vector fij = vec_add(vec_times(force, rij),dudr_spin_aniso);
                //vector fij = vec_times(force, rij);
                #pragma omp atomic
                atom_ptr->f.x += fij.x;
                #pragma omp atomic
                atom_ptr->f.y += fij.y;
                #pragma omp atomic
                atom_ptr->f.z += fij.z;

                #pragma omp atomic
                work_ptr->f.x -= fij.x;
                #pragma omp atomic
                work_ptr->f.y -= fij.y;
                #pragma omp atomic
                work_ptr->f.z -= fij.z;

                #pragma omp atomic
                atom_ptr->stress11 += fij.x*rij.x;
                #pragma omp atomic
                atom_ptr->stress22 += fij.y*rij.y;
                #pragma omp atomic
                atom_ptr->stress33 += fij.z*rij.z;
                #pragma omp atomic
                atom_ptr->stress12 += fij.x*rij.y;
                #pragma omp atomic
                atom_ptr->stress23 += fij.y*rij.z;
                #pragma omp atomic
                atom_ptr->stress31 += fij.z*rij.x;
                #pragma omp atomic
                work_ptr->stress11 += fij.x*rij.x;
                #pragma omp atomic
                work_ptr->stress22 += fij.y*rij.y;
                #pragma omp atomic
                work_ptr->stress33 += fij.z*rij.z;
                #pragma omp atomic
                work_ptr->stress12 += fij.x*rij.y;
                #pragma omp atomic
                work_ptr->stress23 += fij.y*rij.z;
                #pragma omp atomic
                work_ptr->stress31 += fij.z*rij.x;

                atom_ptr->vir += -force*rsq - vec_length(dudr_spin_aniso)*rsq;
                #endif
            }
            work_ptr = work_ptr->next_atom_ptr;
        }
        
        //First Anisotropy energy contribution                      
        double Ki_dot_si = -vec_dot(atom_ptr->s,dphi_rij);          
        //Second Anisotropy energy contribution                     
        double siT_dot_Li_dot_si = -vec_dot(siT_dot_Li,atom_ptr->s);
                                                                    
        #pragma omp atomic                                          
        atom_ptr->me  += C1*Ki_dot_si;                      
        #pragma omp atomic
        atom_ptr->me  += C2*siT_dot_Li_dot_si;
        
    }
        
    #if defined MD || defined SLDH || defined SLDHL
    atom_ptr->pe +=bigf(atom_ptr->rho);
    #endif

    #if defined SDH || defined SDHL || defined SLDH || defined SLDHL
        #ifdef extfield
        atom_ptr->me -= vec_dot(atom_ptr->s, atom_ptr->Hext);
        #endif
        #ifdef SLDHL
        atom_ptr->me += LandauA(atom_ptr->rho)*si_sq
                      + LandauB(atom_ptr->rho)*pow(si_sq,2)
                      + LandauC(atom_ptr->rho)*pow(si_sq,3)
                      + LandauD(atom_ptr->rho)*pow(si_sq,4);
        #endif
        #ifdef SDHL
        atom_ptr->me += LandauA(1)*si_sq
                      + LandauB(1)*pow(si_sq,2)
                      + LandauC(1)*pow(si_sq,3)
                      + LandauD(1)*pow(si_sq,4);
        #endif
    #endif
}

void calculate_force_energy_CPU(){

    #pragma omp parallel for
    for (int i = 0; i < natom; ++i){
        struct atom_struct *atom_ptr;
        atom_ptr = first_atom_ptr + i;
        #if defined MD || defined SLDH || defined SLDHL || defined SLDNC
        atom_ptr->f   = vec_zero();
        atom_ptr->pe  = 0e0;
        atom_ptr->vir = 0e0;
        atom_ptr->stress11 = 0e0;
        atom_ptr->stress22 = 0e0;
        atom_ptr->stress33 = 0e0;
        atom_ptr->stress12 = 0e0;
        atom_ptr->stress23 = 0e0;
        atom_ptr->stress31 = 0e0;
        #endif

        #if defined SDH || defined SDHL || defined SLDH || defined SLDHL
        atom_ptr->me  = 0e0;
        atom_ptr->me0 = 0e0;
        #endif
 
        #ifdef localvol
        atom_ptr->sum_rij_m1 = 0e0; //Sum rij^-1
        atom_ptr->sum_rij_m2 = 0e0; //Sum rij^-2
        #endif
    }


    for (int i = 0 ; i < ngroups ; ++i){
        #pragma omp parallel for
        for (int j = 0 ; j < *(allocate_threads_ptr+i); ++j){
            struct atom_struct *atom_ptr;
            atom_ptr = (*(allocate_cell_ptr_ptr + i*max_no_of_members + j))->head_ptr;
            while(atom_ptr != NULL){
                inner_loop(atom_ptr);   // calculate force and energy; both lattice and spin
                atom_ptr = atom_ptr->next_atom_ptr;
            }
        }
    }
    //#pragma omp parallel for
    //for (int i = 0; i < natom; ++i) inner_loop((first_atom_ptr+i));


    double sum_volume = 0e0;

    #pragma omp parallel for reduction(+:sum_volume)
    for (int i = 0; i < natom; ++i){
        struct atom_struct *atom_ptr;
        atom_ptr = first_atom_ptr + i;
        #ifdef localvol
          double local_radius = 0.5e0*atom_ptr->sum_rij_m1/atom_ptr->sum_rij_m2;
          atom_ptr->local_volume = 4e0*Pi_num/3e0*pow(local_radius, 3e0); //it is only an estimation!!!
          sum_volume += atom_ptr->local_volume;
        #else
          atom_ptr->local_volume = box_volume/natom;
        #endif
    }
    #ifdef localvol
      double volume_factor = box_volume/sum_volume;
      #pragma omp parallel for
      for ( int i = 0; i < natom; ++i) (first_atom_ptr+i)->local_volume *= volume_factor;
    #endif

    #if defined MD || defined SLDH || defined SLDHL || defined SLDNC
      #pragma omp parallel for
      for (int i = 0; i < natom; ++i){
          struct atom_struct *atom_ptr;
          atom_ptr = first_atom_ptr + i;
          atom_ptr->stress11 = (pow(atom_ptr->p.x,2)/atmass + atom_ptr->stress11/2e0)/atom_ptr->local_volume;
          atom_ptr->stress22 = (pow(atom_ptr->p.y,2)/atmass + atom_ptr->stress22/2e0)/atom_ptr->local_volume;
          atom_ptr->stress33 = (pow(atom_ptr->p.z,2)/atmass + atom_ptr->stress33/2e0)/atom_ptr->local_volume;
          atom_ptr->stress12 = ((atom_ptr->p.x*atom_ptr->p.y)/atmass + atom_ptr->stress12/2e0)/atom_ptr->local_volume;
          atom_ptr->stress23 = ((atom_ptr->p.y*atom_ptr->p.z)/atmass + atom_ptr->stress23/2e0)/atom_ptr->local_volume;
          atom_ptr->stress31 = ((atom_ptr->p.z*atom_ptr->p.x)/atmass + atom_ptr->stress31/2e0)/atom_ptr->local_volume;
      }
      virial = 0e0;
      #pragma omp parallel for reduction(+:virial)
      for (int i = 0; i < natom; ++i) virial += (first_atom_ptr+i)->vir;
    #endif
}


void calculate_force_energy(){
    calculate_force_energy_CPU();
}

#endif

