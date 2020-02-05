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

#include "spilady.h"

#if defined SDH || defined SDHL || defined SLDH || defined SLDHL

double phi_ij_gen(double rij){

    double rc = 3.5e0;

    double phi_ij = 0e0;
    if (rij <= rc) phi_ij = pow(1e0-rij/rc,4)*exp(1e0-rij/rc);
    
    return phi_ij;
}

double dphi_ij_gen(double rij){

    double rc = 3.5e0;
    double rc5 = pow(rc,5);

    double dphi_ij = 0e0;
    if (rij <= rc) dphi_ij = -(1e0/rc5)*pow(rij-rc,3)*(rij-5e0*rc)*exp(1e0-rij/rc);
    
    return dphi_ij;
}

double ddphi_ij_gen(double rij){

    double rc = 3.5e0;
    double rc6 = pow(rc,6);
    double rsq = pow(rij,2);
    double rc2 = pow(rc,2);

    double ddphi_ij = 0e0;
    if (rij <= rc) ddphi_ij = (1e0/rc6)*pow(rij-rc,2)*(rsq-10e0*rc*rij+21e0*rc2)*exp(1e0-rij/rc);
    
    return ddphi_ij;
}

double dddphi_ij_gen(double rij){

    double rc = 3.5e0;
    double rc7 = pow(rc,7);
    double rsq = pow(rij,2);
    double rc2 = pow(rc,2); 
    double rcub = pow(rij,3);
    double rc3 = pow(rc,3);
     
    double dddphi_ij = 0e0;
    if (rij <= rc) dddphi_ij = -(1e0/rc7)*(rij-rc)*(rcub-15e0*rc*rsq+63e0*rc2*rij-73e0*rc3)*exp(1e0-rij/rc);
    
    return dddphi_ij;
}
	
#endif
