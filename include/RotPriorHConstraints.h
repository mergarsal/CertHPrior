#pragma once 


#include <Eigen/Core>

#include "RotPriorHTypes.h"



namespace RotPriorH{

       void createHRConstraints(std::vector<Matrix12>& As);
        
       void createRtnConstraints(std::vector<Matrix12>& As);
       
       void createHRtnConstraints(std::vector<Matrix21>& As);

       void createHRtnConstraintsRed(std::vector<Matrix21>& As); 
      
} // end of essential namespace
