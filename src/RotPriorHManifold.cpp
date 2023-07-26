#include "RotPriorHManifold.h" 

#include "RotPriorHUtils.h"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>

namespace RotPriorH{


            Matrix34 RotPriorHManifold::project(const Matrix34 & P) const
            {
                 Matrix34 A;
                 
                 A.setZero();
                 // 1. Project normal vector
                 A.block<3, 1>(0, 3) = (P.block<3, 1>(0, 3)).normalized();
                 
                 // 2. Project rotation                 
                 Matrix2 R = P.block<2,2>(0, 0);  

                // fill the output matrix
                A.block<2,2>(0,0) = obtainRot(R); 
                
                // 3. Copy translation
                A.block<3,1>(0,2) = P.block<3,1>(0,2);                

                return A;
            }
            
            
            Vector3 RotPriorHManifold::ProjEuc(const Vector3 &t, const Vector3 &Vt) const
            { return  Vt; }
            
            
            Vector3 RotPriorHManifold::ProjSphere(const Vector3 &t, const Vector3 &Vt) const
            { Vector3 tret = Vt - t.dot(Vt) * t;
            return  tret; }



             /* ROTATION */
            /** 
            Matrix2 RotPriorHManifold::ProjRotation(const Matrix2& R, const Matrix2 & VR) const
            {
                Matrix2 RtVr = R.transpose() * VR;
                Matrix2 Vret = (0.5 * (RtVr - RtVr.transpose()));
                // Matrix2 symmRtVr = 0.5 * (RtVr + RtVr.transpose()); 
                // Matrix2 Vret = VR - R * symmRtVr;
                    return Vret;
            **/
            
            // NO se de donde sale esta
             /*
             Matrix2 RotPriorHManifold::ProjRotation(const Matrix2& R, const Matrix2 & VR) const
            {
                Matrix2 RtVr = R.transpose() * VR;
                Matrix2 Vret = (R * 0.5 * (RtVr - RtVr.transpose()));
                    return Vret;
            }
            */
            
            /** Stiefel **/
            Matrix2 RotPriorHManifold::ProjRotation(const Matrix2& R, const Matrix2 & VR) const
            {
                Matrix2 RtVr = R.transpose() * VR;
                Matrix2 symmRtVr = 0.5 * (RtVr + RtVr.transpose()); 
                Matrix2 Vret = VR - R * symmRtVr;
                    return Vret;
            }

            Matrix34 RotPriorHManifold::retract(const Matrix34 &Y, const Matrix34 &V) const {

              // We use projection-based retraction, as described in "Projection-Like
              // Retractions on Matrix Manifolds" by Absil and Malick
              return project(Y + V); 
              // return project(Proj(Y, V));
        }
        
        

         // We just call here to ProjSphere & ProjRotation2
         Matrix34 RotPriorHManifold::Proj(const Matrix34 &Y, const Matrix34 & Vy) const{
         
            Matrix34 V_tan;
            V_tan.setZero();

            Matrix2 R, Vr;
            Vector3 t, Vt;
            Vector3 n, Vn; 
           
            n = Y.block<3,1>(0,3);
            t = Y.block<3,1>(0,2);
            R = Y.block<2,2>(0,0); 
           
            Vn = Vy.block<3,1>(0,3);
            Vt = Vy.block<3,1>(0,2);
            Vr = Vy.block<2,2>(0,0); 
           
            // for the normal vector 
            V_tan.block<3, 1>(0, 3) = ProjSphere(n, Vn);   
                      
            // for the translation
            V_tan.block<3, 1>(0, 2) = ProjEuc(t, Vt);

            // for the rotation
            V_tan.block<2, 2>(0, 0) = ProjRotation(R, Vr);
            
            
            return V_tan;
         }


} // end of essential namespace
