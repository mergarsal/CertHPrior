#pragma once

#include <Eigen/Dense>

#include "RotPriorHTypes.h"

/*Define the namespace*/
namespace RotPriorH{

  class RotPriorHManifold{
  

      public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW 
        /* Default constructor */
        RotPriorHManifold(void){};

        /*Delete each component manifold*/
        ~RotPriorHManifold(void){};
        
        /// GEOMETRY
        /** Given a generic matrix A in R^{3 x 3}, this function computes the
        * projection of A onto R (closest point in the Frobenius norm sense).  */
        Matrix34 project(const Matrix34 &A) const;
        
        
        /** Given an element Y in M and a tangent vector V in T_Y(M), this function
       * computes the retraction along V at Y using the QR-based retraction
       * specified in eq. (4.8) of Absil et al.'s  "Optimization Algorithms on
       * Matrix Manifolds").
       */
      Matrix34 retract(const Matrix34 &Y, const Matrix34 &V) const;
            
      /* Projections for each manifold */
      /// Euclidean Y = X
      Vector3 ProjEuc(const Vector3 &t, const Vector3 &Vt) const;
      /// Sphere
      Vector3 ProjSphere(const Vector3 &t, const Vector3 &Vt) const;
      /// Rotation  
      Matrix2 ProjRotation(const Matrix2 & R, const Matrix2 & VR) const;
                                    
                               
   /** Given an element Y in M and a matrix V in T_X(R^{p x kn}) (that is, a (p
   * x kn)-dimensional matrix V considered as an element of the tangent space to
   * the *entire* ambient Euclidean space at X), this function computes and
   * returns the projection of V onto T_X(M), the tangent space of M at X (cf.
   * eq. (42) in the SE-Sync tech report).*/
  Matrix34 Proj(const Matrix34 &Y, const Matrix34 &V) const;
    
      };
} /*end of RotPrior namespace*/

