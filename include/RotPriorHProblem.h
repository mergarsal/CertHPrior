#pragma once


/* RotPrior manifold related */
#include "RotPriorHTypes.h"
#include "RotPriorHManifold.h"
#include "RotPriorHUtils.h"


#include <Eigen/Dense>



namespace RotPriorH {


struct ProblemCachedMatrices{
EIGEN_MAKE_ALIGNED_OPERATOR_NEW 
        Matrix34 NablaF_Y = Matrix34::Identity();
        Matrix93 N_ = Matrix93::Identity(), QN_ = Matrix93::Identity();
        Matrix93 T_ = Matrix93::Identity(), QT_ = Matrix93::Identity();
        Vector9 Qr_ = Vector9::Zero();
        /// DEFAULT CONSTRUCTOR with default values
        ProblemCachedMatrices(){}; 
        
        ProblemCachedMatrices(  const Matrix34 & Nabla,
                                const Matrix93 & N, 
                                const Matrix93 & T, 
                                const Matrix93 & QN, 
                                const Matrix93 & QT, 
                                const Vector9 & Qr ) :
                                NablaF_Y(Nabla), N_(N), T_(T), QN_(QN), QT_(QT), Qr_(Qr){}
};


class RotPriorHProblem{
public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW 
        RotPriorHProblem(){};  // Default
        
        RotPriorHProblem(const Matrix9 & C);  

        // Constructor using two vectors of 3XN corresponding features
        RotPriorHProblem(const Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                         const Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                         const Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations);

        // todo
        void setNumberPoints(const Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0) 
        {number_points_ = bearingVectors0.cols();}

                      
        ~RotPriorHProblem(){};

         Matrix9 getDataMatrixC(void) {return data_matrix_C_;}
         

         void setPointCorrespondences(const Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                                      const Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                                      const Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations) 
                        {p0_ = bearingVectors0; p1_ = bearingVectors1; w_ = weights_observations; 
                        number_points_ = bearingVectors0.cols();}
   

         void setMatrixPrecon(Matrix3 & matrix_precon) {Matrix_precon_ = matrix_precon;}

         // Pseudo jacobi preconditioner based on the three largest eigenvalues of C
         Matrix3 computePseudoJacobiPrecon(void);


        /// ACCESSORS

        /** Get a const pointer to the SO(3) x S(2) product manifold */
         const RotPriorHManifold& getRotPriorManifold() const {
            return domain_;
          }
          
          
           /// OPTIMIZATION AND GEOMETRY

          /** Given a matrix Y, this function computes and returns F(Y), the value of
   * the objective evaluated at X */
  double evaluate_objective(const Matrix34 &Y) const;

  double evaluate_objective(const Matrix34 &Y, ProblemCachedMatrices & problem_matrices) const;

    /** Given a matrix Y, this function computes and returns nabla F(Y), the
   * *Euclidean* gradient of F at Y. */

  Matrix34 Euclidean_gradient(const Matrix34 &Y, const ProblemCachedMatrices & problem_matrices) const;

  /** Given a matrix Y in the domain D of the SE-Sync optimization problem and
   * the *Euclidean* gradient nabla F(Y) at Y, this function computes and
   * returns the *Riemannian* gradient grad F(Y) of F at Y */

   Matrix34 Riemannian_gradient(const Matrix34 &Y, const Matrix34 &nablaF_Y) const;

   Matrix34 Riemannian_gradient(const Matrix34 &Y, const ProblemCachedMatrices & problem_matrices) const;

  /* Preconditioner */
  Matrix34 precondition(const Matrix34& X, const Matrix34 & Xdot) const;



  /** Given a matrix Y in the domain D of the SE-Sync optimization problem, the
   * *Euclidean* gradient nablaF_Y of F at Y, and a tangent vector dotY in
   * T_D(Y), the tangent space of the domain of the optimization problem at Y,
   * this function computes and returns Hess F(Y)[dotY], the action of the
   * Riemannian Hessian on dotY */

   Matrix34 Riemannian_Hessian_vector_product(const Matrix34 &Y,
                                              const ProblemCachedMatrices & problem_matrices,
                                              const Matrix34 &dotY) const;

 

    /** Given a matrix Y in the domain D of the SE-Sync optimization problem and a
  tangent vector dotY in T_Y(E), the tangent space of Y considered as a generic
  matrix, this function computes and returns the orthogonal projection of dotY
  onto T_D(Y), the tangent space of the domain D at Y*/
  Matrix34 tangent_space_projection(const Matrix34 &Y, const Matrix34 &dotY) const;

  /** Given a matrix Y in the domain D of the SE-Sync optimization problem and a
   * tangent vector dotY in T_D(Y), this function returns the point Yplus in D
   * obtained by retracting along dotY */
  Matrix34 retract(const Matrix34 &Y, const Matrix34 &dotY) const;


  // compute Residuals
  Eigen::Matrix<double, 1, Eigen::Dynamic> computeResiduals(const Matrix34 & X) const;

  // update points with weights
  void updateWeights(const Eigen::Matrix<double, 1, Eigen::Dynamic> & new_weights);  // This function modifies the weights



private:

  size_t number_points_;

  Eigen::Matrix<double, 3, Eigen::Dynamic> p0_, p1_;
  Eigen::Matrix<double, 1, Eigen::Dynamic> w_;


  /** The product manifold of SO(2) X S(2) ~ E(3) that is the domain of our method */
  RotPriorHManifold domain_;

  Matrix9 data_matrix_C_;

  Matrix3 Matrix_precon_;




}; // end of RotPrior problem class




}  // end of RotPrior namespace


