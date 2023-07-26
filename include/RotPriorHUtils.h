#pragma once 
 
 
#include "RotPriorHTypes.h"
 
 
namespace RotPriorH{

   Matrix3 computeHfromRtPrior(const Matrix2 & rot, const Vector3 trans, const Vector3 n);
   Matrix3 computeHfromRtPrior(const Matrix34 & rt);
   Matrix3 computeHfromRt(const Matrix3 & R, const Vector3 & t, const Vector3 & n);
   Matrix3 computeHfromRt(const Matrix35& Rt);
   
   void computeHomography(Vector3 & n, double d_plane, 
                               Matrix3 & R1, Matrix3 & R2, 
                               Vector3 & t1, Vector3 & t2, 
                               Matrix3 & H);
                               

   
   Matrix35 convertRt2ExpandRt(const Matrix34 & Rt, double pad_int);
   
   Matrix2 projectRtoPrior(const Matrix3 & R);
   Matrix2 obtainRot(const Matrix2 & R);
   Matrix3 SymProduct(const Matrix3 & Ydot, const Matrix3 & Y, const Matrix3 & nabla_Y);
   void createDataCReduced(const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                                 const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                                 const   Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations, 
                                 Matrix9& C, Matrix3& Qr, Matrix39 & Qtr);
                                 
                                 
   Vector9 vec(Matrix3 & M);
   Eigen::MatrixXd skew(const Eigen::MatrixXd & M);
   
   Matrix9 constructDataMatrix(const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                                     const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                                     const   Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations);
                                     
   void createDataCReduced(const Matrix9& C, Matrix3& Qr, Matrix39 & Qtr);
   
   void computeInitPattern(const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                                 const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                                 const   Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations,                                  
                                 std::vector<Matrix3>* H, std::vector<Matrix2> * rot, std::vector<Vector3>* trans, std::vector<Vector3> * n);
                                 
                                 
   void computeInitPattern2(const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors0,
                                 const   Eigen::Matrix<double, 3, Eigen::Dynamic> & bearingVectors1,
                                 const   Eigen::Matrix<double, 1, Eigen::Dynamic> & weights_observations, 
                                 Matrix3& H, Matrix2 & rot, Vector3& trans, Vector3 & n);
   Matrix3 projectToHomography(const Matrix3 & H_hat, const Matrix2 & rot, Vector3 & trans, Vector3 & n);           
   
   /* From colmap: 
   src/base/homography_matrix
   */
   void PoseFromHomographyMatrix(const Matrix3& H, Matrix3* R, Vector3* t,Vector3* n);                                                                                        
}  // end of namespace RotPrior
