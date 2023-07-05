#include "interpolation.hpp"
#include "basisfunctions.hpp"
#include "linalg.hpp"

#include <numeric>
#include <cmath>
#include <exception>

namespace cie
{
namespace splinekernel
{

//using ControlPoints2D = std::array<std::vector<double>, 2>;
//using ControlPointsAndKnotVector = std::pair<ControlPoints2D, std::vector<double>>;

ControlPointsAndKnotVector interpolateWithBSplineCurve( const ControlPoints2D& interpolationPoints, size_t polynomialDegree )
{
    ControlPointsAndKnotVector Q;
    Q.first[0] = interpolationPoints[0];
    Q.first[1] = interpolationPoints[1];
    std::vector<double> parameterPositions = centripetalParameterPositions( interpolationPoints );
    Q.second = knotVectorUsingAveraging( parameterPositions, polynomialDegree );

    
    // Throw exception if number of x-values is not equal number of y-values
      if (Q.first[0].size() != Q.first[1].size()) {
       throw std::invalid_argument("Arrays do not have the same length");
      // throw Q.first[0].size();
    }

    
    // Define the Matrix of Shape function
    
    std :: size_t size = Q.first[0].size();
    std :: cout << "This is the size :" << size << std :: endl; 
    linalg :: Matrix N = linalg :: Matrix(size, size,0.0);

   
    // evaluateBSplineBasis( double t, size_t i, size_t p, const std::vector<double>& knotVector )
    for (std :: size_t i = 0 ; i < size  ; i++){
        for (std :: size_t j = 0 ; j < size; j++){
       
        N(i,j) = evaluateBSplineBasis(parameterPositions[i], j, polynomialDegree, Q.second);
     
        }
    }   
    
    
    for (std :: size_t i = 0 ; i < size ; i++){
    
        for (std :: size_t j = 0 ; j < size ; j++){
        
         std :: cout << N.operator()(i,j) << " ";
        }
        std :: cout <<"One line is over" << std :: endl;
    }
    

     
    // P
    ControlPointsAndKnotVector P;
    
    
    // Solve and aqcuire new points
    
    P.first[0] = linalg:: solve(N,Q.first[0]);  
    P.first[1] = linalg:: solve(N,Q.first[1]);
    //acquire knotVector
    std::vector<double> parameterPositionsP = centripetalParameterPositions( P.first );
   // P.second = knotVectorUsingAveraging( parameterPositionsP, polynomialDegree );
    P.second = Q.second;

    
    return {P.first,P.second};
}

//first step, we need to create t_k = [...] using centripetal technique


std::vector<double> centripetalParameterPositions( const ControlPoints2D& interpolationPoints )
{
    // first acquire number of Control points then calculate the Euclidean distance between every two point
    std :: size_t size = interpolationPoints[0].size();

    std :: array <std::vector <double> , 2> Temp;
    std :: vector <double> Diff;

    for (std::size_t i = 0 ; i < size-1 ; i++)
    {
        Temp[0].push_back (std :: pow((interpolationPoints[0][i+1] - interpolationPoints[0][i]),2));
        Temp[1].push_back(std :: pow((interpolationPoints[1][i+1] - interpolationPoints[1][i]),2));
        Diff.push_back( std :: sqrt(std :: sqrt ( Temp[0][i]+Temp[1][i])));
    }

    //now according to lecture we have calculated the d_k , next step calculate the d
    double d;
    for (std :: size_t i = 0 ; i < size-1 ; i++){
        d = d + Diff[i];
    }

    std :: vector <double> t(size);
    t[0] = 0;
    t[size- 1] = 1;
    for(std :: size_t i = 1 ; i < size-1 ; i++)
    {
        t[i] = t[i-1] + (Diff[i-1])/d;
    }

    // after calculating the d we must calculate t_k which is something like [0 .... 1]

    return t;
}

std::vector<double> knotVectorUsingAveraging( const std::vector<double>& parameterPositions, size_t polynomialDegree )
{
    
    // Throw exception if polynomial degree is too high for given number of points
    if (polynomialDegree > (parameterPositions.size()-1) || polynomialDegree <= 0) {
        throw std::out_of_range("Invalid polynomial Degree: " + std::to_string(polynomialDegree));
    }


    std :: size_t n = parameterPositions.size();                        
    std :: vector <double> U;

    for (std :: size_t i = 0 ; i < polynomialDegree+1; i++){
        U.push_back(0);
    }
    for (std :: size_t i = 0 ; i < (n -(polynomialDegree+1)) ; i++)
    {
        double temp{0};
       
            for(std :: size_t j = 1 ; j < (polynomialDegree+1);j++){
             temp = temp+parameterPositions[i+j];
            }
        U.push_back(temp/polynomialDegree);

    }

    for (std :: size_t i = 0 ; i < polynomialDegree+1; i++){
        U.push_back(1);
    }

    

    return U;
}

} // namespace splinekernel
} // namespace cie
