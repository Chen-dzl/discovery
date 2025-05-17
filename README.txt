#project introduction
  This project is used to solve nonlinear equations using the difference form of Newton's method. The characteristic is the use of variable step size and center difference or backward difference (instead of the common forward difference).

#Environmental dependence

#Catalog Structure Description
  |——README.txt                                               //Help document
  |——main.cpp                                                   //The initial idea, no need
  |——CenterDifferenceNewtonMethod.cpp   //Available code, including general forward difference for testing and comparing performance
  |——performance_test.xlsx                             //Performance Test Files

#instructions
  For backward differencing, choose a larger initial step size as much as possible. The specific reasons can be found in the performance_test.xlsx

#Version content update
## v2.0 :
     1.Using center differential or backward differential with variable step size
     2.Due to the variable step size, the center difference only needs to be used to calculate the function value twice in each iteration, just like forward difference, or only one function value is required for backward differencing at each iteration
    3.According to the testing of the Laval nozzle formula, the forward difference and center difference are stable for most initial values and initial step sizes, and in most cases, the number of times the central differential calculation function value is less than the forward differential calculation. In a few cases, the backward difference becomes very slow, but in most cases, the number of times the backward difference calculates the function value is half that of the forward difference
   4.Improving the problem that the number of iterations of central difference and backward difference is too large in some cases. The solution is to use forward difference when the number of iterations reaches the expected number. The test shows that this method can significantly reduce the number of iterations of backward difference in some cases