#include <chrono>
#include <iostream>
#include <vector>
#include "twostream.h"

int main()
{
  int nz = 200;
  std::vector<double> tau(nz,0.1);
  std::vector<double> w0(nz,0.9999);
  double u0 = 0.7;
  double Rsfc = 0.25; 
  std::vector<double> amean(nz);
  double surface_radiance;
  
  auto t1 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i<10000; i++)
    two_stream(tau, w0, u0, Rsfc,surface_radiance, amean);
  auto t2 = std::chrono::high_resolution_clock::now();
  
  std::chrono::duration<double, std::milli> ms_double = t2 - t1;
  std::cout << ms_double.count()/1000.0 << " s\n";

}