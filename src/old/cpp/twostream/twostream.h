#pragma once
#include <vector>

void solve_tridiag(std::vector<double>& a, std::vector<double>& b, 
                   std::vector<double>& c, std::vector<double>& d, int n);

void two_stream(std::vector<double>& tau, std::vector<double>& w0, 
              double& u0, double& Rsfc, double& surface_radiance, 
              std::vector<double>& amean);