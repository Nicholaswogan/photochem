#include<iostream>
#include<vector>
#include "yaml-cpp/yaml.h"



struct Photomechanism {
public:
  std::vector<std::string> species_names;

};




int main(){
  
  
  int n = 10;
  int m = 2;
  
  std::vector<double> a;
  
  a.resize(n*m);
  
  for(int i = 0; i < n; i++){
    for(int j= 0; j < m; j++){
      a[i + j*m] = 2.0*i;
      std::cout << std::to_string(a[i+j*m]);
      std::cout << " ";
    }
    std::cout <<"  "<< std::endl;
  }
  
  std::string b = "hi";
  
  
  Photomechanism photomech;
  photomech.species_names.resize(2);

}