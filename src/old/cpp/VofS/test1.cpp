#include <iostream>
#include <vector>

struct ThermoData
{
  int type;
  int npoly;
  double* temps = nullptr;
  double* poly = nullptr;
  ThermoData()
  {
  }
  ThermoData(double* Poly, double* Temps, int NPoly, int Type){
    type = Type;
    npoly = NPoly;
    int len;
    if (type == 1){
      len = 7;
    }
    else if (type == 2){
      len = 7;
    }
    else if (type == 3){
      len = 9;
    }
    temps = new double[npoly+1];
    poly = new double[npoly*len];
    for (int i = 0; i < npoly*len; i++){
      poly[i] = Poly[i];
    }
    for (int i = 0; i < npoly+1; i++){
      temps[i] = Temps[i];
    }
  }  
  ~ThermoData(){
    if (poly){
      delete[] poly;
      delete[] temps;
    }
  }
};

void GibbsShomate(double* poly, double& gibbs){
  gibbs = poly[0]*poly[1]*poly[2]*poly[3]+poly[4]/poly[5]*poly[6];
};

void GibbsNASA9(double* poly, double& gibbs){
  gibbs = poly[0]*poly[1]*poly[2]*poly[3]+poly[4]/poly[5]*poly[6]/poly[7]/poly[8];
};

int main()
{
  int n = 1000;
  std::vector<ThermoData> vec;
  vec.reserve(n);
  
  double poly1[7*2];
  double temps1[3];
  temps1[0] = 200.0;
  temps1[1] = 1000.0;
  temps1[2] = 2000.0;
  for (int i = 0; i < 7; i++){
    poly1[i] = 1.0;
  }
  for (int i = 7; i < 7*2; i++){
    poly1[i] = 2.0;
  }
  
  double poly3[9*2];
  for (int i = 0; i < 9; i++){
    poly3[i] = 1.0;
  }
  for (int i = 9; i < 9*2; i++){
    poly3[i] = 3.0;
  }
  
  for (int i = 0; i < 500; i++){
    vec.emplace_back(poly1, temps1, 2, 1);
  }
  for (int i = 500; i < n; i++){
    vec.emplace_back(poly3, temps1, 2, 3);
  }
  
  std::vector<double> gibbs(n);
  for (int i = 0; i < n; i++){
    if (vec[i].type == 1){
      
    }
    else if (vec[i].type == 2){
      GibbsNASA9(vec[i].poly, gibbs[i]);
    }
    
  }
  
  
  std::cout << vec[1].temps[1] << std::endl;

}


