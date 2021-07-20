#include <iostream>
#include <vector>

struct ThermoData
{
  int type;
  double* p1;
  double* p2;
};

struct ShomateData : public ThermoData
{
  ShomateData(double P1[7], double P2[7]){
    type = 1;
    p1 = new double[7];
    p2 = new double[7];
    for (int i = 0; i < 7; i++){
      p1[i] = P1[i];
      p2[i] = P2[i];
    }
  }  
  ~ShomateData(){
    delete[] p1;
    delete[] p2;
  }
};

struct NASA7Data : public ThermoData
{
  double* p3;
  
  NASA7Data(double P1[7], double P2[7]){
    type = 2;
    p1 = new double[7];
    p2 = new double[7];
    for (int i = 0; i < 7; i++){
      p1[i] = P1[i];
      p2[i] = P2[i];
    }
  }
  ~NASA7Data(){
    delete[] p1;
    delete[] p2;
  }
};


struct NASA9Data : public ThermoData
{
  double* p3;
  
  NASA9Data(double P1[9], double P2[9]){
    type = 3;
    p1 = new double[9];
    p2 = new double[9];
    for (int i = 0; i < 9; i++){
      p1[i] = P1[i];
      p2[i] = P2[i];
    }
  }
  ~NASA9Data(){
    delete[] p1;
    delete[] p2;
  }
};

int main()
{
  int n = 1000;
  std::vector<ThermoData> vec(n);
  
  double P1[7];
  double P2[7];
  for (int i = 0; i < 7; i++){
    P1[i] = 1.0;
    P2[i] = 2.0;
  }
  
  double P1_N[9];
  double P2_N[9];
  for (int i = 0; i < 9; i++){
    P1_N[i] = 1.0;
    P2_N[i] = 2.0;
  }
  ShomateData shom(P1,P2);
  NASA7Data nas7(P1,P2);
  NASA9Data nas9(P1_N,P2_N);

  for (int i = 0; i < 300; i++)
    vec[i] = shom;
  for (int i = 300; i < 700; i++)
    vec[i] = nas7;
  for (int i = 700; i < n; i++)
    vec[i] = nas9;
    
  for (int i = 0; i < n; i++){
    int& type = vec[i].type;
    if (type == 1){
      for (int j = 0; j < 7; j++)
        std::cout << vec[i].p1[j] << " ";
      std::cout << "\n";
    }
    else if (type == 2){
      for (int j = 0; j < 7; j++)
        std::cout << vec[i].p2[j] << " ";
      std::cout << "\n";
    }
    else if (type == 3){
      for (int j = 0; j < 9; j++)
        std::cout << vec[i].p1[j] << " ";
      std::cout << "\n";
    }
  }
}


