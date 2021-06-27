


struct ThermoData
{
  int type;
  double* temps = nullptr;
  int npoly;
  double* p1 = nullptr;
  double* p2 = nullptr;
};

struct ShomateData : public ThermoData
{
  ShomateData(double P1[7], double Temps[2]){
    type = 1;
    npoly = 1;
    p1 = new double[7];
    temps = new double[2];
    for (int i = 0; i < 7; i++){
      p1[i] = P1[i];
    }
    for (int i = 0; i < 2; i++){
      temps[i] = Temps[i];
    }
  }  
  
  ShomateData(double P1[7], double P2[7], double Temps[3]){
    type = 1;
    npoly = 2;
    p1 = new double[7];
    p2 = new double[7];
    temps = new double[3];
    for (int i = 0; i < 7; i++){
      p1[i] = P1[i];
      p2[i] = P2[i];
    }
    for (int i = 0; i < 3; i++){
      temps[i] = Temps[i];
    }
  }  
  
  ~ShomateData(){
    if (p1 != nullptr){
      delete[] p1;
    }
    if (p2 != nullptr){
      delete[] p2;
    }
  }
};