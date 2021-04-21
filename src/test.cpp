
#include <iostream>
#include <vector>

struct PhotoPlanet{
public:
  double gravity;
  double surface_pressure;
};

struct PhotoMolecules{
  std::vector<std::string> atoms_names;
  std::vector<std::string> species_names;
  std::vector<int> species_composition;
  std::vector<int> lowerboundcond;
};

class Photochem{  
public:
  PhotoMolecules molecules;
  PhotoPlanet planet;
};

int main(){
  
  PhotoPlanet photoplanet;
  
  photoplanet.gravity = 981.0;
  photoplanet.surface_pressure = 1.0;
  
  std::cout << "hello world" << std::endl;
}