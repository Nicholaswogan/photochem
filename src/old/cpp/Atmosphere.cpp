#include <vector>

struct ThermoData
{
  int type;
  double* temps;
  double* p1;
  double* p2;
};

struct Reaction
{
  int rxtype;
  int nreactants;
  int nproducts;
  int* reactants;
  int* products;
};

struct ElementaryReaction : public Reaction
{
  double* rateparams;
};

struct ThreeBodyReaction : public Reaction
{
  double* rateparams;
  int num_efficient;
  double* efficiencies;
  int eff_sp_inds;
  double default_eff;
};

struct FalloffReaction : public Reaction
{
  double* rateparams;
  int num_efficient;
  double* efficiencies;
  int eff_sp_inds;
  double default_eff;
  int falloff_type;
};

class Atmosphere
{
// species and atoms
private:
  std::vector<ThermoData> thermodata;
public:
  int nq;
  int nsl;
  int nsp;
  int natoms;
  std::vector<std::string> atoms_names;
  std::vector<double> atoms_mass;
  std::vector<std::string> species_names;
  std::vector<double> species_mass;
  
// reactions
private:
  bool reverse;
  int max_num_reactants;
  int max_num_products;
  std::vector<int> nump;
  std::vector<int> numl;
  std::vector<int> iprod;
  std::vector<int> iloss;
public:
  int nrF;
  int nrR;
  int nrT;
  std::vector<Reaction> rxns; // forward reactions
  std::vector<int> reverse_info;
  
// photolysis stuff
public:
  int nw;
  std::vector<double> wavl;
  
};




