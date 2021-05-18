#include <iostream>
#include <vector>
#include <map>
#include "yaml-cpp/yaml.h"

class PhotoMechanism
{
public:
  
  int natoms;
  std::vector<std::string> atoms_names;
  std::vector<double> atoms_mass;
  
  int nsp;
  std::vector<std::string> species_names;
  std::vector<int> species_composition;
  std::vector<double> species_mass;
  std::vector<double> thermo_data;
  std::vector<double> thermo_temps;
  
  bool reverse;
  int nrF;
  int nrR;
  int nrT;
  std::vector<std::string> reactants_names;
  
  
  PhotoMechanism(std::string rxfile)
  {
    YAML::Node mechfile = YAML::LoadFile("../zahnle.yaml");
    
    // atoms
    YAML::Node atoms_node = mechfile["atoms"];
    natoms = atoms_node.size();
    atoms_names.resize(natoms);
    atoms_mass.resize(natoms);
    int i = 0;
    for (YAML::const_iterator it=atoms_node.begin(); it != atoms_node.end(); ++it)
    {
      atoms_names[i] = it->first.as<std::string>();
      atoms_mass[i] = it->second.as<double>();
      i++;
    }
    
    // species
    YAML::Node species_node = mechfile["species"];
    nsp = species_node.size();
    species_names.resize(nsp);
    species_mass.resize(nsp);
    std::fill(species_mass.begin(),species_mass.end(),0.0);
    species_composition.resize(nsp * natoms);
    thermo_temps.resize(3 * nsp);
    thermo_data.resize(7 * 2 * nsp);
    std::fill(thermo_temps.begin(),thermo_temps.end(),-1.0);
    for (int i = 0; i < species_node.size(); i++)
    {
      YAML::Node spnode = species_node[i];
      species_names[i] = spnode["name"].as<std::string>();
      YAML::Node comp = spnode["composition"];
      for(int j = 0; j < natoms; j++)
      {
        try
        {
          species_composition[j + i*natoms] = comp[atoms_names[j]].as<int>();
        }
        catch(YAML::BadConversion& e)
        {
          species_composition[j + i*natoms] = 0;
        }
        species_mass[i] += species_composition[j + i*natoms] * atoms_mass[j];
      }
      
      YAML::Node therm = spnode["thermo"];
      if (therm["model"].as<std::string>() != "Shomate")
        throw "Shomate is the only allowed data type";      
      for (int j = 0; j < therm["temperature-ranges"].size(); j++)
        thermo_temps[j + i*3] = therm["temperature-ranges"][j].as<double>();
      
      if (therm["temperature-ranges"].size()-1 != therm["data"].size())
        throw "Problem here";
      for (int j = 0; j < therm["data"].size(); j++)
      {
        if (therm["data"][j].size() != 7)
          throw "Problem here!!";
        for (int k = 0; k < 7; k++)
        {
          thermo_data[k + 7*(j + 2*i)] = therm["data"][j][k].as<double>();
        }
      }  
    }
    
    // reactions
    nrF = mechfile["reactions"].size();
    for (int i = 0; i < nrF; i++)
    {
      std::cout << mechfile["reactions"][i]["equation"].as<std::string>() << std::endl;
    }
  }
};

int main()
{

  // PhotoMechanism photomech("../zahnle.yaml");
  // std::cout << photomech.atoms_names.size() << std::endl;
  
  std::string a = "HO2 + HO2 (+ M) <=> H2O2 + O2 (+ M)";
  
  if (a.find("<=>") != std::string::npos)
  {
      
  }
  else if (a.find(" =>") != std::string::npos)
  {
    
  }
  
  
  

  
}