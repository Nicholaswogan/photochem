#include <iostream>
#include <vector>
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


std::vector<std::string> SplitString(std::string& s, std::string& delim)
{
  std::vector<std::string> out(0);
  int start = 0;
  int end = s.find(delim);
  while(end != std::string::npos)
  {
    out.push_back(s.substr(start, end - start));
    start = end + delim.length();
    end = s.find(delim,start);
  }
  out.push_back(s.substr(start, end));
  return out;
}

void RemoveSubstrings(std::string& s, std::string& p) { 
  std::string::size_type n = p.length();
  for (int i = s.find(p); i != std::string::npos; i = s.find(p))
      s.erase(i, n);
}

// trim from start (in place)
void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}


void parse_reaction(std::string& rx, 
                    std::vector<std::string>& reacts, 
                    std::vector<std::string>& prods)
{
  std::string leftp = "(";
  std::string rightp = ")";
  RemoveSubstrings(rx, leftp);
  RemoveSubstrings(rx, rightp);
  
  std::string delim1 = "<=>";
  std::string delim2 = " =>";
  std::vector<std::string> out;
  
  if (rx.find(delim1) != std::string::npos)
  {
    out = SplitString(rx, delim1);
    if (out.size() != 2)
      std::cout << "Problem" << std::endl;
  }
  else if (rx.find(delim2) != std::string::npos)
  {
    out = SplitString(rx, delim2);
    if (out.size() != 2)
      std::cout << "Problem" << std::endl;
  }
  else
  {
    std::cout << "Problem" << std::endl;
    // return;
  }
  
  std::string plus = "+";
  reacts = SplitString(out[0], plus);
  prods = SplitString(out[1], plus);
  
  for(int i = 0; i < reacts.size(); i++)
    trim(reacts[i]);
  for(int i = 0; i < prods.size(); i++)
    trim(prods[i]);
}


int main()
{

  PhotoMechanism photomech("../zahnle.yaml");
  std::cout << photomech.atoms_names.size() << std::endl;
  
  
  std::string a = "HO2 + HO2 (+ M) <=> H2O2 + O2 (+ M)";
  std::vector<std::string> reacts, prods;
  
  parse_reaction(a, reacts, prods);
  
  // std::cout << prods[0] << std::endl;
}




