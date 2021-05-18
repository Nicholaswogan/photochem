#include <iostream>
#include "yaml-cpp/yaml.h"

void Log(std::string message)
{
  std::cout << message << std::endl;
}

int main()
{
    YAML::Node config = YAML::LoadFile("config.yaml");

    const std::string username = config["username"].as<std::string>();
    const std::string password = config["password"].as<std::string>();
    
    Log(username);
    Log(password);
    
    
}