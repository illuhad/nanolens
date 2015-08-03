/*
 * This file is part of nanolens, a free program to calculate microlensing 
 * magnification patterns.
 * Copyright (C) 2015  Aksel Alpay
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <iomanip>
#include "render_engine.hpp"
#include "timer.hpp"
#include "system.hpp"
#include "star_generator.hpp"
#include "input.hpp"
#include "status_handler.hpp"
#include "launch.hpp"




int main(int argc, char** argv) 
{ 
  boost::mpi::environment env;
  boost::mpi::communicator world;

  nanolens::util::master_ostream master_cout(std::cout, 0);

  master_cout << "nanolens Copyright (C) 2015 Aksel Alpay\n"
    "This program comes with ABSOLUTELY NO WARRANTY; It is free software,\n"
    "and you are welcome to redistribute it under the conditions of the\n"
    "GNU General Public License v3.\n\n"
"                                                                           \n" 
"                                                                           \n"    
"                             ````                                          \n"    
"                             `````                                         \n"    
"                             ``..``                                        \n"    
"                              ``/.`                                        \n"    
"                               `.o`                                        \n"    
"                                `hs:`                                      \n"    
"                                -o::+:.                                    \n"    
"                                +-....-/:-.`               ````````        \n"    
"                               /-.........-:::----....--::-:-.```          \n"    
"                              /-..`````````.......--:+o/.`                 \n"    
"                            `/-..`````````````.....::.                     \n"    
"                           ::..```````````````...::`                       \n"    
"                         -/-..```````````````..-/`                         \n"    
"                      `:/-....``````````````..:-                           \n"    
"                 ``./oo/::--......````````.../.                            \n"    
"           ````.----.`````...---::-.........+.                             \n"    
"         ````````                 .-:/-...-/:                              \n"    
"                                      .:+:/o                               \n"    
"                                        `/y+`                              \n"    
"                                         `.o``                             \n"    
"                                         ``:-``                            \n"    
"                                          ``.```                           \n"    
"                                          ``````                           \n"    
"                                            ```                            \n";
                                                                                
                                                                                

  
  master_cout << "nanolens version 2.9e-19 (beta) launching..." << std::endl;
  master_cout << "Using " << world.size() << " process(es).\n";
  
  std::vector<nanolens::star> stars;
  
  master_cout << "Preparing system...\n";
  
  nanolens::configuration config(world, 0);
  config.load_from_file("nanolens.xml");
  
  // Generate stars
  nanolens::star_generator star_gen(world);
  
  for(const std::string& filename : config.get_star_files())
  {
    master_cout << "Star genesis: Processing file: " << filename << std::endl;
    std::vector<nanolens::star> generated_stars;
    star_gen.from_file(filename, stars); 
    
    for(const nanolens::star& s : generated_stars)
      stars.push_back(s);
  }
  
  for(const nanolens::configuration::random_star_generator_descriptor& descr : 
      config.get_random_star_generators())
  {
    master_cout << "Star genesis: Generating " << descr.num_stars << " random stars...\n";
    
    std::vector<nanolens::star> generated_stars;
    
    star_gen.from_random_distribution(descr.num_stars,
                                      generated_stars,
                                      descr.x_distribution,
                                      descr.y_distribution,
                                      descr.mass_distribution);
    
    for(const nanolens::star& s : generated_stars)
      stars.push_back(s);
  }
  master_cout << "Star genesis: Created " << stars.size() << " stars." << std::endl;
  
  star_gen.save_generated_stars("nanolens_star_log.dat");

  std::shared_ptr<nanolens::lens_plane> deflector(new nanolens::lens_plane(stars,
                                                                           config.get_shear(),
                                                                           config.get_sigma_smooth()));
  
  nanolens::system<nanolens::lens_plane> lensing_system(deflector);
  
  master_cout << "Lensing system initialized:\n";
  
  std::map<std::string, nanolens::util::scalar> lens_plane_statistics;
  lensing_system.get_deflector().obtain_properties_set(lens_plane_statistics);
  
  for(const auto& element : lens_plane_statistics)
    master_cout << "  " << element.first << " = " << element.second << std::endl;
  
  master_cout << "Starting computation..." << std::endl;

  nanolens::util::timer timer;
  timer.start();

  nanolens::standard_launcher<nanolens::system<nanolens::lens_plane>> launcher;
  launcher.run_configured_method(world, config, lensing_system);

  double time = timer.stop();
  
  master_cout << std::endl;
  master_cout << "Total elapsed time: " << time << "s\n";
  world.barrier();
  
  
  return 0;
}

