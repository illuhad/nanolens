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
                                                                                
                                                                                

  
  master_cout << "nanolens version 2.96e-19 launching..." << std::endl;
  master_cout << "Using " << world.size() << " process(es).\n";
  
  std::vector<nanolens::star> stars;
  
  master_cout << "Preparing system...\n";
  
  nanolens::configuration config(world, 0);
  config.load_from_file("nanolens.xml");
  
  nanolens::standard_launcher launcher(world);
  
  launcher.execute_configuration(config, master_cout);
  
  
  return 0;
}

