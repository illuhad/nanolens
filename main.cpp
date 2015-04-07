/* 
 * File:   main.cpp
 * Author: aksel
 *
 * Created on 7. Dezember 2014, 14:55
 */

#include <iostream>
#include <iomanip>
#include "render_engine.hpp"
#include "timer.hpp"

/*
 * 
 */
int main(int argc, char** argv) 
{
  boost::mpi::environment env;
  boost::mpi::communicator world;

  nanolens::util::master_ostream master_cout(std::cout, 0);

  master_cout << "nanolens version 1.0e-22 launching..." << std::endl;
  master_cout << "Using " << world.size() << " process(es).\n";

  nanolens::system lensing_system("stars.dat", {1.0, 1.0});

  nanolens::util::vector2 npixels = {512, 512};
  //nanolens::util::vector2 physical_size = {0.07, 0.07};
  //nanolens::util::vector2 screen_position = {0.37, 0.63};
  
  nanolens::util::vector2 physical_size = {1.6, 1.6};
  nanolens::util::vector2 screen_position = {0.0, 0.0};

  master_cout << "Starting adaptive meshtracing..." << std::endl;

  nanolens::util::timer timer;
  timer.start();

  nanolens::render_engine tracer(npixels,
                                                      physical_size,
                                                      screen_position, 4, 30);

  std::cout.width(6);
  std::cout.precision(3);
  tracer.run(lensing_system, world, [&](double progress)
  { 
    master_cout << "\r" << progress * 100 << "%";
    std::cout.flush();
  });

  double time = timer.stop();
  
  master_cout << std::endl;
  world.barrier();

  std::cout << "Process " << world.rank() << " traced " << tracer.get_num_traced_rays() << " photons via "
          << tracer.get_num_domes() << " self-refining ray domes" << std::endl;

  world.barrier();
  
  std::cout.precision(6);
  std::cout << "Process "<< world.rank() << "'s average performance: " << 
          tracer.get_num_traced_rays() / time << " rays/s" << std::endl;

  world.barrier();
  
  master_cout << "Saving output..." << std::endl;

  if(world.rank() == 0)
    tracer.save_pixels("output.map");

  return 0;
}

