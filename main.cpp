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
#include "system.hpp"
/*
 * 
 */
int main(int argc, char** argv) 
{
  boost::mpi::environment env;
  boost::mpi::communicator world;

  nanolens::util::master_ostream master_cout(std::cout, 0);

  master_cout << "nanolens version 2.0e-22 launching..." << std::endl;
  master_cout << "Using " << world.size() << " process(es).\n";

  nanolens::system lensing_system("stars.dat", {1.0, 1.0});

  std::array<std::size_t, 2> npixels = {512, 512};
  //nanolens::util::vector2 physical_size = {0.07, 0.07};
  //nanolens::util::vector2 screen_position = {0.37, 0.63};
  
  nanolens::util::vector2 physical_size = {1.6, 1.6};
  nanolens::util::vector2 screen_position = {0.0, 0.0};

  master_cout << "Starting computation..." << std::endl;

  nanolens::util::timer timer;
  timer.start();

  nanolens::render_engine engine(npixels,
                                 physical_size,
                                 screen_position,
                                 1.e-7);

  std::cout.width(6);
  std::cout.precision(5);
  
  std::size_t max_line_length = 0;
  
  auto status_handler = [&](const nanolens::status_info& status)
  {
    if(status.get_current_notification() != "")
    {
      std::cout << "Message from process " << world.rank() << ": " << status.get_current_notification() << std::endl;
    }
    
    if(world.rank() == 0)
    {
      if(status.get_schedule() != nullptr)
      {
        master_cout << "Performance metric of current task schedule:\n";

        std::vector<double> performances = status.get_schedule()->get_relative_performances();

        for(std::size_t i = 0; i < performances.size(); ++i)
          master_cout << "Process " << i << ": " << performances[i]/performances[0] << "x\n";
        
        master_cout << "(relative to process 0)\n";
      }

      std::string clear_string = "\r";
      for(std::size_t i = 0; i < max_line_length; ++i)
        clear_string += " ";

      master_cout << clear_string;
      std::stringstream sstr;
      sstr << "\r" << "Status: " << status.get_task() << " | " << progress * 100 << "%";
      std::string msg =  sstr.str();

      if(msg.length() > max_line_length)
        max_line_length = msg.length();

      master_cout << "\r" << msg;
      std::cout.flush();
    }
  }
  
  
  engine.run(lensing_system, world, status_handler);

  double time = timer.stop();
  
  master_cout << std::endl;
  master_cout << "Total elapsed time: " << time << "s\n";
  world.barrier();

  
  master_cout << "Saving output..." << std::endl;

  if(world.rank() == 0)
    engine.save_pixels("output.map");

  return 0;
}

