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

#ifndef LAUNCH_HPP
#define	LAUNCH_HPP

#include "input.hpp"
#include "method_irs.hpp"
#include "star_generator.hpp"
#include "lens_plane.hpp"
#include "fragmented_lens_plane.hpp"

namespace nanolens {

template<class System_type>
class standard_method_launcher
{
public:
  standard_method_launcher(util::master_ostream* ostr)
  : _ostr(ostr)
  {}
  
  standard_method_launcher()
  : _ostr(nullptr) {}
  
  template<class Method_type>
  void run_method(const boost::mpi::communicator& comm,
                           const configuration& config,
                           System_type& sys) const
  {
    if(_ostr != nullptr)
    {
      *_ostr << "Using lensing system:\n";

      std::map<std::string, util::scalar> lens_plane_statistics;
      sys.get_deflector().obtain_properties_set(lens_plane_statistics);

      for(const auto& element : lens_plane_statistics)
        *_ostr << "  " << element.first << " = " << element.second << std::endl;
      
      *_ostr << "\n";
    }
    
    typename Method_type::settings method_settings(config);

    nanolens::render_engines::standard_renderer
    <
      System_type, 
      Method_type,
      util::scalar
    > engine(comm,
             method_settings,
             config.get_resolution(),
             config.get_physical_screen_size(),
             config.get_screen_pos());

    nanolens::standard_terminal_output status_handler(comm, 0);

    engine.run(sys, status_handler);
    
    if(_ostr != nullptr)
      *_ostr << "Saving output...\n";
    
    if(!config.get_fits_output().empty())
      engine.save_as_fits(config.get_fits_output());
  
    if(!config.get_raw_output().empty())
      engine.save_as_raw(config.get_raw_output());

  }


  void run_configured_method(const boost::mpi::communicator& comm,
                               const configuration& config,
                               System_type& sys) const
  {
    switch(config.get_method_type())
    {
    case nanolens::configuration::INVERSE_RAY_SHOOTING:
      run_method
      <
        magnification_pattern_generation::inverse_ray_shooting<System_type>
      >(comm, config, sys);
      break;
    }
  }
private:
  util::master_ostream* _ostr;
};

class standard_launcher
{
public:
  standard_launcher(const boost::mpi::communicator& comm)
  : _comm(comm)
  {}
  
  void execute_configuration(const configuration& config,
                             util::master_ostream& ostr) const
  {
    util::timer timer;
    timer.start();
    
    switch(config.get_lens_plane_type())
    {
    case configuration::MICROLENSING:
      {

        configuration::deflection_engine_type defl_engine
          = config.get_deflection_engine_type();
        switch(defl_engine)
        {
        case configuration::TREE:
          run_with_microlensing_lens_plane<tree_deflector>(config, ostr);
          break;
        case configuration::EXACT:
          run_with_microlensing_lens_plane<exact_deflector>(config, ostr);
          break;
        }
      }
      break;
    case configuration::FRAGMENTED_MICROLENSING:
      {
        configuration::deflection_engine_type defl_engine
          = config.get_deflection_engine_type();
        switch(defl_engine)
        {
        case configuration::TREE:
          run_with_fragmented_microlensing_lens_plane<tree_deflector>(config, ostr);
          break;
        case configuration::EXACT:
          run_with_fragmented_microlensing_lens_plane<exact_deflector>(config, ostr);
          break;
        }
      }
      break;
    }
    
    double time = timer.stop();

    ostr << std::endl;
    ostr << "Total elapsed time: " << time << "s\n";
    _comm.barrier();
  }
  
  
private:
  template<class Deflection_engine_type>
  void run_with_microlensing_lens_plane(const configuration& config,
                                        util::master_ostream& ostr) const
  {
    star_generator star_gen(_comm);
    
    std::vector<star> stars;
  
    for(const std::string& filename : config.get_star_files())
    {
      ostr << "Star genesis: Processing file: " << filename << std::endl;
      std::vector<star> generated_stars;
      star_gen.from_file(filename, stars); 

      for(const star& s : generated_stars)
        stars.push_back(s);
    }
  
    for(const configuration::random_star_generator_descriptor& descr : 
        config.get_random_star_generators())
    {
      ostr << "Star genesis: Generating " << descr.num_stars << " random stars...\n";

      std::vector<star> generated_stars;

      star_gen.from_random_distribution(descr.num_stars,
                                        generated_stars,
                                        descr.x_distribution,
                                        descr.y_distribution,
                                        descr.mass_distribution,
                                        descr.circularize,
                                        descr.circularization_radius);

      for(const star& s : generated_stars)
        stars.push_back(s);
    }
    ostr << "Star genesis: Created " << stars.size() << " stars." << std::endl;
  
    star_gen.save_generated_stars("nanolens_star_log.dat");

    std::shared_ptr<microlensing_lens_plane<Deflection_engine_type>> deflector_plane(
      new microlensing_lens_plane<Deflection_engine_type>(stars,
                                                  typename Deflection_engine_type::settings(config),
                                                  config.get_shear(),
                                                  config.get_sigma_smooth()));

    system<microlensing_lens_plane<Deflection_engine_type>> lensing_system(deflector_plane);
    
    standard_method_launcher<system<microlensing_lens_plane<Deflection_engine_type>>> launcher(&ostr);
    launcher.run_configured_method(_comm, config, lensing_system);

  }
  
  template<class Deflection_engine_type>
  void run_with_fragmented_microlensing_lens_plane(const configuration& config,
                                        util::master_ostream& ostr) const
  {
    star_generator star_gen(_comm);
    
    util::vector2 db_cell_size = config.get_lens_plane_vector2_property("star_db_cell_size", 
                                                                        util::vector2({10.0, 10.0}));
                                                                        
    std::size_t num_cells = config.get_lens_plane_property("star_db_expected_num_cells", 40000);
    
    typename horizontally_fragmented_microlensing_lens_plane<Deflection_engine_type>::star_db_type
      star_db("star_db", db_cell_size, num_cells, _comm, 0);
    
    std::size_t star_count = 0;
  
    for(const std::string& filename : config.get_star_files())
    {
      ostr << "Iterative star genesis: Processing file: " << filename << std::endl;
      std::vector<star> loaded_stars;
      star_gen.from_file(filename, loaded_stars); 

      for(const star& s : loaded_stars)
      {
        star_db.insert(s.get_position(), s.get_mass());
        ++star_count;
      }
    }
  
    for(const configuration::random_star_generator_descriptor& descr : 
        config.get_random_star_generators())
    {
      ostr << "Iterative star genesis: Generating " << descr.num_stars << " random stars...\n";

      for(std::size_t i = 0; i < descr.num_stars; ++i)
      {
        star current_star = star_gen.from_random_distribution(
                                          descr.x_distribution,
                                          descr.y_distribution,
                                          descr.mass_distribution,
                                          descr.circularize,
                                          descr.circularization_radius);
        
        star_db.insert(current_star.get_position(), current_star.get_mass());
        ++star_count;
      }

    }
  
    star_db.commit();

    ostr << "Committed " << star_count << " stars to the spatial database.\n";
    
    util::scalar lens_plane_y_center = config.get_lens_plane_property("lens_plane_y_center", 0.0);
    
    util::scalar fragment_size = config.get_lens_plane_property("fragment_size", 1.0);
    
    util::scalar fragment_star_distribution_radius = config.get_lens_plane_property(
                                                              "fragment_star_inclusion_radius", 
                                                              10.0);
    
    std::shared_ptr<horizontally_fragmented_microlensing_lens_plane<Deflection_engine_type>> deflector_plane(
      new horizontally_fragmented_microlensing_lens_plane<Deflection_engine_type>(
                                                  star_db,
                                                  typename Deflection_engine_type::settings(config),
                                                  lens_plane_y_center,
                                                  config.get_shear(),
                                                  config.get_sigma_smooth(),
                                                  fragment_size,
                                                  fragment_star_distribution_radius));

    system<horizontally_fragmented_microlensing_lens_plane<Deflection_engine_type>> lensing_system(deflector_plane);
    
    standard_method_launcher
    <
      system
      <
        horizontally_fragmented_microlensing_lens_plane
        <
          Deflection_engine_type
        >
      >
    > launcher(&ostr);
    
    launcher.run_configured_method(_comm, config, lensing_system);

  }
  
  boost::mpi::communicator _comm;
};


}

#endif	/* LAUNCH_HPP */

