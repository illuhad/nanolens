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

#ifndef POSTPROCESSING_HPP
#define	POSTPROCESSING_HPP

#include <boost/mpi/communicator.hpp>
#include "input.hpp"
#include "status_handler.hpp"
#include "convolution.hpp"
#include "model.hpp"
#include "fits.hpp"

namespace nanolens {

class post_processing_launcher
{
public:
  post_processing_launcher(const boost::mpi::communicator& comm)
  : _comm(comm){}
  
  void execute_configuration(const configuration& config,
                              util::master_ostream& ostr) const
  
  {
    ostr << "\nStarting postprocessing.\n";
    if(config.get_post_processing_steps().empty())
    {
      ostr << "No postprocessing steps have been specified.\n";
    }
    
    std::size_t step_index = 0;
    for(const configuration::post_processing_step_descriptor& step 
      : config.get_post_processing_steps())
    {
      
      ostr << "Processing step " << step_index << std::endl;
      
      standard_terminal_output out(_comm, ostr.get_master_rank());
      
      std::string input = step.get_fits_input();
      std::string output = step.get_fits_output();
      
      if(output.empty())
        ostr << "Warning: No output specified for postprocessing step "
             << step_index << std::endl;
      
      if(input.empty())
      {
        // Try using the output of the lensing code as input
        ostr << "No input specified, using the lensing output file.\n";
        input = config.get_fits_output();
      }
      
      if(input.empty())
      {
        ostr << "No input file available. Skipping step.\n";
      }
      else
      {
        // Use the values from the lensing screen, if no input size has been specified
        // here
        util::vector2 input_center = config.get_config_node_vector2_property(
              step.get_config_node(),
              "input_position",
              config.get_screen_pos());
        
        util::vector2 input_sidelengths = config.get_config_node_vector2_property(
              step.get_config_node(),
              "input_physical_size",
              config.get_physical_screen_size());
        
        util::vector2 half_sidelength = input_sidelengths;
        util::scale(half_sidelength, 0.5);
        
        util::vector2 min_extent = input_center;
        util::vector2 max_extent = input_center;
        util::sub(min_extent, half_sidelength);
        util::add(max_extent, half_sidelength);
        
        switch(step.get_type())
        {
        case configuration::DIRECT_CONVOLUTION:
          run_convolution<direct_convolution<util::scalar,util::vector2>>(
            input, output, step, config, min_extent, max_extent, ostr, out);
          break;
        }
      }
      ++step_index;
    }
  }
  
private:
  template<class Convolution_type>
  void run_convolution(const std::string& input, 
                       const std::string& output,
                       const configuration::post_processing_step_descriptor& step,
                       const configuration& config,
                       const util::vector2& min_extent,
                       const util::vector2& max_extent,
                       util::master_ostream& ostr,
                       standard_terminal_output& out) const
  {
              
    Convolution_type conv(_comm, out);
    std::shared_ptr<models::basic_model<util::scalar, util::vector2>> model
      = get_kernel_model(config, step.get_config_node());

    // Execute convolution
    util::fits<util::scalar> input_file(input);
    util::multi_array<util::scalar> input_data;
    if(_comm.rank() == ostr.get_master_rank())
    {
      input_file.load<2>(input_data);
    }
    boost::mpi::broadcast(_comm, input_data, ostr.get_master_rank());

    util::multi_array<util::scalar> output_data;
    conv.run(input_data, min_extent, max_extent, *model, 
             model->get_evaluation_diameter(), output_data);

    if(!output.empty())
    {
      if(_comm.rank() == ostr.get_master_rank())
      {
        ostr << "Saving output...\n";
        util::fits<util::scalar> output_file(output);
        output_file.save(output_data);
      }
    }
  }
  
  std::shared_ptr<models::basic_model<util::scalar, util::vector2>> get_kernel_model(
                                          const configuration& config,
                                          const configuration::config_node& ppstep_node) const
  {
    auto kernel_node = ppstep_node.second.find("kernel");
    if(kernel_node != ppstep_node.second.not_found())
    {
      std::string id = config.get_config_node_type(*kernel_node, "");
      if(id == "")
        throw std::invalid_argument("No kernel type specified.");
      
      if(id == "shakura_sunyaev")
      {
        return std::shared_ptr<models::basic_model<util::scalar,util::vector2>>(
          new models::shakura_sunyaev(config, *kernel_node));
      }
      else throw std::invalid_argument("Unknwon convolution kernel: " + id);
    }
    else throw std::invalid_argument("No kernel specified.");
    
    return std::shared_ptr<models::basic_model<util::scalar, util::vector2>>();
  }
  
  boost::mpi::communicator _comm;
};

}

#endif	/* POSTPROCESSING_HPP */

