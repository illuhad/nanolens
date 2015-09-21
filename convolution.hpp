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

#ifndef CONVOLUTION_HPP
#define	CONVOLUTION_HPP

#include "util.hpp"
#include "grid.hpp"
#include "status.hpp"
#include "scheduler.hpp"
#include "fits.hpp"

namespace nanolens{

template<typename Scalar_type, class Vector_type>
class convolution
{
public:
  convolution(const boost::mpi::communicator& comm, status_handler_type handler)
  : _comm(comm), _handler(handler) {}
  
  virtual ~convolution(){}
  
  typedef util::multi_array<Scalar_type>                  discrete_convolution_kernel_type;
  typedef util::multi_array<Scalar_type> array_type;
  
  template<class Functional_convolution_kernel_type>
  void run(const array_type& input, 
           const util::vector2& input_min_extent,
           const util::vector2& input_max_extent,
           const Functional_convolution_kernel_type& kernel,
           util::scalar kernel_evaluation_diameter,
           array_type& output)
  {
    // Start by filling a grid with evaluations of the kernel function
    util::vector2 pixel_size = input_max_extent;
    util::sub(pixel_size, input_min_extent);
    for(std::size_t i = 0; i < pixel_size.size(); ++i)
      pixel_size[i] /= static_cast<util::scalar>(input.get_extent_of_dimension(i));
    
    std::array<std::size_t, 2> num_kernel_pixels;
    for(std::size_t i = 0; i < pixel_size.size(); ++i)
      num_kernel_pixels[i] = static_cast<std::size_t>(kernel_evaluation_diameter / pixel_size[i]);

    util::vector2 expected_kernel_diameter = {num_kernel_pixels[0] * pixel_size[0],
                                              num_kernel_pixels[1] * pixel_size[1]};
    
    array_type kernel_array(num_kernel_pixels[0], num_kernel_pixels[1]);
    
    util::vector2 kernel_min_extent = {-expected_kernel_diameter[0] / 2, 
                                       -expected_kernel_diameter[1] / 2};
    
    util::vector2 kernel_max_extent = { expected_kernel_diameter[0] / 2, 
                                        expected_kernel_diameter[1] / 2};
    
    util::grid_translator<util::scalar, 2, false> translator(kernel_min_extent,
                                                             kernel_max_extent, 
                                                             num_kernel_pixels);
    
    for(std::size_t x = 0; x < kernel_array.get_extent_of_dimension(0); ++x)
      for(std::size_t y = 0; y < kernel_array.get_extent_of_dimension(1); ++y)
      {
        util::grid_translator<util::scalar, 2, false>::index_type pixel_index = {x,y};
        util::vector2 pixel_position = translator.get_central_position_of_bucket(pixel_index);
        
        kernel_array[pixel_index.data()] = kernel(pixel_position);
      }
    
    assert(kernel_array.get_extent_of_dimension(0) 
            == kernel_array.get_extent_of_dimension(1));
    
    run(input, kernel_array, output);
  }
  
  virtual void run(const array_type& input, 
                   const discrete_convolution_kernel_type& kernel,
                   array_type& output) = 0;
  
protected:
  boost::mpi::communicator _comm;
  status_handler_type _handler;
};

template<typename Scalar_type, class Vector_type>
class direct_convolution : public convolution<Scalar_type, Vector_type>
{
public:
  using convolution<Scalar_type, Vector_type>::run;
  using typename convolution<Scalar_type, Vector_type>::discrete_convolution_kernel_type;
  using typename convolution<Scalar_type, Vector_type>::array_type;
  typedef std::array<std::size_t, 2> pixel_index_type;
  
  direct_convolution(const boost::mpi::communicator& comm,
                     status_handler_type handler)
  : convolution<Scalar_type, Vector_type>(comm, handler)
  {}
  
  virtual void run(const array_type& input, 
                   const discrete_convolution_kernel_type& kernel,
                   array_type& output)
  {
    if(this->_comm.rank() == 1)
    {
      util::fits<util::scalar> kernel_file("kernel.fits");
      kernel_file.save(kernel);
    }
    assert(input.get_dimension() == 2);
    assert(kernel.get_dimension() == 2);
    
    output = array_type(input.get_extent_of_dimension(0), input.get_extent_of_dimension(1));
    std::fill(output.begin(), output.end(), 0.0);
   
    this->_handler(status_info("Scheduling direct convolution"));
    
    scheduler job_scheduler(this->_comm);
    util::scalar dummy_sum = 0.0;
    job_scheduler.autosized_run(input.get_extent_of_dimension(0), 
      [&](std::size_t test_ray_index, std::size_t num_tests)
      {
        pixel_index_type current_test_pixel = 
                    {test_ray_index % input.get_extent_of_dimension(0), 0};
        
        dummy_sum += get_convolved_pixel_value(input, kernel, current_test_pixel);
      }, 
    2.0);
      
    this->_handler(status_info("Scheduling complete", &job_scheduler, ""));
    
    std::size_t num_jobs_completed = 0;
    for(std::size_t x = job_scheduler.get_ownership_range_begin();
            x <= job_scheduler.get_ownership_range_end(); ++x)
    {
      util::scalar progress = static_cast<util::scalar>(num_jobs_completed) /
        static_cast<util::scalar>(job_scheduler.get_num_assigned_jobs());
      
      this->_handler(status_info("Convolving", nullptr, "", progress));
      
      for(std::size_t y = 0; y < input.get_extent_of_dimension(1); ++y)
      {
        pixel_index_type current_pixel = {x,y};
        output[current_pixel.data()] = 
                get_convolved_pixel_value(input, kernel, current_pixel);
      }
      
      ++num_jobs_completed;
    }
    this->_handler(status_info("", nullptr, "Finished computation"));
    this->_handler(status_info("Waiting for processes...\n"));
    
    auto add =  [](util::scalar& a, util::scalar b){ a += b;};
    output.reduce_parallel_array(this->_comm, 0, add);
  }
  
private:
  util::scalar get_convolved_pixel_value(const array_type& input,
                                         const discrete_convolution_kernel_type& kernel,
                                         const pixel_index_type& px) const
  {
    util::scalar sum = 0.0;
    
    pixel_index_type kernel_center = {kernel.get_extent_of_dimension(0) / 2,
                                      kernel.get_extent_of_dimension(1) / 2};
    
    for(std::size_t i = 0; i < kernel.get_extent_of_dimension(0); ++i)
    {
      std::size_t i_flipped = kernel.get_extent_of_dimension(0) - i - 1;
      for(std::size_t j = 0; j < kernel.get_extent_of_dimension(1); ++j)
      {
        std::size_t j_flipped = kernel.get_extent_of_dimension(1) - j - 1;
        
        long long int input_pixel0 = 
              static_cast<long long int>(px[0] + i - kernel_center[0]);
        long long int input_pixel1 =
              static_cast<long long int>(px[1] + j - kernel_center[1]);
        
        pixel_index_type kernel_pixel = {i_flipped, j_flipped};
        
        if(input_pixel0 >= 0 && 
           input_pixel0 < static_cast<long long int>(input.get_extent_of_dimension(0)) &&
           input_pixel1 >= 0 && 
           input_pixel1 < static_cast<long long int>(input.get_extent_of_dimension(1)))
        {
          pixel_index_type input_pixel = {static_cast<std::size_t>(input_pixel0), 
                                          static_cast<std::size_t>(input_pixel1)};
          sum += input[input_pixel.data()] * kernel[kernel_pixel.data()];
        }
      }
    }
    
    return sum;
  }
};

template<typename Scalar_type, class Vector_type>
class fftw_convolution : public convolution<Scalar_type, Vector_type>
{
public:
  using typename convolution<Scalar_type, Vector_type>::discrete_convolution_kernel_type;
  using typename convolution<Scalar_type, Vector_type>::array_type;
  
  fftw_convolution(const boost::mpi::communicator& comm,
                   status_handler_type handler)
  : convolution<Scalar_type, Vector_type>(comm, handler)
  {}
  
  virtual void run(const array_type& input, 
                   const discrete_convolution_kernel_type& kernel,
                   array_type& output)
  {
    assert(input.get_dimension() == 2);
    assert(kernel.get_dimension() == 2);
    
    output = array_type(input.get_extent_of_dimension(0), input.get_extent_of_dimension(1));
  }
  

};

}


#endif	/* CONVOLUTION_HPP */

