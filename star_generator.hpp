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

#ifndef STAR_GENERATOR_HPP
#define	STAR_GENERATOR_HPP

#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include "star.hpp"
#include "input.hpp"

namespace nanolens{

class star_generator
{
public:
  star_generator(const boost::mpi::communicator& comm, int master_rank = 0)
  : _rand_generator(_rd()), _comm(comm), _master_rank(master_rank)
  {}
  
  void from_file(const std::string& filename, std::vector<star>& out)
  {
    out.clear();
    
    if(_comm.rank() == _master_rank)
    {
      std::ifstream input_file;

      input_file.open(filename.c_str());

      if(input_file.is_open())
      {
        while(input_file.good())
        {
          util::vector2 position = {0.0, 0.0};
          util::scalar mass = 0.0;

          input_file >> position[0];
          input_file >> position[1];
          input_file >> mass;

          star new_star(position, mass);
          
          if(input_file.good())
            out.push_back(new_star);
        }
      }
    }
    
    boost::mpi::broadcast(_comm, out, _master_rank);
    
    _star_list = out;
  }
  
  
  template<class Random_distribution_pos_x,
           class Random_distribution_pos_y,
           class Random_distribution_mass>
  void from_random_distribution(std::size_t num_stars,
                                std::vector<star>& out,
                                const Random_distribution_pos_x& rand_x,
                                const Random_distribution_pos_y& rand_y,
                                const Random_distribution_mass& rand_mass)
  {
    if(_comm.rank() == _master_rank)
    {
      out.clear();
      
      Random_distribution_pos_x x = rand_x;
      Random_distribution_pos_y y = rand_y;
      Random_distribution_mass m = rand_mass;
      
      for(std::size_t i = 0; i < num_stars; ++i)
      {
        util::vector2 position;
        position[0] = x(_rand_generator);
        position[1] = y(_rand_generator);

        util::scalar mass = std::abs(m(_rand_generator));

        out.push_back(star(position, mass));
      }
    }
    
    boost::mpi::broadcast(_comm, out, _master_rank);
    
    for(const star& s : out)
      _star_list.push_back(s);
  }
  
  void from_random_distribution(std::size_t num_stars,
                                std::vector<star>& out,
                                const configuration::random_distribution_descriptor& rand_x,
                                const configuration::random_distribution_descriptor& rand_y,
                                const configuration::random_distribution_descriptor& rand_mass)
  {
    if(_comm.rank() == _master_rank)
    {
      out.clear();
      
      std::vector<util::scalar> x_values;
      std::vector<util::scalar> y_values;
      std::vector<util::scalar> masses;
      
      generate_random_numbers(num_stars, rand_x, x_values);
      generate_random_numbers(num_stars, rand_y, y_values);
      generate_random_numbers(num_stars, rand_mass, masses);
      
      for(std::size_t i = 0; i < num_stars; ++i)
      {
        util::vector2 position;
        position[0] = x_values[i];
        position[1] = y_values[i];

        util::scalar mass = std::abs(masses[i]);

        out.push_back(star(position, mass));
      }
    }
    
    boost::mpi::broadcast(_comm, out, _master_rank);
    
    for(const star& s : out)
      _star_list.push_back(s);
  }
  
  // only works on the master process
  void save_generated_stars(const std::string& filename, bool append=false)
  {
    if(_comm.rank() == _master_rank)
    {
      std::ofstream file;
      if(!append)
        file.open(filename.c_str(), std::ios::out | std::ios::trunc);
      else file.open(filename.c_str(), std::ios::out | std::ios::app);
      
      if(file.is_open())
      {
        for(const star& s : _star_list)
        {
          file << s.get_position()[0] << " " 
               << s.get_position()[1] << " " 
               << s.get_mass() << "\n";
        }
        
        file.close();
      }
    }
  }
  
  const std::vector<star>& get_generated_stars() const
  { return _star_list; }
  
  void clear_generated_stars()
  { _star_list.clear(); }
private:
  
  inline void generate_random_numbers(std::size_t num_values,
                               const configuration::random_distribution_descriptor& descr,
                               std::vector<util::scalar>& out)
  {
    if(descr.get_type() == configuration::random_distribution_descriptor::NORMAL)
      generate_random_numbers(num_values,
                              std::normal_distribution<>(descr.get_center(), descr.get_width()),
                              out);
    else if(descr.get_type() == configuration::random_distribution_descriptor::UNIFORM)
      generate_random_numbers(num_values,
                          std::uniform_real_distribution<>(descr.get_center()-descr.get_width(), 
                                                           descr.get_center()+descr.get_width()),
                          out);  
      
  }
  
  template<class Distribution>
  void generate_random_numbers(std::size_t num_values,
                               const Distribution& distr,
                               std::vector<util::scalar>& out)
  {
    out.clear();
    
    Distribution d = distr;
    for(std::size_t i = 0; i < num_values; ++i)
    {
      out.push_back(d(_rand_generator));
    }
  }
  
  std::vector<star> _star_list;
  
  std::random_device _rd;
  std::mt19937 _rand_generator;
  boost::mpi::communicator _comm;
  int _master_rank;
};

}

#endif	/* STAR_GENERATOR_HPP */

