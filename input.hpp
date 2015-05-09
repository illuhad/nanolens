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


#ifndef INPUT_HPP
#define	INPUT_HPP

#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/mpi.hpp>
#include <boost/program_options.hpp>
#include "util.hpp"

namespace nanolens{

class configuration
{
  
public:
  class random_distribution_descriptor
  {
    
  public:
    
    enum distribution_type
    {
      UNIFORM,
      NORMAL
    };
   
    random_distribution_descriptor() = default;
    random_distribution_descriptor(distribution_type t,
                                   util::scalar c,
                                   util::scalar w)
      : _type(t), _center(c), _width(w){}
    
    distribution_type get_type() const{return _type; }
    util::scalar get_center() const{return _center; }
    util::scalar get_width() const{return _width; }
  private:
    distribution_type _type;
    util::scalar _center;
    util::scalar _width;
    
    
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & _type;
      ar & _center;
      ar & _width;
    }
  };
  
  struct random_star_generator_descriptor
  {
    random_distribution_descriptor x_distribution;
    random_distribution_descriptor y_distribution;
    random_distribution_descriptor mass_distribution;
    std::size_t num_stars;
        
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & x_distribution;
      ar & y_distribution;
      ar & mass_distribution;
      ar & num_stars;
    }
  };
  
  configuration(const boost::mpi::communicator& comm, int master_rank)
  : _comm(comm), _master_rank(master_rank) {}
  
  void load_from_file(const std::string& filename)
  {
    if(_comm.rank() == _master_rank)
    {
      boost::property_tree::xml_parser::read_xml(filename, _tree);

      _screen_pos = get_vector("nanolens.system.source_plane.position", util::vector2({0.0, 0.0}));
      _screen_size = get_vector("nanolens.system.source_plane.physical_size", util::vector2({10.0, 10.0}));
      _observer_pos = get_vector("nanolens.system.observer_plane.position", util::vector2({0.0, 0.0}));
      _dL = get<util::scalar>("nanolens.system.dL", 1.0);
      _dLS = get<util::scalar>("nanolens.system.dLS", 1.0);
      
      _resolution = get_vector("nanolens.system.source_plane.num_pixels", std::array<std::size_t, 2>({100, 100}));
      
      BOOST_FOREACH(boost::property_tree::ptree::value_type &v,
              _tree.get_child("nanolens.system.lens_plane"))
      {
        try
        {
          if(v.first == "star_generator")
          {
            std::string type = v.second.get<std::string>("<xmlattr>.type");

            if(type == "from_random_distribution")
            {
              random_star_generator_descriptor star_gen_descr;
              star_gen_descr.num_stars = v.second.get<std::size_t>("<xmlattr>.num_stars");
              star_gen_descr.x_distribution = get_random_distribution(v, "x");
              star_gen_descr.y_distribution = get_random_distribution(v, "y");
              star_gen_descr.mass_distribution = get_random_distribution(v, "mass");

              _random_star_generators.push_back(star_gen_descr);
            }
            else if(type == "from_file")
            {
              std::string filename = v.second.get<std::string>("<xmlattr>.filename");
              _star_files.push_back(filename);
            }
          }
        }
        catch(...)
        {}
      }
    }
    boost::mpi::broadcast(_comm, *this, _master_rank);
  }
  
  /*
  void load_from_command_line_args(int argc, char** argv)
  {
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("dL", boost::program_options::value<util::scalar>(), "distance between observer and lens plane")
        ("dLS", boost::program_options::value<util::scalar>(), "distance between lens and source plane")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }

    if (vm.count("compression")) {
        cout << "Compression level was set to " 
     << vm["compression"].as<int>() << ".\n";
    } else {
        cout << "Compression level was not set.\n";
    }
  }
   * */
  
  const std::array<std::size_t, 2>& get_resolution() const
  { return _resolution; }
  
  const util::vector2& get_screen_pos() const
  { return _screen_pos; }
  
  const util::vector2& get_physical_screen_size() const
  { return _screen_size; }
  
  const util::vector2& get_observer_pos() const
  { return _observer_pos; }
  
  util::scalar get_dLS() const
  { return _dLS; }
  
  util::scalar get_dL() const
  { return _dL; }
  
  const std::vector<std::string>& get_star_files() const
  { return _star_files; }
  
  const std::vector<random_star_generator_descriptor>&
  get_random_star_generators() const
  { return _random_star_generators; }
  
  
private:
  template<class T>
  T get(const std::string& identifier, const T& default_parameter)
  {
    try
    {
      return _tree.get<T>(identifier);
    }
    catch(...)
    {
      return default_parameter;
    }
  }
  
  template<class T>
  std::array<T,2> get_vector(const std::string& id, const std::array<T,2>& default_val)
  {
    std::array<T,2> result;
    result[0] = get<T>(id + ".<xmlattr>.x", default_val[0]);
    result[1] = get<T>(id + ".<xmlattr>.y", default_val[1]);
    return result;
  }
  
  util::vector2 get_vector(const boost::property_tree::ptree::value_type& v,
                           const std::string& id)
  {
    util::vector2 result;
    result[0] = v.second.get<util::scalar>(id + ".<xmlattr>.x");
    result[1] = v.second.get<util::scalar>(id + ".<xmlattr>.y");
    return result;
  }
  
  random_distribution_descriptor get_random_distribution(
                           const boost::property_tree::ptree::value_type& v,
                           const std::string& id)
  {
    random_distribution_descriptor::distribution_type type = 
      random_distribution_descriptor::NORMAL;
    
    util::scalar center = 0.0;
    util::scalar width = 1.0;
    
    std::string type_string = v.second.get<std::string>(id + ".<xmlattr>.distribution");
    if(type_string == "gaussian")
    {
      type = random_distribution_descriptor::NORMAL;
      center = v.second.get<util::scalar>(id + ".<xmlattr>.mean");
      width = v.second.get<util::scalar>(id + ".<xmlattr>.stddev");
    }
    else if(type_string == "uniform")
    {
      type = random_distribution_descriptor::UNIFORM;
      
      util::scalar min = v.second.get<util::scalar>(id + ".<xmlattr>.min");
      util::scalar max = v.second.get<util::scalar>(id + ".<xmlattr>.max");
      
      center = 0.5 * (min + max);
      width = std::abs(0.5 * (max - min));
    }
    else throw std::invalid_argument("invalid random distribution");
    
    return random_distribution_descriptor(type, center, width);
  }
  
  std::vector<random_star_generator_descriptor> _random_star_generators;
  
  std::vector<std::string> _star_files;
  
  boost::property_tree::ptree _tree;
  
  std::array<std::size_t, 2> _resolution;
  util::vector2 _screen_pos;
  util::vector2 _screen_size;
  util::vector2 _observer_pos;
  util::scalar _dL;
  util::scalar _dLS;
  boost::mpi::communicator _comm;
  int _master_rank;
  
  
  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & _resolution;
    ar & _screen_pos;
    ar & _screen_size;
    ar & _observer_pos;
    ar & _dL;
    ar & _dLS;
    ar & _star_files;
    ar & _random_star_generators;
  }
};

}

#endif	/* INPUT_HPP */

