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
#include <boost/property_tree/ptree_serialization.hpp>
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
  typedef const boost::property_tree::ptree::value_type config_node;
  
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
    
  };
  
  struct random_star_generator_descriptor
  {
    random_distribution_descriptor x_distribution;
    random_distribution_descriptor y_distribution;
    random_distribution_descriptor mass_distribution;
    std::size_t num_stars;
    
    bool circularize;
    util::scalar circularization_radius;
        
  };
  
  enum post_processing_step_type
  {
    DIRECT_CONVOLUTION
  };
  
  class post_processing_step_descriptor
  {
  public:
    typedef const boost::property_tree::ptree::value_type config_node;
    
    post_processing_step_descriptor()
    {}
    
    post_processing_step_descriptor(const config_node& v, const configuration& config)
    : _config_node(v)
    {
      _fits_input = config.get_property(v, "fits_input", std::string());
      _fits_output = config.get_property(v, "fits_output", std::string());
      
      std::string type = config.get_config_node_type(v, "");
      
      if(type == "direct_convolution")
        _type = DIRECT_CONVOLUTION;
      else throw std::invalid_argument("Invalid postprocessing type: "+type);
    }
    
    const std::string& get_fits_input() const
    { return _fits_input; }
    
    const std::string& get_fits_output() const
    { return _fits_output; }
    
    post_processing_step_type get_type() const
    { return _type; }
    
    const config_node& get_config_node() const
    { return _config_node; }
  private:
    std::string _fits_input;
    std::string _fits_output;
    post_processing_step_type _type;
    config_node _config_node;
  };
  
  enum method_type
  {
    INVERSE_RAY_SHOOTING = 0
  };
  
  enum lens_plane_type
  {
    MICROLENSING = 0,
    FRAGMENTED_MICROLENSING
  };
  
  enum deflection_engine_type
  {
    EXACT = 0,
    TREE
  };
  
  configuration(const boost::mpi::communicator& comm, int master_rank)
  : _comm(comm), _master_rank(master_rank) {}
  
  void load_from_file(const std::string& filename)
  {
    boost::property_tree::xml_parser::read_xml(filename, _tree);
    boost::mpi::broadcast(_comm, _tree, _master_rank);

    _screen_pos = get_vector2_property("nanolens.screen.position", util::vector2({0.0, 0.0}));
    _screen_size = get_vector2_property("nanolens.screen.physical_size", util::vector2({10.0, 10.0}));

    _resolution = get_vector2_property("nanolens.screen.num_pixels", std::array<std::size_t, 2>({100, 100}));

    _fits_output = get_property<std::string>("nanolens.fits_output", std::string());
    _raw_output = get_property<std::string>("nanolens.raw_output", std::string());

    load_method();
    load_lens_plane_type();

    BOOST_FOREACH(boost::property_tree::ptree::value_type &v,
            _tree.get_child("nanolens.system.lens_plane"))
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

          if(v.second.find("circularize") != v.second.not_found())
          {
            star_gen_descr.circularize = true;
            star_gen_descr.circularization_radius = v.second.get<util::scalar>("circularize.<xmlattr>.radius");
          }
          else
          {
            star_gen_descr.circularize = false;
            star_gen_descr.circularization_radius = 0.0;
          }

          _random_star_generators.push_back(star_gen_descr);
        }
        else if(type == "from_file")
        {
          std::string filename = v.second.get<std::string>("<xmlattr>.filename");
          _star_files.push_back(filename);
        }
      }
    }

    _shear = get_property<util::scalar>("nanolens.system.lens_plane.shear", 0.0);
    _sigma_smooth = get_property<util::scalar>("nanolens.system.lens_plane.sigma_smooth", 0.0);
    _shear_rotation_angle = get_property<util::scalar>("nanolens.system.lens_plane.shear_rotation_angle", 0.0);

    if(_tree.find("nanolens.post_processing") != _tree.not_found())
    {
      BOOST_FOREACH(boost::property_tree::ptree::value_type &v,
              _tree.get_child("nanolens.post_processing"))
      {
        if(v.first == "step")
        {
          post_processing_step_descriptor step(v, *this);
          this->_post_processing_steps.push_back(step);
        }
      }
    }
  }
  
  const std::array<std::size_t, 2>& get_resolution() const
  { return _resolution; }
  
  const util::vector2& get_screen_pos() const
  { return _screen_pos; }
  
  const util::vector2& get_physical_screen_size() const
  { return _screen_size; }
  
  const std::vector<std::string>& get_star_files() const
  { return _star_files; }
  
  const std::vector<random_star_generator_descriptor>&
  get_random_star_generators() const
  { return _random_star_generators; }
  
  util::scalar get_shear() const
  { return _shear; }
  
  util::scalar get_shear_rotation_angle() const
  { return _shear_rotation_angle; }
  
  util::scalar get_sigma_smooth() const
  { return _sigma_smooth; }
  
  const std::string& get_fits_output() const
  { return _fits_output; }
  
  const std::string& get_raw_output() const
  { return _raw_output; }
  
  template<class T>
  T get_method_property(const std::string& identifier, const T& default_value) const
  { return get_namespace_property("method", identifier, default_value); }
  
  template<class T>
  std::array<T,2> get_method_vector2_property(const std::string& id, 
                                              const std::array<T,2>& default_val) const
  { return get_namespace_vector2_property("method", id, default_val); }
  
  template<class T>
  T get_lens_plane_property(const std::string& identifier, const T& default_value) const
  { return get_namespace_property("system.lens_plane", identifier, default_value); }
  
  template<class T>
  std::array<T,2> get_lens_plane_vector2_property(const std::string& id,
                                                        const std::array<T,2>& default_val) const
  { return get_namespace_vector2_property("system.lens_plane", id, default_val); }
  
  
  
  template<class T>
  T get_deflection_engine_property(const std::string& identifier, const T& default_value) const
  { return get_namespace_property("system.lens_plane.deflection_engine", identifier, default_value); }
  
  template<class T>
  std::array<T,2> get_deflection_engine_vector2_property(const std::string& id,
                                                        const std::array<T,2>& default_val) const
  { return get_namespace_vector2_property("system.lens_plane.deflection_engine", id, default_val); }
  
  deflection_engine_type get_deflection_engine_type() const
  {
    std::string type 
      = get_namespace_property<std::string>("system.lens_plane.deflection_engine", "<xmlattr>.type", "exact");
    
    if(type == "exact")
      return EXACT;
    else if(type == "tree")
      return TREE;
    else throw std::invalid_argument("Invalid deflection engine: "+type);
  }
  
  template<class T>
  T get_namespace_property(const std::string& nspace, const std::string& property, const T& default_value) const
  {
    return get_property("nanolens."+nspace+"."+property, default_value);
  }
  
  template<class T>
  std::array<T,2> get_namespace_vector2_property(const std::string& nspace, 
                                                 const std::string& property,
                                                 const std::array<T,2>& default_val) const
  {
    return get_vector2_property("nanolens."+nspace+"."+property, default_val);
  }
  
  template<class T>
  T get_config_node_property(const config_node& node, const std::string& property, const T& default_value) const
  {
    return get_property(node, property, default_value);
  }
  
  std::string get_config_node_type(const config_node& node, const std::string& default_value) const
  {
    return get_property(node, "<xmlattr>.type", default_value);
  }
  
  template<class T>
  std::array<T,2> get_config_node_vector2_property(const config_node& node, 
                                                 const std::string& property,
                                                 const std::array<T,2>& default_val) const
  {
    return get_vector2_property(node, property, default_val);
  }
  
  method_type get_method_type() const
  {
    return _method;
  }
  
  lens_plane_type get_lens_plane_type() const
  {
    return _lens_plane_type;
  }
  
  const std::vector<post_processing_step_descriptor>& get_post_processing_steps() const
  {
    return _post_processing_steps;
  }
private:
  void load_method()
  {
    std::string method_string = get_property<std::string>(
      "nanolens.method.<xmlattr>.type",
      "inverse_ray_shooting");
    
    if(method_string == "inverse_ray_shooting")
      _method = INVERSE_RAY_SHOOTING;
    else
      throw std::invalid_argument("Invalid method: " + method_string);
  }
  
  void load_lens_plane_type()
  {
    std::string type_string = get_property<std::string>("nanolens.system.lens_plane.<xmlattr>.type",
                                               "microlensing_lens_plane");
    
    if(type_string == "microlensing_lens_plane")
      _lens_plane_type = MICROLENSING;
    else if(type_string == "fragmented_microlensing_lens_plane")
      _lens_plane_type = FRAGMENTED_MICROLENSING;
    else
      throw std::invalid_argument("Invalid lens plane type: " + type_string);
  }
  
  template<class T>
  T get_property(const std::string& identifier, const T& default_parameter) const
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
  T get_property(const config_node& node, const std::string& identifier, const T& default_parameter) const
  {
    try
    {
      return node.second.get<T>(identifier);
    }
    catch(...)
    {
      return default_parameter;
    }
  }
  
  template<class T>
  std::array<T,2> get_vector2_property(const std::string& id, const std::array<T,2>& default_val) const
  {
    std::array<T,2> result;
    result[0] = get_property<T>(id + ".<xmlattr>.x", default_val[0]);
    result[1] = get_property<T>(id + ".<xmlattr>.y", default_val[1]);
    return result;
  }
  
  template<class T>
  util::vector2 get_vector2_property(const config_node& v,
                           const std::string& id,
                           const std::array<T,2>& default_val) const
  {
    util::vector2 result;
    result[0] = v.second.get<util::scalar>(id + ".<xmlattr>.x", default_val[0]);
    result[1] = v.second.get<util::scalar>(id + ".<xmlattr>.y", default_val[1]);
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

  util::scalar _shear;
  util::scalar _sigma_smooth;
  util::scalar _shear_rotation_angle;
  boost::mpi::communicator _comm;
  int _master_rank;
  
  std::string _fits_output;
  std::string _raw_output;
  
  method_type _method;
  lens_plane_type _lens_plane_type;
  
  std::vector<post_processing_step_descriptor> _post_processing_steps;
  
  friend class boost::serialization::access;

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & _resolution;
    ar & _screen_pos;
    ar & _screen_size;
    ar & _shear;
    ar & _sigma_smooth;
    ar & _star_files;
    ar & _random_star_generators;
    ar & _fits_output;
    ar & _raw_output;
    ar & _method;
    ar & _lens_plane_type;
    ar & _post_processing_steps;
    ar & _tree;
  }
};

}

#endif	/* INPUT_HPP */

