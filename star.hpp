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

#ifndef STAR_HPP
#define	STAR_HPP

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>
#include <boost/mpi.hpp>

#include "util.hpp"
#include "numeric.hpp"

namespace nanolens {

class star
{
public:
  star() = default;
  
  star(const util::vector2& pos, util::scalar mass)
  : _mass(mass), _position(pos), _deflection_constant(0.)
  {
    init_deflection_constant();
  }

  inline void calculate_deflection_angle(const util::vector2& position,
                                        util::vector2& result) const
  {
    //TODO Think about a solution for this case
    assert(position != _position);
    
    result = this->_position;
    util::sub(result, position);

    util::scalar squared_norm = util::dot(result, result);


    util::scalar factor = _deflection_constant / squared_norm;

    util::scale(result, factor);
  }

  const util::vector2& get_position() const
  { return _position; }

  util::scalar get_mass() const
  { return _mass; }
  
private:
  util::vector2 _position;
  util::scalar _mass;
  util::scalar _deflection_constant;
  
  inline void init_deflection_constant() 
  {
    _deflection_constant = 4 * util::G * _mass / util::square(util::c);
  }
  
  friend class boost::serialization::access;
  

  template<class Archive>
  void save(Archive& ar, const unsigned int version) const
  {
    ar & _position;
    ar & _mass;
  }

  template<class Archive>
  void load(Archive& ar, const unsigned int version)
  {
    ar & _position;
    ar & _mass;

    init_deflection_constant();
    
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER()
};


namespace impl_{
template<unsigned Multipole_order>
struct pseudo_star_data
{};

template<>
struct pseudo_star_data<1>
{
  
  util::scalar defl_constant;
  util::scalar total_mass;
  util::vector2 monopole_position;
  
  template<unsigned Init_multipole_order>
  static void init(pseudo_star_data<Init_multipole_order>& e,
                   const std::vector<star>& star_list)
  {
    static_assert(Init_multipole_order >= 1, "");
    
    e.monopole_position = {0.0, 0.0};
    e.total_mass = 0.0;
    for(const star& s : star_list)
    {
      util::scale_add(e.monopole_position, s.get_position(), s.get_mass());
      e.total_mass += s.get_mass();
    }
    util::scale(e.monopole_position, 1.0 / e.total_mass);
    
    e.defl_constant = 4 * util::G / util::square(util::c);  
  }
  
  template<class Expansion_type>
  static util::scalar get_squared_expansion(const Expansion_type& e,
                                      const util::vector2& R,
                                      util::scalar R_squared) // R is vector from cms
  {
    return e.total_mass / R_squared;
  }
};

template<>
struct pseudo_star_data<2>
{
  util::scalar defl_constant;
  util::scalar total_mass;
  util::vector2 monopole_position;
  util::vector2 dipole;
  
  template<unsigned Init_multipole_order>
  static void init(pseudo_star_data<Init_multipole_order>& e,
                   const std::vector<star>& star_list)
  {
    static_assert(Init_multipole_order >= 2, "");
    
    pseudo_star_data<1>::init(e, star_list);
    
    e.dipole = {0.0, 0.0};
    for(const star& s : star_list)
    {
      util::vector2 contribution = s.get_position();
      util::sub(contribution, e.monopole_position);
      util::scale(contribution, s.get_mass());
      util::add(e.dipole, contribution);
    }
  }
  
  template<class Expansion_type>
  static util::scalar get_squared_expansion(const Expansion_type& e,
                                      const util::vector2& R,
                                      util::scalar R_squared) // R is vector from cms
  {
    util::scalar result = pseudo_star_data<1>::get_squared_expansion(e, R, R_squared);
    
    result += 1.0 / util::square(R_squared) * util::dot(e.dipole, R);
    return result;
  }
};

template<>
struct pseudo_star_data<3>
{
  util::scalar defl_constant;
  util::scalar total_mass;
  util::vector2 monopole_position;
  util::vector2 dipole;
  util::matrix_nxn<util::scalar, 2> quadrupole;
  
  template<unsigned Init_multipole_order>
  static void init(pseudo_star_data<Init_multipole_order>& e,
                   const std::vector<star>& star_list)
  {
    static_assert(Init_multipole_order >= 3, "");
    
    pseudo_star_data<2>::init(e, star_list);
    
    for(std::size_t i = 0; i < 2; ++i)
      for(std::size_t j = 0; j < 2; ++j)
      {
        e.quadrupole[i][j] = 0.0;
        for(const star& s : star_list)
        {
          const util::vector2& r = s.get_position();
          util::vector2 R = r;
          util::sub(R, e.monopole_position);

          e.quadrupole[i][j] += s.get_mass() * (3. * R[i] * R[j] - util::dot(R,R) * numeric::kronecker_delta(i,j));
        }
      }
  }
  
  template<class Expansion_type>
  static util::scalar get_squared_expansion(const Expansion_type& e,
                                      const util::vector2& R,
                                      util::scalar R_squared) // R is vector from cms
  {
    util::scalar result = pseudo_star_data<2>::get_squared_expansion(e, R, R_squared);
    
    util::scalar quadrupole_factor = 1.0 / (6.0 * util::square(R_squared) * R_squared);
    for(std::size_t i = 0; i < 2; ++i)
      for(std::size_t j = 0; j < 2; ++j)
      {
        result += quadrupole_factor * 
          e.quadrupole[i][j] * (3. * R[i] * R[j] - numeric::kronecker_delta(i,j) * R_squared);
      }
    return result;
  }
};

}

template<unsigned Multipole_order>
class pseudo_star
{
public:
  static_assert(Multipole_order > 0, "Order must be at least 1 for monopole");
  static_assert(Multipole_order <= 3, "No moments higher than quadrupole are supported");
  
  pseudo_star() = default;
  
  pseudo_star(const std::vector<star>& star_list)
  {
    impl_::pseudo_star_data<Multipole_order>::init(this->_data, star_list);
  }  
  
  inline void calculate_deflection_angle(const util::vector2& position,
                                        util::vector2& result) const
  {
    //TODO Think about a solution for this case
    assert(position != _data.monopole_position);
    
    util::vector2 R = this->_data.monopole_position;
    util::sub(R, position);

    util::scalar squared_norm = util::dot(R, R);

    util::vector2 defl_direction = R;
    util::scale(defl_direction, _data.defl_constant);
    
    util::scalar multipole_expansion = 
      impl_::pseudo_star_data<Multipole_order>::get_squared_expansion(_data, R, squared_norm);
      
    util::scale(defl_direction, multipole_expansion);
    
    result = defl_direction;
  }
  
  const util::vector2& center_of_mass() const
  {
    return _data.monopole_position;
  }
private:
 
  
  impl_::pseudo_star_data<Multipole_order> _data;

};

}

#endif	/* STAR_HPP */

