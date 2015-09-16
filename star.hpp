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
    result = this->_position;
    util::sub(result, position);

    util::scalar squared_norm = util::dot(result, result);

    if(squared_norm == 0.0)
    {
      result = {0.0, 0.0};
      return;
    }

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
    _deflection_constant = _mass;
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

template<unsigned Multipole_order>
class pseudo_star
{
public:
  static_assert(Multipole_order <= 6, "No moments higher than 64-pole are supported");
  
  pseudo_star() = default;
  
  pseudo_star(const std::vector<star>& star_list)
  {
    _center_of_mass = {0.0, 0.0};
    _total_mass = 0.0;
    for(std::size_t i = 0; i < Multipole_order - 1; ++i)
    {
      _multipole_moments[i] = {0.0, 0.0};
    }
    
    for(const star& s : star_list)
    {
      _total_mass += s.get_mass();
      
      util::vector2 weighted_position = s.get_position();
      util::scale(weighted_position, s.get_mass());
      
      util::add(_center_of_mass, weighted_position);
    }
    
    if(!star_list.empty())
    {
      util::scale(_center_of_mass, 1.0 / _total_mass);

      for(const star& s : star_list)
      {
        util::vector2 delta = s.get_position();
        util::sub(delta, _center_of_mass);

        std::array<util::scalar, 7> delta_x_power;
        std::array<util::scalar, 7> delta_y_power;
        delta_x_power[0] = 1.;
        delta_y_power[0] = 1.;
        delta_x_power[1] = delta[0];
        delta_y_power[1] = delta[1];

        // For the multipole expansion of the deflection angle,
        // see Wambsganß, J.: Gravitational microlensing., 
        // Max-Planck-Institut für Physik und Astrophysik, Garching (Germany). Inst. für Astrophysik, 
        // Oct 1990
        for(std::size_t i = 2; i < 7; ++i)
        {
          delta_x_power[i] = delta[0] * delta_x_power[i - 1];
          delta_y_power[i] = delta[1] * delta_y_power[i - 1];
        }

        if(Multipole_order >= 2)
        {
          // Quadrupole
          _multipole_moments[0][0] += s.get_mass() * (delta_x_power[2] - delta_y_power[2]);
          _multipole_moments[0][1] += 2 * s.get_mass() * delta_x_power[1] * delta_y_power[1];
        }

        if(Multipole_order >= 3)
        {
          // Octopole
          _multipole_moments[1][0] += s.get_mass() * (delta_x_power[3] - 3. * delta_x_power[1] * delta_y_power[2]);
          _multipole_moments[1][1] += s.get_mass() * (3 * delta_x_power[2] * delta_y_power[1] - delta_y_power[3]);     
        }

        if(Multipole_order >= 4)
        {
          // 16-pole
          _multipole_moments[2][0] += 
            s.get_mass() * (  delta_x_power[4] 
                            - 6. * delta_x_power[2] * delta_y_power[2] 
                            + delta_y_power[4]);

          _multipole_moments[2][1] += 
            s.get_mass() * (  4 * delta_x_power[3] * delta_y_power[1] 
                            - 4 * delta_x_power[1] * delta_y_power[3]);       
        }

        if(Multipole_order >= 5)
        {
          // 32-pole
          _multipole_moments[3][0] += 
            s.get_mass() * (  delta_x_power[5] 
                            - 10. * delta_x_power[3] * delta_y_power[2] 
                             + 5. * delta_x_power[1] * delta_y_power[4]);
          _multipole_moments[3][1] += 
            s.get_mass() * (   5 * delta_x_power[4] * delta_y_power[1] 
                            - 10 * delta_x_power[2] * delta_y_power[3]
                             + delta_y_power[5]);         
        }

        if(Multipole_order >= 6)
        {
          // 64-pole
          _multipole_moments[4][0] +=
            s.get_mass() * (  delta_x_power[6] 
                            - 15 * delta_x_power[4] * delta_y_power[2] 
                            + 15 * delta_x_power[2] * delta_y_power[4] 
                            - delta_y_power[6]);

          _multipole_moments[4][1] +=
            s.get_mass() * (6 * delta_x_power[5] * delta_y_power[1] 
                         - 20 * delta_x_power[3] * delta_y_power[3] 
                          + 6 * delta_x_power[1] * delta_y_power[5]);
        }
      }
    }
  }  
  
  inline void calculate_deflection_angle(const util::vector2& position,
                                        util::vector2& result) const
  {
    //TODO Think about a solution for this case
    assert(position != this->_center_of_mass);
    
    util::vector2 R = position;
    util::sub(R, this->_center_of_mass);

    util::scalar squared_norm = util::dot(R, R);

    // monopole
    result = R;
    util::scale(result, -_total_mass / squared_norm);

    std::array<util::vector2, 8> r_power;
    r_power[0] = {1., 1.};
    r_power[1] = R;
    for(std::size_t i = 2; i < r_power.size(); ++i)
    {
      r_power[i][0] = R[0] * r_power[i - 1][0];
      r_power[i][1] = R[1] * r_power[i - 1][1];
    }
    
    std::array<util::vector2, Multipole_order - 1> multipole_evaluation_coefficients;
    
    if(Multipole_order >= 2)
    {
      multipole_evaluation_coefficients[0][0] = r_power[3][0] - 3 * r_power[1][0] * r_power[2][1];
      multipole_evaluation_coefficients[0][1] = 3 * r_power[2][0]*r_power[1][1] - r_power[3][1];
    }
    
    if(Multipole_order >= 3)
    {
      multipole_evaluation_coefficients[1][0] =
        r_power[4][0] - 6 * r_power[2][0] * r_power[2][1] + r_power[4][1];
      
      multipole_evaluation_coefficients[1][1] =
        4 * r_power[3][0] * r_power[1][1] - 4 * r_power[1][0] * r_power[3][1];
    }
    
    if(Multipole_order >= 4)
    {
      multipole_evaluation_coefficients[2][0] =
        r_power[5][0] - 10 * r_power[3][0] * r_power[2][1] + 5 * r_power[1][0] * r_power[4][1];
      
      multipole_evaluation_coefficients[2][1] =
        5 * r_power[4][0] * r_power[1][1] - 10 * r_power[2][0] * r_power[3][1] + r_power[5][1];
    }
    
    if(Multipole_order >= 5)
    {
      multipole_evaluation_coefficients[3][0] =
        r_power[6][0] - 15 * r_power[4][0] * r_power[2][1] + 15 * r_power[2][0] * r_power[4][1] - r_power[6][1];
      
      multipole_evaluation_coefficients[3][1] =
        6 * r_power[5][0] * r_power[1][1] - 20 * r_power[3][0] * r_power[3][1] + 6 * r_power[1][0] * r_power[5][1];
    }
    
    if(Multipole_order >= 6)
    {
      multipole_evaluation_coefficients[4][0] =
        r_power[7][0] - 21 * r_power[5][0] * r_power[2][1] 
        + 35 * r_power[3][0] * r_power[4][1] - 7 * r_power[1][0] * r_power[6][1];
      
      multipole_evaluation_coefficients[4][1] =
        7 * r_power[6][0] * r_power[1][1] - 35 * r_power[4][0] * r_power[3][1] 
        + 21 * r_power[2][0] * r_power[5][1] - r_power[7][1];
    }
    
    // Add multipole corrections
    util::scalar r_power2i_plus2 = squared_norm * squared_norm * squared_norm;
    for(std::size_t i = 0; i < Multipole_order - 1; ++i)
    {
      util::scalar m0c0 = _multipole_moments[i][0] * multipole_evaluation_coefficients[i][0];
      util::scalar m1c1 = _multipole_moments[i][1] * multipole_evaluation_coefficients[i][1];
      
      util::scalar m1c0 = _multipole_moments[i][1] * multipole_evaluation_coefficients[i][0];
      util::scalar m0c1 = _multipole_moments[i][0] * multipole_evaluation_coefficients[i][1];
      
      result[0] -= (m0c0 + m1c1) / r_power2i_plus2;
      result[1] -= (m0c1 - m1c0) / r_power2i_plus2;

      r_power2i_plus2 *= squared_norm;
    }
    
  }
  
  const util::vector2& center_of_mass() const
  {
    return _center_of_mass;
  }
private:
  util::scalar _total_mass;
  util::vector2 _center_of_mass;

  std::array<util::vector2, Multipole_order - 1> _multipole_moments;

};

}

#endif	/* STAR_HPP */

