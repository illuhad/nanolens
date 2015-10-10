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

#ifndef LENS_PLANE_HPP
#define	LENS_PLANE_HPP

#include <vector>
#include <cmath>
#include "util.hpp"
#include "plane.hpp"
#include "star.hpp"
#include "tree.hpp"
#include "input.hpp"

namespace nanolens
{

/// An implementation of the deflection engine concept that uses a tree to
/// accelerate the calculation of the deflection angle.
class tree_deflector
{
public:
  class settings 
  {
  public:
    /// Construct object
    /// \param config A configuration object from which the deflector will extract
    /// all relevant settings
    explicit settings(const configuration& config)
    {
      _opening_angle = config.get_deflection_engine_property("opening_angle", 1.0);
    }
    
    /// Construct object
    /// \param opening_angle The opening angle of the tree, i.e. the ratio between
    /// cell diameter and distance above which the deflection of a cell can
    /// be approximated
    explicit settings(util::scalar opening_angle = 1.0)
    : _opening_angle(opening_angle) {}
    
    /// @return The opening angle
    util::scalar get_opening_angle() const
    {
      return _opening_angle;
    }
  private:
    util::scalar _opening_angle;
  };
  
  /// Construct object
  /// @param stars A list containing the stars to be put into the tree
  /// @param config A settings object from which all settings will be read
  tree_deflector(const std::vector<star>& stars, const settings& config)
  : _tree(stars, config.get_opening_angle())
  {
  }
  
  /// Calculates the deflection angle
  /// \param position The position in the lens plane at which the deflection
  /// shall be calculated
  /// \param result A two dimensional vector into which the result of the
  /// calculation will be written.
  inline void operator()(const util::vector2& position, util::vector2& result) const
  {
    result = _tree.get_deflection(position);
  }
  
private:
  barnes_hut_tree _tree;
};

/// An implemtation of the deflection engine concept that calculates the deflection
/// of the stars by simple, direct summation
class exact_deflector
{
public:
  /// The definition of a settings object is required by the deflection engine concept,
  /// however the exact deflector does not need any configuration. Hence, the settings
  /// object is empty.
  struct settings 
  {
    explicit settings(const configuration&){}
    settings(){}
  };
  
  /// Construct object
  /// \param stars A list of stars that shall be used to calculate the deflection
  /// angle
  /// \param settings An instance of a settings object (which can be empty for
  /// the exact deflector)
  exact_deflector(const std::vector<star>& stars, const settings&)
  : _stars(stars)
  {}
  
  /// Calculates the deflection angle
  /// \param position The position in the lens plane at which the deflection
  /// shall be calculated
  /// \param result A two dimensional vector into which the result of the
  /// calculation will be written.
  inline void operator()(const util::vector2& position, util::vector2& result) const
  {
    result = {0.0, 0.0};
    
    util::vector2 current_deflection = {0.0, 0.0};
    for(std::size_t i = 0; i < _stars.size(); ++i)
    {
      _stars[i].calculate_deflection_angle(position, current_deflection);
      util::add(result, current_deflection);
    }
  }
  
private:
  std::vector<star> _stars;
};

/// Implementation of a lens plane for microlensing, i.e. the deflection will be
/// due to many individual stars. This implementation uses the normalized lensing
/// equation: \f$ y=Ax-\sigma_cx-\sum_i\alpha_i(x) \f$ with the shear matrix A,
/// the deflection due to the stars \f$\alpha_i\f$, and the smooth matter contribution
/// \f$\sigma_c\f$.
/// \tparam Deflector_type The deflection engine that shall be used, e.g a tree or
/// direct summation of the deflection angles of all individual stars.
template<class Deflector_type>
class microlensing_lens_plane : public plane
{
public:
  typedef typename Deflector_type::settings deflector_settings_type;

  typedef std::vector<star>::iterator star_iterator;
  typedef std::vector<star>::const_iterator const_star_iterator;

  /// Construct object
  /// \param deflectors A list of all individual stars that shall be used. The
  /// positions of the stars are expected to be given in units of Einstein radii.
  /// \param settings The settings object of the selected deflection engine.
  /// This object will be used to initialize the deflection engine, and can
  /// thus be used to pass any necessary configuration to the deflection engine.
  /// \param shear The external shear
  /// \param sigma_smooth The normalized surface density due to smooth matter
  /// \param shear_rotation_angle An angle in degrees which s pecifies the 
  /// direction of the shear. 0 degrees means that the shear will be along the
  /// y-axis.
  explicit microlensing_lens_plane(const std::vector<star>& deflectors,
                      const deflector_settings_type& settings,
                      util::scalar shear = 0.,
                      util::scalar sigma_smooth = 0.,
                      util::scalar shear_rotation_angle = 0.)
  : _deflectors(deflectors), 
  _deflection_engine(deflectors, settings),
  _shear(shear),
  _sigma_smooth(sigma_smooth),
  _shear_rotation_angle(shear_rotation_angle)
  {
    // Get mean mass
    util::scalar mean_mass = 0.0;

    for(const star& s : deflectors)
      mean_mass += s.get_mass();
    mean_mass /= deflectors.size();

    // Estimate sigma*
    // First calculate distribution center
    util::vector2 distribution_center = {0.0,0.0};
    for(const star& s : deflectors)
      util::add(distribution_center, s.get_position());
    util::scale(distribution_center, 1.0 / deflectors.size());

    // Estimate average radius
    util::scalar avg_radius_squared = 0.0;
    for(const star& s : deflectors)
    {
      util::vector2 R = s.get_position();
      util::sub(R, distribution_center);
      avg_radius_squared += util::dot(R,R);
    }
    avg_radius_squared /= deflectors.size();

    util::scalar N_stars = 0.0;
    // Count stars within average radius
    for(const star& s: deflectors)
    {
      util::vector2 R = s.get_position();
      util::sub(R, distribution_center);
      if(util::dot(R,R) <= avg_radius_squared)
        ++N_stars;
    }
    //util::scalar area_radius_squared = 
    //  avg_radius_squared / util::square(calculate_einstein_radius(1.0,1.0,2.0));
    if(!deflectors.empty())
      _sigma_star = N_stars * mean_mass / avg_radius_squared;
    else _sigma_star = 0.0;


    this->_smooth_matter_fraction = _sigma_smooth / (_sigma_star + _sigma_smooth);

    // Set up shear matrix
    
    util::scalar shear_angle_rad = 
      _shear_rotation_angle * 2.0 * boost::math::constants::pi<util::scalar>() / 360.0;
    
    util::scalar cos_phi = std::cos(shear_angle_rad);
    util::scalar sin_phi = std::sin(shear_angle_rad);
    util::scalar cos2_phi = util::square(cos_phi);
    util::scalar sin2_phi = util::square(sin_phi);
    util::scalar A = 1.0 - _shear;
    util::scalar B = 1.0 + _shear;
    
    _shear_matrix = {{{A * cos2_phi + B * sin2_phi, (A - B)* cos_phi * sin_phi}, 
                      {(A -B )* cos_phi * sin_phi, A * sin2_phi + B * cos2_phi}}};

    //_shear_matrix[0][0] -= _sigma_smooth - _sigma_star;
    //_shear_matrix[1][1] -= _sigma_smooth - _sigma_star;

    _shear_matrix[0][0] -= _sigma_smooth;
    _shear_matrix[1][1] -= _sigma_smooth;

  }

  /// Calculates the deflection angle  at a given position in the lens plane.
  /// \param position The position in the lens plane at which the deflection
  /// shall be calculated
  /// \param result A two dimensional vector into which the result of the
  /// calculation will be written.
  inline void get_deflection_angle(const util::vector2& position, util::vector2& result) const
  {
    this->_deflection_engine(position, result);
  }

  /// Estimate the radius of the stellar mass distribution in the plane
  /// @return An estimate for the radius in Einstein radii.
  util::scalar get_radius_estimate() const
  {
    util::scalar result = 0.0;

    for(std::size_t i = 0; i < _deflectors.size(); ++i)
    {
      util::scalar distance = std::sqrt(util::dot(_deflectors[i].get_position(),
                                                  _deflectors[i].get_position()));

      if(distance > result)
        result = distance;
    }

    if(result == 0.0)
      return 1.0;

    return result;
  }

  /// Transforms from a position in the lens plane to a position in the source plane
  /// using the normalized lens equation.
  /// \param lens_plane_pos The position in the lens plane in Einstein radii.
  /// \return The resulting position in the source plane in Einstein radii.
  inline util::vector2 lensing_transformation(const util::vector2& lens_plane_pos) const
  {
    util::vector2 out;
    util::matrix_vector_mult(_shear_matrix, lens_plane_pos, out);

    util::vector2 deflection;
    get_deflection_angle(lens_plane_pos, deflection);
    util::add(out, deflection);

    return out;
  }

  /// Given a position in the source plane and a position in the lens plane,
  /// calculates the corresponding position in the image plane.
  /// \param lens_plane_pos The position in the plane in Einstein radii
  /// \param source_plane_pos The position in the source plane in Einstein radii
  /// \return The position in the image plane
  inline util::vector2 inverse_lensing_transformation(const util::vector2& lens_plane_pos,
                                                  const util::vector2& source_plane_pos) const
  {
    util::vector2 out = source_plane_pos;

    util::vector2 lensing_transformation_result = lensing_transformation(lens_plane_pos);

    util::sub(out, lensing_transformation_result);

    return out;
  }

  /// \return A iterator to the beginning of the star list
  inline star_iterator begin_stars()
  { return _deflectors.begin(); }

  /// \return A iterator to the beginning of the star list
  inline const_star_iterator begin_stars() const
  { return _deflectors.begin(); }

  /// \return A iterator to the end of the star list
  inline star_iterator end_stars()
  { return _deflectors.end(); }

  /// \return A iterator to the end of the star list
  inline const_star_iterator end_stars() const
  { return _deflectors.end(); }

  /// \return The number of stars
  inline std::size_t num_stars() const
  { return _deflectors.size(); }

  /// Calculates the distance to the nearest star from a given position
  /// \param position The position in Einstein radii from which the nearest star
  /// shall be found
  /// \return The squared distance to the nearest star
  util::scalar find_squared_distance_to_nearest_star(const util::vector2& position) const
  {
    auto nearest_star = get_nearest_star(position);

    if(nearest_star == _deflectors.end())
      return std::numeric_limits<util::scalar>::max();

    util::vector2 diff_vector = position;
    util::sub(diff_vector, nearest_star->get_position());

    return util::dot(diff_vector, diff_vector);
  }

  /// Finds the nearest star from a given position
  /// \param position The position in Einstein radii from which the nearest star
  /// shall be found
  /// \return The position of the nearest star in Einstein radii, or a vector
  /// containing the maximum values of util::scalar if there are no stars.
  util::vector2 get_nearest_star_position(const util::vector2& position) const
  {
    auto nearest_star = get_nearest_star(position);
    if(nearest_star != _deflectors.end())
      return nearest_star->get_position();
    else return {std::numeric_limits<util::scalar>::max(),
                 std::numeric_limits<util::scalar>::max()};

  }

  /// \return A star based on its index in the star list
  /// \param idx The index of the star. It is the responsibility of the user
  /// to ensure that \c idx is smaller than the number of stars.
  star& get_star_by_index(std::size_t idx)
  { return _deflectors[idx]; }

  /// \return A star based on its index in the star list
  /// \param idx The index of the star. It is the responsibility of the user
  /// to ensure that \c idx is smaller than the number of stars.
  const star& get_star_by_index(std::size_t idx) const
  { return _deflectors[idx]; }

  /// \return The number of stars
  std::size_t get_num_stars() const
  { return _deflectors.size(); }

  /// \return The external shear
  util::scalar get_shear() const
  {
    return _shear;
  }

  /// \return The fraction between smooth matter density and total density
  util::scalar get_smooth_matter_fraction() const
  {
    return _smooth_matter_fraction;
  }

  /// \return The normalized surface density of the smooth matter.
  util::scalar get_sigma_smooth() const
  {
    return _sigma_smooth;
  }

  /// \return The normalized surface density of the stellar mass. If the
  /// star distribution is not uniform, the mean surface density will be returned.
  util::scalar get_sigma_star() const
  {
    return _sigma_star;
  }

  /// \return The mean total surface density.
  util::scalar get_mean_surface_density() const
  {
    return _sigma_smooth + _sigma_star;
  }

  
  /// Gathers information about the lens plane.
  /// \param out A map where the information will be stored. Each key in the map
  /// will specify the property name and the corresponding value will be the property
  /// value.
  void obtain_properties_set(std::map<std::string, util::scalar>& out) const
  {
    out.clear();

    out["N*"] = this->get_num_stars();
    out["sigma* (estimated)"] = this->get_sigma_star();
    out["sigma smooth"] = this->get_sigma_smooth();
    out["sigma total"] = this->get_mean_surface_density();
    out["shear"] = this->get_shear();
  }
  
  /// Estimates the coordinates of a rectangle in the lens plane that is mapped
  /// to a given rectangle in the source plane
  /// \param source_plane_coordinates A matrix containing the coordinates of the
  /// source plane rectangle
  /// \param out A matrix into which the resulting coordinates in the lens plane
  /// will be written.
  void estimate_mapped_region_coordinates(const util::matrix_nxn<util::vector2, 2>& source_plane_coordinates,
                                        util::matrix_nxn<util::vector2, 2>& out) const
  {
    estimate_mapped_region_coordinates(source_plane_coordinates, out, _sigma_star);
  }
  
  /// Estimates the coordinates of a rectangle in the lens plane that is mapped
  /// to a given rectangle in the source plane. This version uses a user supplied
  /// value for sigma star.
  /// \param source_plane_coordinates A matrix containing the coordinates of the
  /// source plane rectangle
  /// \param out A matrix into which the resulting coordinates in the lens plane
  /// will be written.
  /// \param sigma_star The stellar density that shall be used.
  void estimate_mapped_region_coordinates(const util::matrix_nxn<util::vector2, 2>& source_plane_coordinates,
                                        util::matrix_nxn<util::vector2, 2>& out,
                                        util::scalar sigma_star) const
  {
    // out = (shear_matrix - diag(sigma_total))^-1 * x
    // (The smooth matter term is already contained in the shear matrix)
    util::matrix_nxn<util::scalar, 2> lensing_transformation = _shear_matrix;
    lensing_transformation[0][0] -= sigma_star;
    lensing_transformation[1][1] -= sigma_star;
    
    util::matrix_nxn<util::scalar, 2> inverted_lensing_transformation;
    numeric::invert_matrix(lensing_transformation, inverted_lensing_transformation);
    
    for(std::size_t i = 0; i < 2; ++i)
      for(std::size_t j = 0; j < 2; ++j)
        util::matrix_vector_mult(inverted_lensing_transformation, source_plane_coordinates[i][j],
                                 out[i][j]);

  }
private:
  
  /// Finds the nearest star from a given position
  /// \param position The position from which the nearest star shall be found
  /// \return An iterator to the star or an end-vector if there are no stars.
  std::vector<star>::const_iterator get_nearest_star(const util::vector2& position) const
  {
    auto nearest_star = _deflectors.end();
    util::scalar min_distance = std::numeric_limits<util::scalar>::max();

    for(auto it = _deflectors.begin(); it != _deflectors.end(); ++it)
    {
      util::vector2 diff_vector = position;
      util::sub(diff_vector, it->get_position());

      util::scalar dist = util::dot(diff_vector, diff_vector);

      if(dist < min_distance && dist != 0.0)
      {
        min_distance = dist;
        nearest_star = it;
      }
    }

    return nearest_star;
  }

  /// Finds the nearest star from a given position (index version)
  /// \param position The position from which the nearest star shall be found
  /// \return The index of the nearest star in the star list, or \c get_num_stars()
  /// if there are no stars.
  std::size_t get_nearest_star_index(const util::vector2& position) const
  {
    std::size_t idx = _deflectors.size();
    util::scalar min_distance = std::numeric_limits<util::scalar>::max();

    for(std::size_t i = 0; i < _deflectors.size(); ++i)
    {
      util::vector2 diff_vector = position;
      util::sub(diff_vector, _deflectors[i].get_position());

      util::scalar dist = util::dot(diff_vector, diff_vector);

      if(dist < min_distance && dist != 0.0)
      {
        min_distance = dist;
        idx = i;
      }
    }

    return idx;
  }

  std::vector<star> _deflectors;

  util::scalar _shear;
  util::scalar _smooth_matter_fraction;

  util::scalar _sigma_star;
  util::scalar _sigma_smooth;

  util::scalar _shear_rotation_angle;
  
  util::matrix_nxn<util::scalar, 2> _shear_matrix;

  Deflector_type _deflection_engine;
};
}

#endif	/* LENS_PLANE_HPP */

