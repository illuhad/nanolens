/* 
 * File:   system.hpp
 * Author: aksel
 *
 * Created on 7. April 2015, 05:14
 */

#ifndef SYSTEM_HPP
#define	SYSTEM_HPP

#include <vector>
#include <memory>
#include <fstream>
#include <random>
#include <boost/mpi.hpp>
#include "util.hpp"
#include "ray.hpp"
#include "lens_plane.hpp"
#include "observer_plane.hpp"

namespace nanolens{

class system
{
public:
  static const std::size_t num_planes = 2;

  explicit system(const std::vector<lens_plane::star>& deflectors,
                  const std::array<util::scalar, num_planes>& plane_distances = {1.0, 1.0})
  : _deflector(new lens_plane(deflectors, plane_distances[1])),
    _observer(new observer_plane({0.0, 0.0}, plane_distances[0]))
  {
    init_einstein_radius();
  }


  explicit system(const std::string& star_file_,
                  const std::array<util::scalar, num_planes>& plane_distances = {1.0, 1.0})
  : _observer(new observer_plane({0.0, 0.0}, plane_distances[1]))
  {
    init_einstein_radius(plane_distances[1],
                         plane_distances[0],
                         plane_distances[0] + plane_distances[1]);

    util::scalar lens_plane_einstein_radius = get_einstein_radius() 
            * _observer->distance_to_previous_plane();

    std::ifstream input_file;

    input_file.open(star_file_.c_str());

    std::vector<lens_plane::star> stars;

    if(input_file.is_open())
    {
      while(input_file.good())
      {
        util::vector2 position = {0.0, 0.0};
        util::scalar mass = 0.0;

        input_file >> position[0];
        input_file >> position[1];
        input_file >> mass;

        util::scale(position, lens_plane_einstein_radius);

        lens_plane::star new_star(position, mass);
        stars.push_back(new_star);
      }
    }

    _deflector = std::shared_ptr<lens_plane>(new lens_plane(stars, plane_distances[0]));


  }

  inline const observer_plane& get_observer() const
  { return *_observer; }

  inline const lens_plane& get_deflector() const
  { return *_deflector; }

  inline util::scalar get_distance_from_source_to_observer() const
  {
    return _deflector->distance_to_previous_plane()
            + _observer->distance_to_previous_plane();
  }

  template<typename RayType>
  inline void traverse(RayType& ray) const
  {
    ray.propagate(*_deflector);
    ray.propagate(*_observer);
  }

  inline util::scalar get_einstein_radius() const
  {
    return _einstein_radius;
  }

  system get_empty_system() const
  {
    return system(std::vector<lens_plane::star>(), 
      {_deflector->distance_to_previous_plane(),_observer->distance_to_previous_plane()});
  }
private:
  void init_einstein_radius(util::scalar d_ls, util::scalar d_l, util::scalar d_s)
  {
    _einstein_radius = lens_plane::calculate_einstein_radius(d_ls, d_l, d_s);
  }

  void init_einstein_radius()
  {
    init_einstein_radius(_deflector->distance_to_previous_plane(),
            _observer->distance_to_previous_plane(),
            get_distance_from_source_to_observer());
  }

  std::shared_ptr<observer_plane> _observer;
  std::shared_ptr<lens_plane> _deflector;

  util::scalar _einstein_radius;
};

}

#endif	/* SYSTEM_HPP */

