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

#ifndef TREE_HPP
#define	TREE_HPP

#include "util.hpp"
#include "grid.hpp"
#include "star.hpp"
#include "numeric.hpp"
#include "input.hpp"
#include "interpolator.hpp"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>

#include <boost/geometry/index/rtree.hpp>

namespace nanolens{


namespace tree_impl_{

typedef unsigned cell_id;

template<class Cell_type>
class cell_db
{
public:
  cell_db()
  {
  }
  
  cell_id create_cell()
  {
    _db.push_back(Cell_type());
    return _db.size() - 1;
  }

  inline const Cell_type& operator[](const cell_id id) const
  { return _db[id]; }
  
  inline Cell_type& operator[](const cell_id id)
  { return _db[id]; }
  
  inline Cell_type* data()
  { return _db.data(); }
  
  inline const Cell_type* data() const
  { return _db.data(); }
  
private:
  std::vector<Cell_type> _db;
};

template<unsigned Multipole_order>
class barnes_hut_grid_cell
{
public:
  static const std::size_t max_allowed_stars_per_cell = 1;
  
  static const std::size_t N_subcells_per_dim = 2;
  
  typedef util::matrix_nxn<cell_id, N_subcells_per_dim> grid_layer_type;
  
  barnes_hut_grid_cell()
  : _min_extent({0.0,0.0}), _max_extent({0.0,0.0}) {}

  barnes_hut_grid_cell(const util::vector2& min_extent,
            const util::vector2& max_extent,
            cell_db<barnes_hut_grid_cell<Multipole_order>>* db)
  : _num_exact_stars(0), _db(db), _is_leaf(true)
  {
    set_cell_size(min_extent, max_extent);
  }
  
  class this_descriptor
  {
    cell_id _this_id;
    cell_db<barnes_hut_grid_cell<Multipole_order>>* _db;
  public:
    this_descriptor(cell_id id, cell_db<barnes_hut_grid_cell<Multipole_order>>* db)
    : _this_id(id), _db(db)
    {}
    
    cell_db<barnes_hut_grid_cell<Multipole_order>>* get_db() const
    { return _db; }
    
    barnes_hut_grid_cell<Multipole_order>* get_this()
    { return _db->data() + _this_id; }
    
    cell_id get_cell_id() const
    { return _this_id; }
  };

  inline const util::vector2& get_cell_center() const
  { return _center; }
  
  // only defined if this cell is not a leaf
  inline const util::vector2& get_cell_center_of_mass() const
  { return _approx_pseudo_star.center_of_mass(); }
  
  inline util::scalar get_cell_diameter() const
  { return _diameter; }
  
  inline bool is_leaf() const
  { return _is_leaf; }
  
  inline const grid_layer_type& get_sub_cells() const
  {
    return _sub_grid;
  }
  
  const std::vector<star>& get_interior_stars() const
  {
    return _interior_stars;
  }
  
  const util::vector2& get_cell_size() const
  {
    return _size;
  }
  
  unsigned get_num_directly_computed_stars() const
  {
    return _num_exact_stars;
  }
  
  const std::array<star, max_allowed_stars_per_cell>& get_directly_computed_stars() const
  {
    return _exact_stars;
  }
  
  const pseudo_star<Multipole_order>& get_pseudo_star() const
  {
    return _approx_pseudo_star;
  }
  
  const util::vector2& get_min_extent() const
  {
    return _min_extent;
  }
  
  const util::vector2& get_max_extent() const
  {
    return _max_extent;
  }
  
  inline bool cell_contains_point(const util::vector2& pos) const
  {
    for(std::size_t i = 0; i < N_subcells_per_dim; ++i)
      if(!(_min_extent[i] <= pos[i] && pos[i] < _max_extent[i]))
        return false;
    return true;
  }

  static void add_star(this_descriptor td, const star& s)
  {
    barnes_hut_grid_cell<Multipole_order>* this_ptr = td.get_this();

    assert(this_ptr->cell_contains_point(s.get_position()));
    
    std::array<star, max_allowed_stars_per_cell> existing_stars;

    for(unsigned i = 0;
        i < max_allowed_stars_per_cell;
        ++i)
    {
      if(i < this_ptr->_interior_stars.size())   
        existing_stars[i] = this_ptr->_interior_stars[i];
    }
    this_ptr->_interior_stars.push_back(s);

    if(this_ptr->_interior_stars.size() > max_allowed_stars_per_cell)
    {
      // Subdivide cell
      if(this_ptr->_is_leaf)
      {
        this_ptr->_is_leaf = false;

        
        for(std::size_t x = 0; x < N_subcells_per_dim; ++x)
          for(std::size_t y = 0; y < N_subcells_per_dim; ++y)
          {
            util::vector2 bucket_size = this_ptr->_size;
            util::scale(bucket_size, 1.0 / static_cast<util::scalar>(N_subcells_per_dim));
            
            util::vector2 bucket_min_extent = this_ptr->_min_extent;
            bucket_min_extent[0] += x * bucket_size[0];
            bucket_min_extent[1] += y * bucket_size[1];
            
            util::vector2 bucket_max_extent = bucket_min_extent;
            util::add(bucket_max_extent, bucket_size);
            
            cell_id new_cell_id = td.get_db()->create_cell();
            // Update this_ptr in case the cell_db had to relocate memory
            this_ptr = td.get_this();
            
            this_ptr->_sub_grid[x][y] = new_cell_id;
            (*td.get_db())[new_cell_id] = barnes_hut_grid_cell(bucket_min_extent, bucket_max_extent, td.get_db());

          }

        // Redistribute the one existing star
        for(std::size_t i  = 0; i < max_allowed_stars_per_cell; ++i)
        {
          insert_star_into_subcell(td, existing_stars[i]);
          this_ptr = td.get_this();
        }

      }

      insert_star_into_subcell(td, s);
    }
  }

  void init_pseudo_stars()
  {
    _approx_pseudo_star = pseudo_star<Multipole_order>(this->_interior_stars);
    _center = _approx_pseudo_star.center_of_mass();
    

    if(!_is_leaf)
    {
      for(unsigned x = 0; x < N_subcells_per_dim; ++x)
        for(unsigned y = 0;  y < N_subcells_per_dim; ++y)
          (*_db)[_sub_grid[x][y]].init_pseudo_stars();
      
    }
    
    if(_is_leaf)
      _num_exact_stars = _interior_stars.size();
    
    for(unsigned i = 0; i < _num_exact_stars; ++i)
      _exact_stars[i] = _interior_stars[i];
    
    _interior_stars.clear();
  }
  
private:
  util::vector2 _min_extent;
  util::vector2 _max_extent;
  util::vector2 _size;
  util::scalar _diameter;
  util::vector2 _center;

  grid_layer_type _sub_grid;
  bool _is_leaf;
  std::vector<star> _interior_stars;
  
  std::array<star, max_allowed_stars_per_cell> _exact_stars;
  unsigned _num_exact_stars;
  pseudo_star<Multipole_order> _approx_pseudo_star;
  
  cell_db<barnes_hut_grid_cell<Multipole_order>>* _db;
  
  inline void set_cell_size(const util::vector2& min_extent,
                            const util::vector2& max_extent)
  {
    _min_extent = min_extent;
    _max_extent = max_extent;
    
    _center = min_extent;
    util::add(_center, max_extent);
    util::scale(_center, 0.5);
    
    _diameter = 0.0;
    for(std::size_t i = 0; i < max_extent.size(); ++i)
      _diameter += max_extent[i] - min_extent[i];
    
    _diameter /= static_cast<util::scalar>(max_extent.size());
    
    _size = _max_extent;
    util::sub(_size, _min_extent);
  }

  inline static void insert_star_into_subcell(this_descriptor td, const star& s)
  {
    const util::vector2& pos = s.get_position();
    barnes_hut_grid_cell<Multipole_order>* this_ptr = td.get_this();

    for(unsigned x = 0; x < N_subcells_per_dim; ++x)
      for(unsigned y = 0; y < N_subcells_per_dim; ++y)
      {
        if((*td.get_db())[this_ptr->_sub_grid[x][y]].cell_contains_point(pos))
        {
          add_star(
              this_descriptor(this_ptr->_sub_grid[x][y], td.get_db()),
              s);
          return;
        }
      }
  }
  

};
}

class barnes_hut_tree
{
public:
  typedef tree_impl_::barnes_hut_grid_cell<6> grid_cell;
  
  barnes_hut_tree()
  : _root_descriptor(0, nullptr){}
  
  barnes_hut_tree(const std::vector<star>& stars,
             util::scalar accuracy)
  : _accuracy_squared(util::square(accuracy)),
    _root_descriptor(0, nullptr)
  {

    util::vector2 min_extent;
    util::vector2 max_extent;
    get_extents_of_star_distribution(stars, min_extent, max_extent);
    
    tree_impl_::cell_id root_id = _cell_db.create_cell();
    _root_descriptor = grid_cell::this_descriptor(root_id, &_cell_db);
    _cell_db[root_id] = grid_cell(min_extent, max_extent, &_cell_db);
    _tree_root = _root_descriptor.get_this();

    for(const star& s : stars)
    {
      grid_cell::add_star(_root_descriptor, s);
      _tree_root = _root_descriptor.get_this();
    }
    _tree_root->init_pseudo_stars();
    
    
    util::vector2 cell_size = {0.5, 0.5};
    
    auto interpolation_cell_constructor = 
    [this, cell_size](interpolation_cell& entry, util::vector2 cell_center)
    {
      init_cache_grid_entry(entry, cell_center, cell_size);
    };
    
    
    util::vector2 cache_grid_min_extent = min_extent;
    util::vector2 cache_grid_max_extent = max_extent;
    
    std::array<std::size_t, 2> num_cells = 
          {static_cast<std::size_t>((cache_grid_max_extent[0] - cache_grid_min_extent[0]) / cell_size[0]),
           static_cast<std::size_t>((cache_grid_max_extent[1] - cache_grid_min_extent[1]) / cell_size[1])};
    
    for(std::size_t i = 0; i < num_cells.size(); ++i)
    if(num_cells[i] == 0)
      num_cells[i] = 1;

    _interpolation_cell_cache = std::shared_ptr<interpolation_cell_cache_type>(
                  new interpolation_cell_cache_type(cache_grid_min_extent, 
                                                    cache_grid_max_extent, 
                                                    num_cells, 
                                                    interpolation_cell_constructor,
                                                    2 * num_cells[1]));
    
  }
  
  
  inline util::vector2 get_deflection(const util::vector2& pos) const
  {
    return get_deflection_from_grid(pos);
  }
  
private:
    
  struct interpolation_cell
  {
    typedef standard_vector2_interpolator interpolator_type;
    //typedef standard_vector2_interpolator interpolator_type;
    
    static const std::size_t num_interpolators = 2;
    
    std::array<std::array<interpolator_type, num_interpolators>, num_interpolators> far_cell_interpolators;
    
    std::vector<const star*> close_lenses;
    
    typedef util::grid_translator<util::scalar,2,false> interpolator_grid_type;
    
    interpolator_grid_type far_cell_interpolator_grid;
    
    interpolation_cell() = default;
  
  };
  
  
  inline util::vector2 get_deflection_from_grid(const util::vector2& pos) const
  {
    const interpolation_cell& grid_entry = (*_interpolation_cell_cache)[pos];
    
    util::vector2 result = {0.0, 0.0};
    
    for(unsigned i = 0; i < grid_entry.close_lenses.size(); ++i)
    {
      util::vector2 contribution;
      grid_entry.close_lenses[i]->calculate_deflection_angle(pos, contribution);
      util::add(result, contribution);
    }
    
    interpolation_cell::interpolator_grid_type::index_type interpolator_index
      = grid_entry.far_cell_interpolator_grid(pos);
    
    util::vector2 interpolated_contribution
      = grid_entry.far_cell_interpolators[interpolator_index[0]][interpolator_index[1]].interpolate(pos);
    util::add(result, interpolated_contribution);
    
    return result;
  }
  
  inline void init_cache_grid_entry(interpolation_cell& entry, 
                                    const util::vector2& cell_center,
                                    const util::vector2& cell_size) const
  {

    std::vector<tree_impl_::cell_id> far_cells;
    this->find_contributing_cells(cell_center,
                           this->_root_descriptor.get_cell_id(),
                            entry.close_lenses,
                            far_cells);


    util::vector2 half_bucket_size = cell_size;
    util::scale(half_bucket_size, 0.5);

//      entry.far_cell_interpolator
//                = interpolation_cell::interpolator_type(cell_center,
//                                                        0.01);

    util::vector2 bucket_min_extent = cell_center;
    util::sub(bucket_min_extent, half_bucket_size);

    util::vector2 bucket_max_extent = cell_center;
    util::add(bucket_max_extent, half_bucket_size);

    std::array<std::size_t, 2> num_cells = {interpolation_cell::num_interpolators,
                                            interpolation_cell::num_interpolators};

    entry.far_cell_interpolator_grid = interpolation_cell::interpolator_grid_type(bucket_min_extent,
                                                                                  bucket_max_extent,
                                                                                  num_cells);

    util::vector2 interpolation_cell_size = entry.far_cell_interpolator_grid.get_bucket_size();

    util::vector2 interpolation_cell_min_extent = bucket_min_extent;
    util::vector2 interpolation_cell_max_extent = bucket_min_extent;
    util::add(interpolation_cell_max_extent, interpolation_cell_size);

    for(std::size_t i = 0; i < interpolation_cell::num_interpolators; ++i)
    {
      interpolation_cell_min_extent[1] = bucket_min_extent[1];
      interpolation_cell_max_extent[1] = bucket_min_extent[1] + interpolation_cell_size[1]; 
      for(std::size_t j = 0; j < interpolation_cell::num_interpolators; ++j)
      {

        entry.far_cell_interpolators[i][j]
            = interpolation_cell::interpolator_type(interpolation_cell_min_extent,
                                                    interpolation_cell_max_extent);

        entry.far_cell_interpolators[i][j].init(
        [&far_cells, this](const util::vector2& pos) -> util::vector2
        {
          util::vector2 deflection = {0.0, 0.0};
          for(tree_impl_::cell_id id : far_cells)
          {
            util::vector2 contribution;
            this->_cell_db[id].get_pseudo_star().calculate_deflection_angle(pos, contribution);
            util::add(deflection, contribution);
          }
          return deflection;
        });
        
        interpolation_cell_min_extent[1] += interpolation_cell_size[1];
        interpolation_cell_max_extent[1] += interpolation_cell_size[1];
      }
      interpolation_cell_min_extent[0] += interpolation_cell_size[0];
      interpolation_cell_max_extent[0] += interpolation_cell_size[0];
    }
  }
  
  void get_extents_of_star_distribution(const std::vector<star>& stars,
                                        util::vector2& min_extent,
                                        util::vector2& max_extent) const
  {
    min_extent = {0.0, 0.0};
    max_extent = {0.0, 0.0};
    
    if(!stars.empty())
    {
      util::scalar scalar_max = std::numeric_limits<util::scalar>::max();
      util::scalar scalar_min = std::numeric_limits<util::scalar>::lowest();      
      
      min_extent = {scalar_max, scalar_max};
      max_extent = {scalar_min, scalar_min};
      
      for(const star& s : stars)
      {
        const util::vector2& pos = s.get_position();
        for(std::size_t i = 0; i < 2; ++i)
        {
          if(pos[i] < min_extent[i])
            min_extent[i] = pos[i];
          
          if(pos[i] > max_extent[i])
            max_extent[i] = pos[i];
        }
      }
    }
    
    // Find a square from min and max extent
    util::vector2 center = min_extent;
    util::add(center, max_extent);
    util::scale(center, 0.5);
    
    util::vector2 side_lengths = max_extent;
    util::sub(side_lengths, min_extent);

    // add 0.02 to the side lengths to avoid having stars at the border
    util::scalar max_side_length = std::max(side_lengths[0], side_lengths[1]) + 0.02;
          
    for(std::size_t i = 0; i < 2; ++i)
    {
      max_extent[i] = center[i] + 0.5 * max_side_length;
      min_extent[i] = center[i] - 0.5 * max_side_length;
    }
  }
  
  inline util::vector2 get_deflection_of_cell_by_tree_traversal(const grid_cell& cell, 
                                                                const util::vector2& pos) const
  {
    util::vector2 result = {0.0, 0.0};
    
    if(cell.is_leaf())
    {
      if(cell.get_num_directly_computed_stars() > 0)
      {
        assert(cell.get_num_directly_computed_stars() <= grid_cell::max_allowed_stars_per_cell);
        static_assert(grid_cell::max_allowed_stars_per_cell == 1, 
                      "Currently only one star per leaf cell is supported");

        cell.get_directly_computed_stars()[0].calculate_deflection_angle(pos, result);
      }
      return result;
    }
    else
    {
      // Check if we can use the approximation
      util::vector2 diff = pos;
      util::sub(diff, cell.get_cell_center_of_mass());
      
      util::scalar distance_squared = util::dot(diff, diff);
      
      if(util::square(cell.get_cell_diameter()) / distance_squared < _accuracy_squared)
      {
        // Get approximated deflection
        cell.get_pseudo_star().calculate_deflection_angle(pos, result);
      }
      else
      {
        for(unsigned x = 0; x < grid_cell::N_subcells_per_dim; ++x)
          for(unsigned y = 0;  y < grid_cell::N_subcells_per_dim; ++y)
          {
            util::vector2 contribution = 
              get_deflection_of_cell_by_tree_traversal(_cell_db[cell.get_sub_cells()[x][y]], pos);
            
            util::add(result, contribution);
          }

      }
    }
    
    return result;
  }
  
  void print_cell_coordinates(std::ostream& ostr, tree_impl_::cell_id cell)
  {
    std::array<util::vector2, 4> corners;
    
    util::vector2 min_extent = _cell_db[cell].get_min_extent();
    util::vector2 max_extent = _cell_db[cell].get_max_extent();
    util::vector2 cell_size = _cell_db[cell].get_cell_size();
    
    corners[0] = min_extent;
    corners[1] = {min_extent[0] + cell_size[0], min_extent[1]};
    corners[2] = max_extent;
    corners[3] = {min_extent[0], min_extent[1] + cell_size[1]};
    
    ostr << "\n";
    for(std::size_t i = 0; i < 4; ++i)
    {
      ostr << corners[i][0] << " " << corners[i][1] << std::endl;
    }
    ostr << "\n";
    
    if(!_cell_db[cell].is_leaf())
      for(unsigned x = 0; x < 2; ++x)
        for(unsigned y = 0; y < 2; ++y)
          print_cell_coordinates(ostr, _cell_db[cell].get_sub_cells()[x][y]);
    
  }

  typedef util::buffered_grid_cache<util::scalar, 2, interpolation_cell> interpolation_cell_cache_type;
  
  void find_contributing_cells(const util::vector2& pos,
                               tree_impl_::cell_id cell_idx,
                               std::vector<const star*>& close_lenses,
                               std::vector<tree_impl_::cell_id>& far_cells) const
  {
    const grid_cell& cell = _cell_db[cell_idx];
    
    if(cell.is_leaf())
    {
      for(unsigned i = 0; i < cell.get_num_directly_computed_stars(); ++i)
        close_lenses.push_back(&(cell.get_directly_computed_stars()[i]));
    }
    else
    {
      // Check if we can use the approximation
      util::vector2 diff = pos;
      util::sub(diff, cell.get_cell_center_of_mass());
      
      util::scalar distance_squared = util::dot(diff, diff);
      
      if(util::square(cell.get_cell_diameter()) / distance_squared < _accuracy_squared)
      {
        far_cells.push_back(cell_idx);
      }
      else
      {
        for(unsigned x = 0; x < grid_cell::N_subcells_per_dim; ++x)
          for(unsigned y = 0;  y < grid_cell::N_subcells_per_dim; ++y)
          {
            find_contributing_cells(pos, 
                                    cell.get_sub_cells()[x][y], 
                                    close_lenses,
                                    far_cells);
          }

      }
    }
  }
  
  
  mutable std::shared_ptr<interpolation_cell_cache_type> _interpolation_cell_cache;
  //preevaluation_grid_type _tree_preevaluation_grid;
  grid_cell::this_descriptor _root_descriptor;
  grid_cell* _tree_root;
  tree_impl_::cell_db<grid_cell> _cell_db;
  util::scalar _accuracy_squared;
};




}

#endif

