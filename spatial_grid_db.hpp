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

#ifndef SPATIAL_GRID_DB_HPP
#define	SPATIAL_GRID_DB_HPP

#include <boost/serialization/serialization.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/serialization/unordered_map.hpp>

#include <fstream>
#include "star.hpp"
#include "grid.hpp"
#include <unordered_map>

namespace std {

template<class T, std::size_t N>
struct hash<std::array<T, N>>
{
  typedef std::array<T, N> argument_type;
  typedef std::size_t result_type;

  result_type operator()(argument_type const& s) const
  {
    result_type r = 0;
    for(std::size_t i = 0; i < N; ++i)
    {
      r ^= std::hash<T>()(s[i]);
    }
    return r;
  }
};

}
 

namespace nanolens{

template<class T, std::size_t Dimension>
class spatial_grid_db
{
public:
  typedef std::array<long long int, Dimension> grid_index_type;
  typedef std::array<util::scalar, Dimension> point_type;
  typedef std::vector<std::pair<point_type, T>> query_result_list_type;
  
  class storage_map_entry
  {
  public:
    
    const std::string& get_filename() const
    {
      return _filename;
    }
    
    void set_filename(const std::string& name)
    {
      _filename = name;
    }
    
    const std::shared_ptr<std::fstream>& get_file_handle() const
    {
      return _file_handle;
    }
    
    std::shared_ptr<std::fstream> get_file_handle()
    {
      return _file_handle;
    }
    
    void set_file_handle(const std::shared_ptr<std::fstream>& handle)
    {
      _file_handle = handle;
      _serializer = std::shared_ptr<boost::archive::binary_oarchive>(
          new boost::archive::binary_oarchive(*_file_handle));
    }
    
    const std::shared_ptr<boost::archive::binary_oarchive>& get_serializer() const
    {
      return _serializer;
    }
    
    std::shared_ptr<boost::archive::binary_oarchive> get_serializer()
    {
      return _serializer;
    }
    
  private:
    std::string _filename;
    
    // Filehandles are only for write access, and exist only on the master process
    std::shared_ptr<std::fstream> _file_handle;
    std::shared_ptr<boost::archive::binary_oarchive> _serializer;
    
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive& ar, const unsigned int version) const
    {
      ar & _filename;
    }

    template<class Archive>
    void load(Archive& ar, const unsigned int version)
    {
      ar & _filename;
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
  };
  
  typedef std::unordered_map<grid_index_type, storage_map_entry> storage_map_type;
  
  spatial_grid_db(const std::string& db_name,
                  point_type cell_size,
                  std::size_t expected_num_cells,
                  const boost::mpi::communicator& comm,
                  int master_rank = 0)
  : _db_name(db_name), _comm(comm), _cell_size(cell_size), _storage_map(expected_num_cells),
    _master_rank(master_rank)
  {
    
  }
  
  // Inserts only data from the master process
  void insert(const point_type& position, const T& x)
  {
    if(_comm.rank() == _master_rank)
    {
      grid_index_type corresponding_grid_entry = translate_to_grid_index(position);

      auto it = _storage_map.find(corresponding_grid_entry); 
      if(it != _storage_map.end())
      {
        *(it->second.get_serializer()) << position << x;
      }
      else
      {
        std::string filename = get_file_name(corresponding_grid_entry);

        storage_map_entry entry;
        entry.set_filename(filename);
        entry.set_file_handle(std::shared_ptr<std::fstream>(new std::fstream(filename.c_str(), 
                                                                           std::ios::binary | std::ios::out | std::ios::trunc)));

        assert(entry.get_file_handle()->is_open());

        *(entry.get_serializer()) << position << x;

        _storage_map.insert(std::make_pair(corresponding_grid_entry, entry));
      }
    }
  }
  
  template<class Function>
  void for_each(Function& f) const
  {
    for(const auto& db_file : _storage_map)
    {
      query_result_list_type cell_content;
      load_cell(db_file.first, cell_content);
      
      for(const auto& data_element : cell_content)
      {
        f(data_element.first, data_element.second);
      }
    }
  }
  
  void commit()
  {
    if(_comm.rank() == _master_rank)
    {
      for(const auto& element : _storage_map)
      {
        if(element.second.get_file_handle())
          element.second.get_file_handle()->flush();
      }
    }
    boost::mpi::broadcast(_comm, _storage_map, _master_rank);
  }
  
  void query_within_radius(const point_type& origin, 
                           util::scalar radius,
                           query_result_list_type& out) const
  {
    out.clear();
    
    // Load all necessary cells
    point_type cell_group_min_corner = origin;
    for(std::size_t i = 0; i < Dimension; ++i)
      cell_group_min_corner[i] -= radius;
      
    point_type cell_group_max_corner = origin;
    for(std::size_t i = 0; i < Dimension; ++i)
      cell_group_max_corner[i] += radius;
      
    grid_index_type min_corner_index = translate_to_grid_index(cell_group_min_corner);
    grid_index_type max_corner_index = translate_to_grid_index(cell_group_max_corner);
    
    query_result_list_type candidates;
      
    for(long long int i = min_corner_index[0]; i <= max_corner_index[0]; ++i)
    {
      for(long long int j = min_corner_index[1]; j <= max_corner_index[1]; ++j)
      {
        // Load cell at (i,j)
        grid_index_type current_cell = {i, j};
        load_cell(current_cell, candidates);
      }
    }
    
    for(const auto& element : candidates)
    {
      point_type pos = element.first;
      point_type delta = {pos[0] - origin[0], 
                          pos[1] - origin[1]};
      if(util::dot(delta, delta) <= util::square(radius))
        out.push_back(element);
    }
  }
  
  
private:
  
  std::string get_file_name(const grid_index_type& grid_pos)
  {
    std::string result = "spatial_grid."+_db_name;
    
    for(std::size_t i=0; i < Dimension; ++i)
    {
      result += "_";
      result += std::to_string(grid_pos[i]);
    }
    
    return result;
  }
  
  // Appends the content of a given cell to a result list.
  void load_cell(const grid_index_type& grid_pos, query_result_list_type& out) const
  {
    std::string file = _storage_map[grid_pos].get_filename();
    if(file.length() != 0)
    {
    
      std::ifstream cell_file(file.c_str(), std::ios::binary);

      if(cell_file.is_open())
      {
        boost::archive::binary_iarchive serializer(cell_file);
        
        point_type pos;
        T x;
        try
        {
          while(cell_file)
          {
            serializer >> pos;
            serializer >> x;

            out.push_back(std::make_pair(pos, x));
          }
        }
        catch(...)
        { /* The serializer will throw an error when it has reached the end of the file*/ }
      }
    }
    
  }
  
  grid_index_type translate_to_grid_index(const point_type& position) const
  {
    grid_index_type result;
    
    for(std::size_t i = 0; i < Dimension; ++i)
    {
      result[i] = std::round(position[i] / _cell_size[i]);
    }
    
    return result;
  }
  
  std::string _db_name;
  
  boost::mpi::communicator _comm;
  point_type _cell_size;
  
  mutable storage_map_type _storage_map;
  
  int _master_rank;
};

}

#endif	/* SPATIAL_GRID_DB_HPP */

