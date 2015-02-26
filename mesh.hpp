/* 
 * File:   mesh.hpp
 * Author: aksel
 *
 * Created on 15. Dezember 2014, 13:35
 */

#ifndef MESH_HPP
#define	MESH_HPP

#include <memory>
#include <vector>
#include <random>
#include <initializer_list>
#include "util.hpp"

namespace nanolens
{
  template<class Payload>
  class payload_bearing_type
  {
    Payload _p;
    
  public:
    payload_bearing_type()
    {}
    
    const Payload& get_payload() const {return _p; }
    Payload& get_payload() {return _p; }
    
    void set_payload(const Payload& p) { _p = p; }
  };
  
  template<class VertexPayload,
           class CellPayload>
  class mesh
  {
  public:
    class vertex;
    class edge;
    class cell;

    typedef std::shared_ptr<vertex> vertex_ptr;
    typedef std::shared_ptr<edge> edge_ptr;
    typedef std::shared_ptr<cell> cell_ptr;

    typedef typename std::vector<edge_ptr>::iterator edge_iterator;
    typedef typename std::vector<edge_ptr>::const_iterator const_edge_iterator;
    
    typedef typename std::vector<cell_ptr>::iterator cell_iterator;
    typedef typename std::vector<cell_ptr>::const_iterator const_cell_iterator;
    
    typedef typename std::vector<vertex_ptr>::iterator vertex_iterator;
    typedef typename std::vector<vertex_ptr>::const_iterator const_vertex_iterator;
    
    typedef unsigned long long cell_id_type;
    
    struct cell_address
    {
      cell_address(cell_id_type cell_id, std::size_t cell_refinement_level)
      : id(cell_id), refinement_level(cell_refinement_level) {}
      
      cell_id_type id;
      std::size_t refinement_level;
    };
    
    typedef typename std::vector<cell_address>::iterator cell_address_iterator;
    typedef typename std::vector<cell_address>::const_iterator const_cell_address_iterator;
    
    class mesh_database
    {
    public:
      mesh_database()
      {
        _cells.reserve(30);
      }
      
      const cell_ptr& get_cell(const cell_address& addr) const
      {
        return get_cell(addr.id, addr.refinement_level);
      }
      
      const cell_ptr& get_cell(cell_id_type id, std::size_t refinement_level) const
      {
        assert(refinement_level < _cells.size());
        assert(id < _cells[refinement_level].size());
        return _cells[refinement_level][id];
      }
      
      const cell_ptr& create_new_cell(const edge_ptr& e0,
                                      const edge_ptr& e1,
                                      const edge_ptr& e2,
                                      std::size_t refinement_level)
      {
        cell_id_type id = allocate_new_id(refinement_level);
        cell_ptr new_cell(new cell(e0, e1, e2, this, id, refinement_level));
        
        _cells[refinement_level].push_back(new_cell);
        return _cells[refinement_level].back();
      }
      
      const vertex_ptr& create_new_vertex(const util::vector2& position)
      {
        vertex_ptr new_vertex(new vertex(position));
        _vertices.push_back(new_vertex);
        return _vertices.back();
      }
      
      
      inline bool contains_refinement_level(std::size_t refinement_level) const
      {
        return refinement_level < _cells.size();
      }
      
      inline std::size_t get_num_contained_levels() const
      {
        return _cells.size();
      }
      
      cell_iterator begin_cells(std::size_t refinement_level)
      { 
        assert(refinement_level < _cells.size());
        return _cells[refinement_level].begin();
      }
      
      const_cell_iterator begin_cells(std::size_t refinement_level) const
      { 
        assert(refinement_level < _cells.size());
        return _cells[refinement_level].begin();
      }
      
      cell_iterator end_cells(std::size_t refinement_level)
      { 
        assert(refinement_level < _cells.size());
        return _cells[refinement_level].end();
      }
      
      const_cell_iterator end_cells(std::size_t refinement_level) const
      { 
        assert(refinement_level < _cells.size());
        return _cells[refinement_level].end();
      }
      
      vertex_iterator begin_vertices()
      { 
        return _vertices.begin();
      }
      
      const_vertex_iterator begin_vertices() const
      { 
        return _vertices.begin();
      }
      
      vertex_iterator end_vertices()
      { 
        return _vertices.end();
      }
      
      const_vertex_iterator end_vertices() const
      { 
        return _vertices.end();
      }
      
      std::size_t get_num_vertices() const
      {
        return _vertices.size();
      }
      
      template<class Function>
      void for_each_cell(Function f)
      {
        for(std::size_t ref_level = 0; ref_level < _cells.size(); ++ref_level)
          for_each_cell(ref_level, f);
      }
        
      template<class Function>
      void for_each_cell(std::size_t refinement_level, Function f)
      {
        assert(refinement_level < _cells.size());
        
        
        for(std::size_t i = 0; i < _cells[refinement_level].size(); ++i)
          if(_cells[refinement_level][i])
            if(_cells[refinement_level][i]->is_valid())
              f(_cells[refinement_level][i]);
      }
     
      
      
      //total number of cells, including invalidated cells
      std::size_t get_num_cells() const
      {
        std::size_t size = 0;
        for(auto&& lvl : _cells)
          size += lvl.size();
        
        return size;
      }
      
    private:
      
      cell_id_type allocate_new_id(std::size_t refinement_level)
      { 
        if(refinement_level >= _cells.size())
        {
          _cells.resize(refinement_level + 1);
          /*std::size_t num_new_levels_required = refinement_level - _cells.size() + 1;
          for(std::size_t i=0; i < num_new_levels_required;++i)
          {
            _cells.push_back(std::vector<cell_ptr>());
            // reserve the number of entries of the previous level
            _cells.back().reserve(_cells[_cells.size()-2].size());
          }*/
        }
        
        return _cells[refinement_level].size();
      }
      
      std::vector<std::vector<cell_ptr>> _cells;
      std::vector<vertex_ptr> _vertices;
    };
    
    class vertex : public payload_bearing_type<VertexPayload>
    {
    public:
      vertex()
      : _position({0.0,0.0}){}
      
      explicit vertex(const util::vector2& position)
      : _position(position) {}
      
      
      
      void set_position(const util::vector2& pos)
      { _position = pos; }
      
      
      const util::vector2& get_position() const
      {
        return _position;
      }
      
    private:
      util::vector2 _position;
    };
    
    class edge
    {
    public:
      edge(const vertex_ptr& v0, const vertex_ptr& v1)
      : _vertices({v0, v1})
      {
        _cells.reserve(2);
      }
      
      edge(const vertex_ptr& v0, const vertex_ptr& v1,
           const cell_address& c0, const cell_address& c1)
      : _vertices({v0, v1}), _cells({c0, c1})
      {}
      
      cell_address_iterator begin_cells()
      { return _cells.begin(); }
      
      const_cell_address_iterator begin_cells() const
      { return _cells.begin(); }
      
      cell_address_iterator end_cells()
      { return _cells.end(); }
      
      const_cell_address_iterator end_cells() const
      { return _cells.end(); }
      
      vertex_iterator begin_vertices()
      { return _vertices.begin(); }
      
      const_vertex_iterator begin_vertices() const
      { return _vertices.begin(); }
      
      vertex_iterator end_vertices()
      { return _vertices.end(); }
      
      const_vertex_iterator end_vertices() const
      { return _vertices.end(); }
      
      vertex_ptr get_vertex0() const
      {
        assert(_vertices.size() == 2);
        return _vertices[0];
      }
      
      vertex_ptr get_vertex1() const
      {
        assert(_vertices.size() == 2);
        return _vertices[1];
      }
      
      vertex_ptr get_other_vertex(const vertex_ptr& v) const
      {
        assert(_vertices.size() == 2);
        if(_vertices[0] == v)
          return _vertices[1];
        return _vertices[0];
      }
      
      bool ends_at_vertex(const vertex_ptr& v) const
      {
        assert(_vertices.size() == 2);
        if(_vertices[0] == v || _vertices[1] == v)
          return true;
        return false;
      }
      
      void add_cell(const cell_address& c)
      {
        _cells.push_back(c);
      }
      
      void set_cells(std::initializer_list<cell_address> cells)
      {
        _cells = std::vector<cell_address>(cells);
      }
      
    private:
      std::vector<vertex_ptr> _vertices;
      std::vector<cell_address> _cells;
    };
    
    class cell : public payload_bearing_type<CellPayload>
    {
    public:
      explicit cell(const edge_ptr& e0, const edge_ptr& e1, const edge_ptr& e2,
                    mesh_database* meshdb,
                    cell_id_type id,
                    std::size_t refinement_lvl = 0)
      : _edges({e0, e1, e2}), _refinement_level(refinement_lvl), _valid(true),
        _id(id), _mesh(meshdb)
      {
        assert(e0);
        assert(e1);
        assert(e2);
        
        update_vertices();
        e0->add_cell(get_address());
        e1->add_cell(get_address());
        e2->add_cell(get_address());
      }
      
      edge_iterator begin_edges()
      { return _edges.begin(); }
      
      const_edge_iterator begin_edges() const
      { return _edges.begin(); }
      
      edge_iterator end_edges()
      { return _edges.end(); }
      
      const_edge_iterator end_edges() const
      { return _edges.end(); }
      
      vertex_iterator begin_vertices()
      { return _vertices.begin(); }
      
      const_vertex_iterator begin_vertices() const
      { return _vertices.begin(); }
      
      vertex_iterator end_vertices()
      { return _vertices.end(); }
      
      const_vertex_iterator end_vertices() const
      { return _vertices.end(); }
      
      
      const vertex_ptr& get_vertex(std::size_t local_vertex_id) const
      {
        assert(local_vertex_id < _vertices.size());
        return _vertices[local_vertex_id];
      }
      
      const edge_ptr& get_edge(std::size_t local_edge_id) const
      {
        assert(local_edge_id < _edges.size());
        return _edges[local_edge_id];
      }
      
      std::size_t get_refinement_level() const
      { return _refinement_level; }
      
      cell_id_type get_id() const
      { return _id; }
      
      cell_address get_address() const
      { return cell_address(_id, _refinement_level); }
      
      void invalidate()
      { _valid = false; }
      
      bool is_valid() const
      { return _valid; }
      
      typedef std::vector<cell_ptr> refined_cells_container;
      typedef std::vector<vertex_ptr> new_vertices_container;
      
      // Refine this cell. This cell is obsolete afterwards and will be invalid!
      // Note that the refinement level for this cell is increased by 2,
      // and the refinement level for adjacent bisected cells is increased by 1.
      void refine(refined_cells_container& new_cells, new_vertices_container& new_vertices)
      {
        
        new_vertices.clear();
        new_cells.clear();
        
        std::vector<edge_ptr> half_edges;
        half_edges.reserve(6);

        
        for(auto& current_edge : _edges)
        {
          util::vector2 center_of_edge;
          util::average(current_edge->get_vertex0()->get_position(),
                        current_edge->get_vertex1()->get_position(), 
                        center_of_edge);
          
          vertex_ptr subdivided_vertex = _mesh->create_new_vertex(center_of_edge);
          new_vertices.push_back(subdivided_vertex);
          
          // First, bisect the adjacent cells          
          cell_ptr adjacent_cell = get_adjacent_cell(current_edge);
          
          if(adjacent_cell != nullptr)
          {
            assert(adjacent_cell.get() != this);

            edge_ptr half_edge0, half_edge1;

            adjacent_cell->bisect_cell(current_edge, subdivided_vertex, new_cells,
                                      half_edge0, half_edge1);

            half_edges.push_back(half_edge0);
            half_edges.push_back(half_edge1);              

          }
          else
          {
            edge_ptr half_edge0(new edge(current_edge->get_vertex0(), subdivided_vertex));
            edge_ptr half_edge1(new edge(current_edge->get_vertex1(), subdivided_vertex));
 
            half_edges.push_back(half_edge0);
            half_edges.push_back(half_edge1);
          }
        }
        
        
        std::vector<edge_ptr> new_edges;
        
        this->invalidate();
        
        // create new cells
        for(auto&& current_vertex : _vertices)
        {
          // Find half edges containing this vertex
          std::vector<edge_ptr> adjacent_edges;
          adjacent_edges.reserve(2);
          
          for(auto&& current_edge : half_edges)
            if(current_edge->ends_at_vertex(current_vertex))
              adjacent_edges.push_back(current_edge);
          
          assert(adjacent_edges.size() == 2);
          
          vertex_ptr new_endpoint0 = adjacent_edges[0]->get_other_vertex(current_vertex);
          vertex_ptr new_endpoint1 = adjacent_edges[1]->get_other_vertex(current_vertex);
          edge_ptr new_edge(new edge(new_endpoint0, new_endpoint1));
          
          cell_ptr new_cell = _mesh->create_new_cell(adjacent_edges[0], adjacent_edges[1],
                                                    new_edge, _refinement_level+2);
          
          new_cells.push_back(new_cell);
          new_edges.push_back(new_edge);           
        }
        
        // create central new cell
        assert(new_edges.size() == 3);
        cell_ptr central_cell = _mesh->create_new_cell(new_edges[0], new_edges[1], new_edges[2],
                                       _refinement_level+2);
        
        
        new_cells.push_back(central_cell);
        
      }
      
      

    private:
      void generate_vertex_list(std::vector<vertex_ptr>& out) const
      {
        out.clear();
        out.reserve(3);
        for(auto&& current_edge : _edges)
        {
          vertex_ptr vertex0 = current_edge->get_vertex0();
          vertex_ptr vertex1 = current_edge->get_vertex1();
          
          if(std::find(out.begin(), out.end(), vertex0) == out.end())
            out.push_back(vertex0);
          
          if(std::find(out.begin(), out.end(), vertex1) == out.end())
            out.push_back(vertex1);
        }
      }
      
      void bisect_cell(const edge_ptr& bisected_edge,
                       const vertex_ptr& new_vertex,
                       std::vector<cell_ptr>& new_cells,
                       edge_ptr& half_edge0,
                       edge_ptr& half_edge1)
      {
        // Find opposing vertex
        vertex_ptr opposing_vertex = nullptr;

        for(auto&& current_vertex : _vertices)
        {
          if(current_vertex != bisected_edge->get_vertex0() &&
             current_vertex != bisected_edge->get_vertex1())
          {
            // we have found the opposing vertex!
            opposing_vertex = current_vertex;
          }
        }

        if(opposing_vertex == nullptr)
          throw std::runtime_error("refine(): could not find opposing vertex");

        // Create new edge
        edge_ptr bisecting_edge(new edge(new_vertex, opposing_vertex));
        half_edge0 = edge_ptr(new edge(new_vertex, bisected_edge->get_vertex0()));
        half_edge1 = edge_ptr(new edge(new_vertex, bisected_edge->get_vertex1()));
        
        edge_ptr cell_border0 = find_edge_with_vertices(opposing_vertex, bisected_edge->get_vertex0());
        edge_ptr cell_border1 = find_edge_with_vertices(opposing_vertex, bisected_edge->get_vertex1());
        
        cell_ptr new_cell0 = _mesh->create_new_cell(half_edge0, cell_border0, bisecting_edge, _refinement_level + 1);
        cell_ptr new_cell1 = _mesh->create_new_cell(half_edge1, cell_border1, bisecting_edge, _refinement_level + 1);       
        
        cell_ptr adjacent_cell_border0 = get_adjacent_cell(cell_border0);
        if(adjacent_cell_border0 != nullptr)
          cell_border0->set_cells({new_cell0->get_address(), adjacent_cell_border0->get_address()});
        else
          cell_border0->set_cells({new_cell0->get_address()});
          
        cell_ptr adjacent_cell_border1 = get_adjacent_cell(cell_border1);
        if(adjacent_cell_border1 != nullptr)
          cell_border1->set_cells({new_cell1->get_address(), adjacent_cell_border1->get_address()});
        else
          cell_border1->set_cells({new_cell1->get_address()});
          
        this->invalidate();
        
        new_cells.push_back(new_cell0);
        new_cells.push_back(new_cell1);
        
      }
      
      cell_ptr get_adjacent_cell(const edge_ptr& edge)
      {
        for(auto current_cell_address = edge->begin_cells();
                current_cell_address != edge->end_cells();
                ++current_cell_address)
        {
          if(current_cell_address->id != this->get_id() || 
             current_cell_address->refinement_level != this->get_refinement_level())
            return _mesh->get_cell(*current_cell_address);
        }
        return cell_ptr();
      }
      
      
      edge_ptr find_edge_with_vertices(const vertex_ptr& v0, const vertex_ptr& v1) const
      {
        for(auto&& current_edge : _edges)
        {
          if((current_edge->get_vertex0() == v0 && current_edge->get_vertex1() == v1)
             || (current_edge->get_vertex0() == v1 && current_edge->get_vertex1() == v0))
            return current_edge;
        }
        return nullptr;
      }
      
      void update_vertices()
      {
        generate_vertex_list(_vertices);
      }
      
      std::vector<edge_ptr> _edges;
      std::vector<vertex_ptr> _vertices;
      
      std::size_t _refinement_level;
      cell_id_type _id;
      bool _valid;
      
      mesh_database* _mesh;
    };
    
    mesh_database& get_meshdb()
    {
      return _meshdb;
    }
    
    void create_square_mesh(const util::vector2& pos, 
                            util::scalar radius,
                            std::size_t num_vertices_per_dim,
                            util::scalar wiggles_strength)
    {
      assert(num_vertices_per_dim >= 2);
      
      util::vector2 shift;

      util::multi_array<vertex_ptr> vertices(num_vertices_per_dim, num_vertices_per_dim);
      util::scalar stepwidth = 2 * radius / (num_vertices_per_dim - 1);
      
      std::default_random_engine rng_engine(_rd());
      std::uniform_real_distribution<util::scalar> random_variation(-wiggles_strength * stepwidth,
                                                                    wiggles_strength * stepwidth);
      
      util::multi_array<edge_ptr> horizontal_edges(num_vertices_per_dim - 1, num_vertices_per_dim);
      util::multi_array<edge_ptr> vertical_edges(num_vertices_per_dim, num_vertices_per_dim - 1);
      util::multi_array<edge_ptr> diagonal_edges(num_vertices_per_dim - 1, num_vertices_per_dim - 1);
      
      // Create new vertices
      for(std::size_t x = 0; x < num_vertices_per_dim; ++x)
      {
        for(std::size_t y = 0; y < num_vertices_per_dim; ++y)
        {
          shift[0] = -radius + stepwidth * x + random_variation(rng_engine);
          shift[1] = -radius + stepwidth * y + random_variation(rng_engine);
          
          util::vector2 vertex_position = pos;
          util::add(vertex_position, shift);

          std::size_t index[] = {x, y};

          vertex_ptr new_vertex = _meshdb.create_new_vertex(vertex_position);
          vertices[index] = new_vertex;

        }
      }
      
      // Create edges and cells from the new vertices
      for(std::size_t x = 0; x < num_vertices_per_dim - 1; ++x)
      {
        for(std::size_t y = 0; y < num_vertices_per_dim - 1; ++y)
        {
          std::size_t horizontal_edge_idx0 [] = {x, y};
          std::size_t horizontal_edge_idx1 [] = {x, y + 1};
          std::size_t vertical_edge_idx0 [] = {x, y};
          std::size_t vertical_edge_idx1 [] = {x + 1, y};
          std::size_t diagonal_edge_idx [] = {x, y};
          
          std::size_t vertex_idx0 [] = {x, y};
          std::size_t vertex_idx1 [] = {x + 1, y};
          std::size_t vertex_idx2 [] = {x + 1, y + 1};
          std::size_t vertex_idx3 [] = {x, y + 1};
          
          vertex_ptr vertex0 = vertices[vertex_idx0];
          vertex_ptr vertex1 = vertices[vertex_idx1];
          vertex_ptr vertex2 = vertices[vertex_idx2];
          vertex_ptr vertex3 = vertices[vertex_idx3];
          
          if(!horizontal_edges[horizontal_edge_idx0])
            horizontal_edges[horizontal_edge_idx0] = edge_ptr(new edge(vertex0, vertex1));
          
          if(!horizontal_edges[horizontal_edge_idx1])
            horizontal_edges[horizontal_edge_idx1] = edge_ptr(new edge(vertex2, vertex3));
          
          if(!vertical_edges[vertical_edge_idx0])
            vertical_edges[vertical_edge_idx0] = edge_ptr(new edge(vertex3, vertex0));
          
          if(!vertical_edges[vertical_edge_idx1])
            vertical_edges[vertical_edge_idx1] = edge_ptr(new edge(vertex1, vertex2));
          
          // It is impossible that the diagonal edge has been created already
          edge_ptr diagonal_edge;
          
          if(x % 2)
            // Alternate the direction of the diagonal
            diagonal_edge = edge_ptr(new edge(vertex0, vertex2));
          else
            diagonal_edge = edge_ptr(new edge(vertex1, vertex3));           

          diagonal_edges[diagonal_edge_idx] = diagonal_edge;
          
          edge_ptr horizontal_edge_top = horizontal_edges[horizontal_edge_idx0];
          edge_ptr horizontal_edge_bottom = horizontal_edges[horizontal_edge_idx1];
          edge_ptr vertical_edge_left = vertical_edges[vertical_edge_idx0];
          edge_ptr vertical_edge_right = vertical_edges[vertical_edge_idx1];
          
          // Create the cells. The cell constructor will take care of
          // making the edges point to the new cells
          if(x % 2)
          {
            _meshdb.create_new_cell(horizontal_edge_top,
                                   vertical_edge_right,
                                   diagonal_edge, 0);

            _meshdb.create_new_cell(diagonal_edge,
                                    horizontal_edge_bottom,
                                    vertical_edge_left, 0);
          }
          else
          {
            _meshdb.create_new_cell(horizontal_edge_top,
                                   vertical_edge_left,
                                   diagonal_edge, 0);

            _meshdb.create_new_cell(diagonal_edge,
                                    horizontal_edge_bottom,
                                    vertical_edge_right, 0);
          }
          
        }
      }

    }
    
    void create_circle_mesh(const util::vector2& pos, 
                            util::scalar radius,
                            std::size_t num_vertices_per_dim,
                            util::scalar wiggles_strength)
    {
      // Create square mesh, then shift vertices
      create_square_mesh(pos, radius, num_vertices_per_dim, wiggles_strength);
      
      for(auto current_vertex = _meshdb.begin_vertices();
              current_vertex != _meshdb.end_vertices();
              ++current_vertex)
      {
        util::vector2 vertex_position = (*current_vertex)->get_position();
        util::scalar vertex_radius = std::max(std::abs(pos[0] - vertex_position[0]),
                                              std::abs(pos[1] - vertex_position[1]));
        
        util::vector2 new_position = vertex_position;
        util::normalize(new_position);
        util::scale(new_position, vertex_radius);
        util::add(new_position, pos);
        
        (*current_vertex)->set_position(new_position);
      }
    }
  private:
    mesh_database _meshdb;
    
    std::random_device _rd;

  };
  
  template<class VertexPayload, class CellPayload>
  using mesh_ptr = std::shared_ptr<mesh<VertexPayload, CellPayload>>;
}

#endif	/* MESH_HPP */

