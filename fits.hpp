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

#ifndef FITS_HPP
#define	FITS_HPP

#include <fitsio.h>
#include <string>
#include "util.hpp"

namespace nanolens{

namespace util{

template<typename Scalar_type>
struct fits_datatype
{};

template<>
struct fits_datatype<float>
{
  static int image_type()
  { return FLOAT_IMG; }
  
  static int datatype()
  { return TFLOAT; }
};

template<>
struct fits_datatype<double>
{
  static int image_type()
  { return DOUBLE_IMG; }
  
  static int datatype()
  { return TDOUBLE; }
};


template<typename T>
class fits
{
public:
  fits(const std::string& filename)
  : _filename(filename) {}
  
  void save(const util::multi_array<T>& data) const
  {
    fitsfile* file;
    int status = 0;

    std::vector<long> naxes;
    
    for(std::size_t dim = 0; dim < data.get_dimension(); ++dim)
      naxes.push_back(data.get_extent_of_dimension(dim));

    // cfitsio will only overwrite files when their names are preceded by an
    // exclamation mark...
    std::string fitsio_filename = "!"+_filename;

    if (!fits_create_file(&file, fitsio_filename.c_str(), &status))
    {
      std::vector<long> fpixel(data.get_dimension(), 1);
      
      if (!fits_create_img(file, fits_datatype<T>::image_type(), 
                           naxes.size(), naxes.data(), &status))
      {
        fits_write_pix(file, fits_datatype<T>::datatype(), fpixel.data(), 
                       data.get_num_elements(), const_cast<T*>(data.data()), &status);
        
        fits_close_file(file, &status);
      }
      
    }
  }
  
  void load(const std::string& filename, std::size_t dimensions, util::multi_array<T>& out) const
  {
    fitsfile* file;
    int status = 0;
    int bitpix, naxis_flag;
    std::vector<long> naxes(dimensions, 0);
    
    if(!fits_open_file(&file, filename.c_str(), READONLY, &status))
    {
      if(!fits_get_img_param(file, dimensions, &bitpix, &naxis_flag, naxes.data(), 
                             &status))
      {
        out = util::multi_array<T>(naxes);
        
        // TODO
      }
    }
    
  }
  
private:
  std::string _filename;
};

}

}

#endif	/* FITS_HPP */

