/*
 *  slice_ring_buffer.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "slice_ring_buffer.h"

#include <limits>
#include <cmath>

using namespace nest;

mynest::SliceRingBuffer::SliceRingBuffer()
  : refract_( std::numeric_limits< long_t >::max(), 0, 0 )
{
  //  resize();  // sets up queue_
}

void
mynest::SliceRingBuffer::resize()
{
  long_t newsize = static_cast< long_t >(
    std::ceil( static_cast< double >( Scheduler::get_min_delay() + Scheduler::get_max_delay() )
      / Scheduler::get_min_delay() ) );
  if ( queue_.size() != static_cast< ulong_t >( newsize ) )
  {
    queue_.resize( newsize );
    clear();
  }

#ifndef HAVE_STL_VECTOR_CAPACITY_BASE_UNITY
  // create 1-element buffers
  for ( size_t j = 0; j < queue_.size(); ++j )
    queue_[ j ].reserve( 1 );
#endif
}

void
mynest::SliceRingBuffer::clear()
{
  for ( size_t j = 0; j < queue_.size(); ++j )
    queue_[ j ].clear();
}

void
mynest::SliceRingBuffer::prepare_delivery()
{
  // vector to deliver from in this slice
  deliver_ = &( queue_[ Scheduler::get_slice_modulo( 0 ) ] );

  // sort events, first event last
  std::sort( deliver_->begin(), deliver_->end(), std::greater< SpikeInfo >() );
}

void
mynest::SliceRingBuffer::discard_events()
{
  // vector to deliver from in this slice
  deliver_ = &( queue_[ Scheduler::get_slice_modulo( 0 ) ] );

  deliver_->clear();
}
