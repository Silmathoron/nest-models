/*
 *  nngt_module.cpp
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

// include necessary NEST headers
#include "config.h"
#include "network.h"
#include "model.h"
#include "dynamicloader.h"
#include "genericmodel.h"
#include "booldatum.h"
#include "integerdatum.h"
#include "tokenarray.h"
#include "exceptions.h"
#include "sliexceptions.h"
#include "nestmodule.h"
#include "connector_model_impl.h"
#include "target_identifier.h"

// include headers with your own stuff
#include "nngt_module.h"
#include "aeif/aeif_cond_alpha_mod.h"
#include "aeif/aeif_psc_alpha.h"
#include "aeif/aeif_psc_exp.h"
#include "gp_aeif/gp_aeif_cond_alpha.h"
#include "gp_aeif/gp_aeif_cond_exp.h"
#include "gp_aeif/gp_aeif_psc_alpha.h"
#include "gp_aeif/gp_aeif_psc_exp.h"
#include "ps_aeif/ps_aeif_cond_alpha.h"
#include "ps_aeif/ps_aeif_cond_exp.h"
#include "ps_aeif/ps_aeif_psc_alpha.h"
#include "ps_aeif/ps_aeif_psc_exp.h"
#include "ps_iaf/ps_iaf_cond_alpha.h"
#include "ps_iaf/ps_iaf_psc_alpha.h"

// -- Interface to dynamic module loader ---------------------------------------

/*
 * The dynamic module loader must be able to find your module.
 * You make the module known to the loader by defining an instance of your
 * module class in global scope. This instance must have the name
 *
 * <modulename>_LTX_mod
 *
 * The dynamicloader can then load modulename and search for symbol "mod" in it.
 */

mynest::NngtModule nngt_module_LTX_mod;

// -- DynModule functions ------------------------------------------------------

mynest::NngtModule::NngtModule()
{
#ifdef LINKED_MODULE
  // register this module at the dynamic loader
  // this is needed to allow for linking in this module at compile time
  // all registered modules will be initialized by the main app's dynamic loader
  nest::DynamicLoaderModule::registerLinkedModule( this );
#endif
}

mynest::NngtModule::~NngtModule()
{
}

const std::string
mynest::NngtModule::name( void ) const
{
  return std::string( "NNGT module" ); // Return name of the module
}

const std::string
mynest::NngtModule::commandstring( void ) const
{
  // Instruct the interpreter to load nngt_module-init.sli
  return std::string( "(nngt_module-init) run" );
}

//-------------------------------------------------------------------------------------

void
mynest::NngtModule::init( SLIInterpreter* i )
{
  /* Register a neuron or device model.
   Give node type as template argument and the name as second argument.
   The first argument is always a reference to the network.
  */
  nest::register_model< aeif_cond_alpha_mod >( nest::NestModule::get_network(), "aeif_cond_alpha_mod" );
  nest::register_model< aeif_psc_alpha >( nest::NestModule::get_network(), "aeif_psc_alpha" );
  nest::register_model< aeif_psc_exp >( nest::NestModule::get_network(), "aeif_psc_exp" );
  nest::register_model< gp_aeif_cond_alpha >( nest::NestModule::get_network(), "gp_aeif_cond_alpha" );
  nest::register_model< gp_aeif_cond_exp >( nest::NestModule::get_network(), "gp_aeif_cond_exp" );
  nest::register_model< gp_aeif_psc_alpha >( nest::NestModule::get_network(), "gp_aeif_psc_alpha" );
  nest::register_model< gp_aeif_psc_exp >( nest::NestModule::get_network(), "gp_aeif_psc_exp" );
  nest::register_model< ps_aeif_cond_alpha >( nest::NestModule::get_network(), "ps_aeif_cond_alpha" );
  nest::register_model< ps_aeif_cond_exp >( nest::NestModule::get_network(), "ps_aeif_cond_exp" );
  nest::register_model< ps_aeif_psc_alpha >( nest::NestModule::get_network(), "ps_aeif_psc_alpha" );
  nest::register_model< ps_aeif_psc_exp >( nest::NestModule::get_network(), "ps_aeif_psc_exp" );
  nest::register_model< ps_iaf_cond_alpha >( nest::NestModule::get_network(), "ps_iaf_cond_alpha" );
  nest::register_model< ps_iaf_psc_alpha >( nest::NestModule::get_network(), "ps_iaf_psc_alpha" );

} // NngtModule::init()
