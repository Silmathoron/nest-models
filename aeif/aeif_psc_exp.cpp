/*
 *  aeif_psc_exp.cpp
 *
 */

#include "aeif_psc_exp.h"
#include "nest_names.h" // in aeif but not in example

#ifdef HAVE_GSL_1_11

#include "universal_data_logger_impl.h"

#include "exceptions.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include <limits>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <cstdio>

//~ #include "lockptrdatum.h" // in example but not in aeif

using namespace nest;


/* ----------------------------------------------------------------
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap< mynest::aeif_psc_exp > mynest::aeif_psc_exp::recordablesMap_;

namespace nest // template specialization must be placed in namespace
{
  // Override the create() method with one call to RecordablesMap::insert_()
  // for each quantity to be recorded.
  template <>
  void
  RecordablesMap< mynest::aeif_psc_exp >::create()
  {
    insert_( names::V_m, &mynest::aeif_psc_exp::get_y_elem_< mynest::aeif_psc_exp::State_::V_M > );
    insert_( names::I_ex, &mynest::aeif_psc_exp::get_y_elem_< mynest::aeif_psc_exp::State_::I_EXC > );
    insert_( names::I_in, &mynest::aeif_psc_exp::get_y_elem_< mynest::aeif_psc_exp::State_::I_INH > );
    insert_( names::w, &mynest::aeif_psc_exp::get_y_elem_< mynest::aeif_psc_exp::State_::W > );
  }
}


/* ----------------------------------------------------------------
 * Dynamics for gsl_odeiv
 * ---------------------------------------------------------------- */

extern "C" int
mynest::aeif_psc_exp_dynamics( double, const double y[], double f[], void* pnode )
{
  // a shorthand
  typedef mynest::aeif_psc_exp::State_ S;

  // get access to node so we can almost work as in a member function
  assert( pnode );
  const mynest::aeif_psc_exp& node = *( reinterpret_cast< mynest::aeif_psc_exp* >( pnode ) );

  // y[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.y[].

  // The following code is verbose for the sake of clarity. We assume that a
  // good compiler will optimize the verbosity away ...

  // shorthand for state variables
  const double_t& V = y[ S::V_M ];
  const double_t& I_ex = y[ S::I_EXC ];
  const double_t& I_in = y[ S::I_INH ];
  const double_t& w = y[ S::W ];

  // We pre-compute the argument of the exponential
  const double_t exp_arg = ( V - node.P_.V_th ) / node.P_.Delta_T;

  // Upper bound for exponential argument to avoid numerical instabilities
  const double_t MAX_EXP_ARG = 10.;

  // If the argument is too large, we clip it.
  const double_t I_spike = node.P_.Delta_T * std::exp( std::min( exp_arg, MAX_EXP_ARG ) );

  // dv/dt
  f[ S::V_M ] = ( -node.P_.g_L * ( ( V - node.P_.E_L ) - I_spike ) + I_ex + I_in - w
            + node.P_.I_e + node.B_.I_stim_ ) / node.P_.C_m;

  f[ S::I_EXC ] = -I_ex / node.P_.tau_syn_ex; // Synaptic Conductance (nS)

  f[ S::I_INH ] = -I_in / node.P_.tau_syn_in;// Synaptic Conductance (nS)

  // Adaptation current w.
  f[ S::W ] = ( node.P_.a * ( V - node.P_.E_L ) - w ) / node.P_.tau_w;

  return GSL_SUCCESS;
}


/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

mynest::aeif_psc_exp::Parameters_::Parameters_()
  : V_peak_( 0.0 )   // mV, should not be larger that V_th+10
  , V_reset_( -60.0 ) // mV
  , t_ref_( 0.0 )    // ms
  , g_L( 30.0 )     // nS
  , C_m( 281.0 )    // pF
  , E_L( -70.6 )    // mV
  , Delta_T( 2.0 )   // mV
  , tau_w( 144.0 )   // ms
  , a( 4.0 )       // nS
  , b( 80.5 )      // pA
  , V_th( -50.4 )    // mV
  , tau_syn_ex( 0.2 ) // ms
  , tau_syn_in( 2.0 ) // ms
  , I_e( 0.0 )      // pA
  , gsl_error_tol( 1e-6 )
{
}

mynest::aeif_psc_exp::State_::State_( const Parameters_& p )
  : r_( 0 )
{
  y_[ 0 ] = p.E_L;
  for ( size_t i = 1; i < STATE_VEC_SIZE; ++i )
   y_[ i ] = 0;
}

mynest::aeif_psc_exp::State_::State_( const State_& s )
  : r_( s.r_ )
{
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
   y_[ i ] = s.y_[ i ];
}

mynest::aeif_psc_exp::State_& mynest::aeif_psc_exp::State_::operator=( const State_& s )
{
  assert( this != &s ); // would be bad logical error in program

  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
   y_[ i ] = s.y_[ i ];
  r_ = s.r_;
  return *this;
}


/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
mynest::aeif_psc_exp::Parameters_::get( DictionaryDatum& d ) const
{
  def< double >( d, names::C_m, C_m );
  def< double >( d, names::V_th, V_th );
  def< double >( d, names::t_ref, t_ref_ );
  def< double >( d, names::g_L, g_L );
  def< double >( d, names::E_L, E_L );
  def< double >( d, names::V_reset, V_reset_ );
  def< double >( d, names::tau_syn_ex, tau_syn_ex );
  def< double >( d, names::tau_syn_in, tau_syn_in );
  def< double >( d, names::a, a );
  def< double >( d, names::b, b );
  def< double >( d, names::Delta_T, Delta_T );
  def< double >( d, names::tau_w, tau_w );
  def< double >( d, names::I_e, I_e );
  def< double >( d, names::V_peak, V_peak_ );
  def< double >( d, names::gsl_error_tol, gsl_error_tol );
}

void
mynest::aeif_psc_exp::Parameters_::set( const DictionaryDatum& d )
{
  updateValue< double >( d, names::V_th, V_th );
  updateValue< double >( d, names::V_peak, V_peak_ );
  updateValue< double >( d, names::t_ref, t_ref_ );
  updateValue< double >( d, names::E_L, E_L );
  updateValue< double >( d, names::V_reset, V_reset_ );

  updateValue< double >( d, names::C_m, C_m );
  updateValue< double >( d, names::g_L, g_L );

  updateValue< double >( d, names::tau_syn_ex, tau_syn_ex );
  updateValue< double >( d, names::tau_syn_in, tau_syn_in );

  updateValue< double >( d, names::a, a );
  updateValue< double >( d, names::b, b );
  updateValue< double >( d, names::Delta_T, Delta_T );
  updateValue< double >( d, names::tau_w, tau_w );

  updateValue< double >( d, names::I_e, I_e );

  updateValue< double >( d, names::gsl_error_tol, gsl_error_tol );

  if ( V_peak_ <= V_th )
    throw BadProperty( "V_peak must be larger than threshold." );

  if ( V_reset_ >= V_peak_ )
    throw BadProperty( "Ensure that: V_reset < V_peak ." );

  if ( C_m <= 0 )
    throw BadProperty( "Capacitance must be strictly positive." );

  if ( t_ref_ < 0 )
    throw BadProperty( "Refractory time cannot be negative." );

  if ( tau_syn_ex <= 0 || tau_syn_in <= 0 || tau_w <= 0 )
    throw BadProperty( "All time constants must be strictly positive." );

  if ( gsl_error_tol <= 0. )
    throw BadProperty( "The gsl_error_tol must be strictly positive." );
}

void
mynest::aeif_psc_exp::State_::get( DictionaryDatum& d ) const
{
  def< double >( d, names::V_m, y_[ V_M ] );
  def< double >( d, names::I_ex, y_[ I_EXC ] );
  def< double >( d, names::I_in, y_[ I_INH ] );
  def< double >( d, names::w, y_[ W ] );
}

void
mynest::aeif_psc_exp::State_::set( const DictionaryDatum& d, const Parameters_& )
{
  updateValue< double >( d, names::V_m, y_[ V_M ] );
  updateValue< double >( d, names::I_ex, y_[ I_EXC ] );
  updateValue< double >( d, names::I_in, y_[ I_INH ] );
  updateValue< double >( d, names::w, y_[ W ] );
  
  if ( y_[ I_EXC ] < 0 || y_[ I_INH ] < 0 )
   throw BadProperty( "Currents are positive by definition." );
}

mynest::aeif_psc_exp::Buffers_::Buffers_( aeif_psc_exp& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

mynest::aeif_psc_exp::Buffers_::Buffers_( const Buffers_&, aeif_psc_exp& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}


/* ----------------------------------------------------------------
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

mynest::aeif_psc_exp::aeif_psc_exp()
  : Archiving_Node()
  , P_()
  , S_( P_ )
  , B_( *this )
{
  recordablesMap_.create();
}

mynest::aeif_psc_exp::aeif_psc_exp( const aeif_psc_exp& n )
  : Archiving_Node( n )
  , P_( n.P_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
mynest::aeif_psc_exp::init_state_( const Node& proto )
{
  const aeif_psc_exp& pr = downcast< aeif_psc_exp >( proto );
  S_ = pr.S_;
}

void
mynest::aeif_psc_exp::init_buffers_()
{
  B_.spike_exc_.clear(); // includes resize
  B_.spike_inh_.clear(); // includes resize
  B_.currents_.clear();  // includes resize
  Archiving_Node::clear_history();

  B_.logger_.reset();

  B_.step_ = Time::get_resolution().get_ms();

  // We must integrate this model with high-precision to obtain decent results
  B_.IntegrationStep_ = std::min( 0.01, B_.step_ );

  if ( B_.s_ == 0 )
    B_.s_ = gsl_odeiv_step_alloc( gsl_odeiv_step_rkf45, State_::STATE_VEC_SIZE );
  else
    gsl_odeiv_step_reset( B_.s_ );

  if ( B_.c_ == 0 )
    B_.c_ = gsl_odeiv_control_yp_new( P_.gsl_error_tol, P_.gsl_error_tol );
  else
    gsl_odeiv_control_init( B_.c_, P_.gsl_error_tol, P_.gsl_error_tol, 0.0, 1.0 );

  if ( B_.e_ == 0 )
    B_.e_ = gsl_odeiv_evolve_alloc( State_::STATE_VEC_SIZE );
  else
    gsl_odeiv_evolve_reset( B_.e_ );

  B_.sys_.function = aeif_psc_exp_dynamics;
  B_.sys_.jacobian = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params = reinterpret_cast< void* >( this );

  B_.I_stim_ = 0.0;
}

void
mynest::aeif_psc_exp::calibrate()
{
  B_.logger_.init(); // ensures initialization in case mm connected after Simulate

  V_.I0_ex_ = 1.0 * numerics::e / P_.tau_syn_ex;
  V_.I0_in_ = 1.0 * numerics::e / P_.tau_syn_in;
  V_.RefractoryCounts_ = Time( Time::ms( P_.t_ref_ ) ).get_steps();
  assert( V_.RefractoryCounts_ >= 0 ); // since t_ref_ >= 0, this can only fail in error
}

/* ----------------------------------------------------------------
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void
mynest::aeif_psc_exp::update( const Time& origin, const nest::long_t from, const nest::long_t to )
{
  assert( to >= 0 && ( delay ) from < Scheduler::get_min_delay() );
  assert( from < to );
  assert( State_::V_M == 0 );
  
  for ( long_t lag = from; lag < to; ++lag )
  {
    double t = 0.;

    if ( S_.r_ > 0 )
      --S_.r_;

    // numerical integration with adaptive step size control:
    // ------------------------------------------------------
    // gsl_odeiv_evolve_apply performs only a single numerical
    // integration step, starting from t and bounded by step;
    // the while-loop ensures integration over the whole simulation
    // step (0, step] if more than one integration step is needed due
    // to a small integration step size;
    // note that (t+IntegrationStep > step) leads to integration over
    // (t, step] and afterwards setting t to step, but it does not
    // enforce setting IntegrationStep to step-t

    while ( t < B_.step_ )
    {
      // propagate the ODE
      const int status = gsl_odeiv_evolve_apply( B_.e_,
      B_.c_,
      B_.s_,
      &B_.sys_,         // system of ODE
      &t,             // from t
      B_.step_,         // to t <= step
      &B_.IntegrationStep_, // integration step size
      S_.y_ );          // neuronal state

      if ( status != GSL_SUCCESS )
        throw GSLSolverFailure( get_name(), status );

      // check for unreasonable values; we allow V_M to explode
      if ( S_.y_[ State_::V_M ] < -1e3 || S_.y_[ State_::W ] < -1e6 || S_.y_[ State_::W ] > 1e6 )
        throw NumericalInstability( get_name() );

      // spikes are handled inside the while-loop
      // due to spike-driven adaptation
      if ( S_.r_ > 0 )
        S_.y_[ State_::V_M ] = P_.V_reset_; // only V_m is frozen
      else if ( S_.y_[ State_::V_M ] >= P_.V_peak_ )
      {
        S_.y_[ State_::V_M ] = P_.V_reset_;
        S_.y_[ State_::W ] += P_.b;
        S_.r_ = V_.RefractoryCounts_;
        
        set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );
        SpikeEvent se;
        network()->send( *this, se, lag );
      }
    }

    // influence of received spikes on post-synaptic conductances
    S_.y_[ State_::I_EXC ] += B_.spike_exc_.get_value( lag ) * V_.I0_ex_;
    S_.y_[ State_::I_INH ] -= B_.spike_inh_.get_value( lag ) * V_.I0_in_;

    // set new input current
    B_.I_stim_ = B_.currents_.get_value( lag );

    // log state data
    B_.logger_.record_data( origin.get_steps() + lag );
  }
}

void
mynest::aeif_psc_exp::handle( SpikeEvent& e )
{
  assert( e.get_delay() > 0 );

  if ( e.get_weight() > 0.0 )
    B_.spike_exc_.add_value( e.get_rel_delivery_steps( network()->get_slice_origin() ),
      e.get_weight() * e.get_multiplicity() );
  else
    B_.spike_inh_.add_value( e.get_rel_delivery_steps( network()->get_slice_origin() ),
      -e.get_weight() * e.get_multiplicity() ); // keep conductances positive
}

void
mynest::aeif_psc_exp::handle( CurrentEvent& e )
{
  assert( e.get_delay() > 0 );

  const double_t c = e.get_current();
  const double_t w = e.get_weight();

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value( e.get_rel_delivery_steps( network()->get_slice_origin() ), w * c );
}

// Do not move this function as inline to h-file. It depends on
// universal_data_logger_impl.h being included here.
void
mynest::aeif_psc_exp::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}

#endif // HAVE_GSL_1_11
