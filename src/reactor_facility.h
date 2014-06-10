#ifndef CYCLUS_REACTORS_REACTOR_FACILITY_H_
#define CYCLUS_REACTORS_REACTOR_FACILITY_H_

#include <string>

#include "cyclus.h"

namespace reactor {

/// @class ReactorFacility
///
/// This Facility is intended
/// as a skeleton to guide the implementation of new Facility
/// agents.
/// The ReactorFacility class inherits from the Facility class and is
/// dynamically loaded by the Agent class when requested.
///
/// @section intro Introduction
/// Place an introduction to the agent here.
///
/// @section agentparams Agent Parameters
/// Place a description of the required input parameters which define the
/// agent implementation.
///
/// @section optionalparams Optional Parameters
/// Place a description of the optional input parameters to define the
/// agent implementation.
///
/// @section detailed Detailed Behavior
/// Place a description of the detailed behavior of the agent. Consider
/// describing the behavior at the tick and tock as well as the behavior
/// upon sending and receiving materials and messages.
class ReactorFacility : public cyclus::Facility  {
 public:  
  /// Constructor for ReactorFacility Class
  /// @param ctx the cyclus context for access to simulation-wide parameters
  explicit ReactorFacility(cyclus::Context* ctx);

  /// The Prime Directive
  /// Generates code that handles all input file reading and restart operations
  /// (e.g., reading from the database, instantiating a new object, etc.).
  /// @warning The Prime Directive must have a space before it! (A fix will be
  /// in 2.0 ^TM)
  
  #pragma cyclus

  #pragma cyclus note {"doc": "A reactor facility is provided as a skeleton " \
                              "for the design of new facility agents."}

  /// A verbose printer for the ReactorFacility
  virtual std::string str();
  
  /// The handleTick function specific to the ReactorFacility.
  /// @param time the time of the tick  
  virtual void Tick();

  /// The handleTick function specific to the ReactorFacility.
  /// @param time the time of the tock
  virtual void Tock();

  std::set<cyclus::RequestPortfolio<cyclus::Material>::Ptr> GetMatlRequests();

  /// @brief Place accepted trade Materials in inventory
  virtual void AcceptMatlTrades(const std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                                cyclus::Material::Ptr> >& responses);

  /// This facility has one output commodity and one input commodity
  #pragma cyclus var {"tooltip": "input commodity", \
                      "doc": "commodity that the brightlite reactor consumes"}
  std::string in_commod;

  #pragma cyclus var {"tooltip": "output commodity", \
                      "doc": "commodity that the brightlite supplies"}
  std::string out_commod;

  #pragma cyclus var {"tooltip": "reactor libraries to load", \
                      "doc": "the reactor's burnup & criticality behavior to use"}
  std::vector<std::string> libraries; 
  
  #pragma cyclus var {"tooltip": "number of batches", \
                      "default": 3}
  int batches;

  #pragma cyclus var {"tooltip": "Non-leakage probability", \
                      "doc": "Varies from 0 - 1.", \
                      "default": 0.98}
  double nonleakage;

  #pragma cyclus var {"default": 1e299, \
                      "tooltip": "reactor maximum inventory size", \
                      "doc": "total maximum inventory size of the reactor"}
  double max_inv_size;

  #pragma cyclus var {"tooltip": "Target burnup", \
                      "units": "MWd/kgIHM"}
  double target_burnup;

  #pragma cyclus var {"default": 1000.0, \
                      "units": "MWe", \
                      "tooltip": "Electrical production."}
  double generated_power;

  #pragma cyclus var {"units": "kgIHM", \
                      "tooltip": "Total mass of the core."}
  double core_mass;

  #pragma cyclus var {"default": 0.33, \
                      "tooltip": "Thermal to electric conversion rate."}
  double efficiency;

  #pragma cyclus var {"capacity": "max_inv_size"}
  cyclus::toolkit::ResourceBuff inventory;

 private:
  int cycle_end_;
  fuelBundle fuel_library_;
};

}  // namespace reactor

#endif  // CYCLUS_REACTORS_REACTOR_FACILITY_H_