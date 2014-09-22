#ifndef CYCLUS_REACTORS_REACTOR_FACILITY_H_
#define CYCLUS_REACTORS_REACTOR_FACILITY_H_

#include <string>
#include <sstream>

#include "cyclus.h"
#include "burnupcalc.h"

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

  virtual std::set<cyclus::BidPortfolio<cyclus::Material>::Ptr> GetMatlBids(cyclus::CommodMap<cyclus::Material>::type& commod_requests);

  void GetMatlTrades(
    const std::vector< cyclus::Trade<cyclus::Material> >& trades,
    std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                          cyclus::Material::Ptr> >& responses);
  /// @brief Place accepted trade Materials in inventory
  virtual void AcceptMatlTrades(const std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                                cyclus::Material::Ptr> >& responses);

  double burnup_test(cyclus::Material::Ptr new_batch);

  void batch_reorder();
  /// This facility has one output commodity and one input commodity
  #pragma cyclus var {"tooltip": "input commodity", \
                      "doc": "commodity that the brightlite reactor consumes", \
                      "schematype": "token"}
  std::vector<std::string> in_commods;

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
  
  #pragma cyclus var {"default": 89197, \
                      "units": "cm2", \
                      "tooltip": "Total area of the fuel in core. Used for cylindrical flux calculation."}
  double fuel_area;
  
  #pragma cyclus var {"default": 5, \
                      "tooltip": "Delta to be used for cylindrical flux calculation."}
  double cylindrical_delta;
  
  #pragma cyclus var {"default": 0.0222, \
                      "units": "cm-1", \
                      "tooltip": "Macroscopic absorption cross section of the moderator."}
  double mod_Sig_a;
  
  #pragma cyclus var {"default": 3.46, \
                      "units": "cm-1", \
                      "tooltip": "Macroscopic transport cross section of the moderator."}
  double mod_Sig_tr;
  
  #pragma cyclus var {"default": 0.0, \
                      "units": "cm-1", \
                      "tooltip": "Macroscopic fission cross section of the moderator."}
  double mod_Sig_f;
  
  #pragma cyclus var {"default": 100, \
                      "units": "cm", \
                      "tooltip": "Radial thickness of the moderator used for cylindrical flux calculation."}
  double mod_thickness;
  
  #pragma cyclus var {"default": 3.94, \
                      "units": "cm-1", \
                      "tooltip": "Macroscopic transport cross section of the fuel."}
  double fuel_Sig_tr;
  
  #pragma cyclus var {"default": 10, \
                      "tooltip": "Timestep [days] for the burnup calculation."}
  double burnupcalc_timestep;
  
  #pragma cyclus var {"default": 1, \
                      "tooltip": "Flux calculation method. 1:Uniform, 2:Inv.Neut.Prod, 3:Cylindrical"}
  double flux_mode;

  #pragma cyclus var {"capacity": "max_inv_size"}
  cyclus::toolkit::ResourceBuff inventory;

  #pragma cyclus var{"default": 0.001, \
                     "tooltip": "The convergence requirement for the code"}
  double tolerence;

 private:
  int cycle_end_;
  fuelBundle fuel_library_;
  fuelBundle core_;
};

}  // namespace reactor

#endif  // CYCLUS_REACTORS_REACTOR_FACILITY_H_
