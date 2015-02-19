#ifndef CYCLUS_REPROCESSES_REPROCESS_FACILITY_H_
#define CYCLUS_REPROCESSES_REPROCESS_FACILITY_H_

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include "cyclus.h"
#include "burnupcalc.h"

namespace reprocess {

/// @class ReprocessFacility
///
/// This Facility is intended
/// as a skeleton to guide the implementation of new Facility
/// agents.
/// The ReprocessFacility class inherits from the Facility class and is
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
class ReprocessFacility : public cyclus::Facility  {
 public:
  /// Constructor for ReprocessFacility Class
  /// @param ctx the cyclus context for access to simulation-wide parameters
  explicit ReprocessFacility(cyclus::Context* ctx);

  /// The Prime Directive
  /// Generates code that handles all input file reading and restart operations
  /// (e.g., reading from the database, instantiating a new object, etc.).
  /// @warning The Prime Directive must have a space before it! (A fix will be
  /// in 2.0 ^TM)

  #pragma cyclus

  #pragma cyclus note {"doc": "A reprocess facility is provided as a skeleton " \
                              "for the design of new facility agents.", \
                       "niche": "reprocessing"}

  /// A verbose printer for the ReprocessFacility
  virtual std::string str();

  /// The handleTick function specific to the ReprocessFacility.
  /// @param time the time of the tick
  virtual void Tick();

  /// The handleTick function specific to the ReprocessFacility.
  /// @param time the time of the tock
  virtual void Tock();

  virtual void Build(cyclus::Agent* parent) {
    Facility::Build(parent);
    if (lifetime() >= 0) {
      context()->SchedDecom(this, exit_time());
    }
  }

  std::set<cyclus::RequestPortfolio<cyclus::Material>::Ptr> GetMatlRequests();

  std::vector<std::map<int, double>> out_eff; //isotopic efficiencies of each output
  double pi;

  virtual std::set<cyclus::BidPortfolio<cyclus::Material>::Ptr> GetMatlBids(cyclus::CommodMap<cyclus::Material>::type& commod_requests);

  void GetMatlTrades(const std::vector< cyclus::Trade<cyclus::Material> >& trades,
                     std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                     cyclus::Material::Ptr> >& responses);

  /// @brief Place accepted trade Materials in inventory
  virtual void AcceptMatlTrades(const std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                                cyclus::Material::Ptr> >& responses);

  /// This facility has two output commodity and one input commodity
  #pragma cyclus var {"tooltip": "input commodity", \
                      "doc": "commodity that the brightlite reprocess facility consumes", \
                      "schematype": "token", \
                      "uitype": ["oneOrMore", "incommodity"]}
  std::vector<std::string> in_commod;

  #pragma cyclus var {"tooltip": "output commodity", \
                      "uitype": ["oneOrMore", "outcommodity"]}
  std::vector<std::string> commod_out;

  #pragma cyclus var {"tooltip": "Efficiency input file"}
  std::string repro_input_path;

  #pragma cyclus var {"default": 1e299, \
                      "tooltip": "reactor maximum inventory size", \
                      "doc": "total maximum inventory size of the reactor"}
  double max_inv_size;

  #pragma cyclus var {"default": 150.1, \
                      "tooltip": "Conversion efficiency."}
  double input_capacity;

  #pragma cyclus var {"default": 10, \
                      "tooltip": "Total conversion capacity."}
  double output_capacity;

  #pragma cyclus var
  std::vector<cyclus::toolkit::ResourceBuff> out_inventory;

  #pragma cyclus var {"capacity": "max_inv_size"}
  cyclus::toolkit::ResourceBuff input_inventory;

  #pragma cyclus var {"capacity": "max_inv_size"}
  cyclus::toolkit::ResourceBuff waste_inventory;

};

}  // namespace reactor

#endif  // CYCLUS_REPROCESSES_REPROCESS_FACILITY_H_
