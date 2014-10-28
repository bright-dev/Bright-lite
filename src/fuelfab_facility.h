#ifndef CYCLUS_FUELFABS_FUELFAB_FACILITY_H_
#define CYCLUS_FUELFABS_FUELFAB_FACILITY_H_

#include <string>

#include "cyclus.h"
#include "reactor_facility.h"

namespace fuelfab {


/// @class FuelfabFacility
///
/// This Facility is intended
/// as a skeleton to guide the implementation of new Facility
/// agents.
/// The FuelfabFacility class inherits from the Facility class and is
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
class FuelfabFacility : public cyclus::Facility  {
 public:
  /// Constructor for FuelfabFacility Class
  /// @param ctx the cyclus context for access to simulation-wide parameters
  explicit FuelfabFacility(cyclus::Context* ctx);

  /// The Prime Directive
  /// Generates code that handles all input file reading and restart operations
  /// (e.g., reading from the database, instantiating a new object, etc.).
  /// @warning The Prime Directive must have a space before it! (A fix will be
  /// in 2.0 ^TM)

  #pragma cyclus

  #pragma cyclus note {"doc": "A fuelfab facility is provided as a skeleton " \
                              "for the design of new facility agents.",\
		       "niche": "fuel fabrication"}

    /// A verbose printer for the FuelfabFacility
    virtual std::string str();

    /// The handleTick function specific to the FuelfabFacility.
    /// @param time the time of the tick
    virtual void Tick();

    /// The handleTick function specific to the FuelfabFacility.
    /// @param time the time of the tock
    virtual void Tock();

    std::set<cyclus::RequestPortfolio<cyclus::Material>::Ptr> GetMatlRequests();

    std::set<cyclus::BidPortfolio<cyclus::Material>::Ptr> GetMatlBids(
        cyclus::CommodMap<cyclus::Material>::type& commod_requests);

    void AcceptMatlTrades(const std::vector< std::pair<cyclus::Trade<cyclus::Material>,
                                            cyclus::Material::Ptr> >& responses);

    void GetMatlTrades(const std::vector< cyclus::Trade<cyclus::Material> >& trades,
                                        std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                                        cyclus::Material::Ptr> >& responses);

    double start_time;
    // and away we go!
    std::vector<cyclus::toolkit::ResourceBuff> inventory;

    #pragma cyclus var {"tooltip": "maximum storage capacity of the facility", \
                        "default": 1E60}
    double maximum_storage;

    #pragma cyclus var {"tooltip": "input commodities", \
                        "doc": "A list of the commodities that this facility can recieve", \
                        "uitype": ["oneOrMore", "incommodity", " "]}
    std::vector<std::string> in_commods;

    #pragma cyclus var {"tooltip": "output commodity", \
                      "doc": "commodity that the brightlite supplies", \
                      "uitype": "outcommodity"}
    std::string out_commod;

};

}  // namespace fuelfab

#endif  // CYCLUS_FUELFABS_FUELFAB_FACILITY_H_
