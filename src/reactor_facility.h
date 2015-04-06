#ifndef CYCLUS_REACTORS_REACTOR_FACILITY_H_
#define CYCLUS_REACTORS_REACTOR_FACILITY_H_

#include <string>
#include <sstream>

#include "cyclus.h"
#include "burnupcalc.h"

namespace reactor {

/// @class ReactorFacility
///
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
/// The facility starts by reading the available isotope library databases.
/// A database is specific to the reactor and fuel type, therefore one reactor
/// may have more than one database describing it in cases where specific database
/// is not available, requiring an interpolation
/// to calculate a single isotope library database.
///
/// During the first tick the reactor performs the database interpolation (if
/// necessary) and builds the base parameters from the input file to fuel_library_.
/// All isotopes available for input fuel is placed in fuel_library_.all_iso
/// with correct names and isoinformation. This is independent of batches.
/// At the end of the first tick fuel_library_ has correct number of batches in
/// fuel_library_.batch[] all with zero fluence and no other information.
///
/// During the first material exchange the reactor receives fuel. The composition
/// of fuel is stored in resourcebuff inventory. Fuel may be blended, or a
/// composition may be forced to the reactor. The upstream facilities and
/// target burnup determines how the reactor will receive fuel. If more than one
/// in_commods is defined, then in_commods[0] will be used for refueling and
/// the others will be used for startup. If there are less in_commods than
/// there are batches to define the startup, in_commods[0] will used as a default.
///
/// At the beginning of first tock the reactor inventory is used to update the
/// composition of each batch in fuel_library_.iso[].fraction. Next, batch_reorder
/// is called if reactor is in startup. In this case the batcg libraries are built
/// for all batches the the core ordered.




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
                              "for the design of new facility agents.", \
                        "niche": "reactor"}

  /// A verbose printer for the ReactorFacility
  virtual std::string str();

  /// The handleTick function specific to the ReactorFacility.
  /// @param time the time of the tick
  virtual void Tick();

  /// The handleTick function specific to the ReactorFacility.
  /// @param time the time of the tock
  virtual void Tock();

  virtual void Build(cyclus::Agent* parent) {
    Facility::Build(parent);
    if (lifetime() >= 0) {
      context()->SchedDecom(this, exit_time());
    }
  }

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

  fuelBundle comp_function(cyclus::Material::Ptr mat1, fuelBundle &fuel_library_);

  fuelBundle comp_trans(cyclus::Material::Ptr mat1, fuelBundle &fuel_library_);

  double blend_next(cyclus::toolkit::ResourceBuff fissle,
                                   cyclus::toolkit::ResourceBuff non_fissle,
                                   std::vector<cyclus::toolkit::ResourceBuff> inventory,
                                   std::map<std::string, double> incommods);

  double start_up(cyclus::toolkit::ResourceBuff fissle,
                                   cyclus::toolkit::ResourceBuff non_fissle,
                                   std::vector<cyclus::toolkit::ResourceBuff> inventory,
                                   std::map<std::string, double> incommods);
  int refuels;

  void CoreBuilder();
  void batch_reorder();

  double SS_enrich;
  double ss_fraction;

  #pragma cyclus var {"default": 0.01, \
                      "userlevel": 3, \
                     "tooltip": "The blending convergence requirement for the code"}
  double tolerance;

  #pragma cyclus var {"default": 0.01, \
                      "userlevel": 3, \
                     "tooltip": "The steady state calculation tolerance for the code"}
  double SS_tolerance;

  #pragma cyclus var {"default": 0.4095, \
                     "units": "cm", \
                      "userlevel": 3, \
                     "tooltip": "Fuel pin radius used for thermal disadvantage calculation."}
  double disadv_a;

  #pragma cyclus var {"default": 0.70749, \
                     "units": "cm", \
                      "userlevel": 3, \
                     "tooltip": "Moderator radius used for thermal disadvantage calculation."}
  double disadv_b;

  #pragma cyclus var {"default": 0.222, \
                     "units": "cm-1", \
                      "userlevel": 3, \
                     "tooltip": "Moderator macroscopic abs cross-section for thermal disadvantage calculation."}
  double disadv_mod_siga;

  #pragma cyclus var {"default": 3.44, \
                     "units": "cm-1", \
                      "userlevel": 3, \
                     "tooltip": "Moderator macroscopic scattering cross-section for thermal disadvantage calculation."}
  double disadv_mod_sigs;

  #pragma cyclus var {"default": 0.43, \
                     "units": "cm-1", \
                      "userlevel": 3, \
                     "tooltip": "Fuel macroscopic scattering cross-section for thermal disadvantage calculation."}
  double disadv_fuel_sigs;

  #pragma cyclus var {"tooltip": "input commodity", \
                      "doc": "commodity that the brightlite reactor consumes", \
                      "schematype": "token", \
                      "uitype": ["oneOrMore", "incommodity"]}
  std::vector<std::string> in_commods;

  #pragma cyclus var {"tooltip": "output commodity", \
                      "doc": "commodity that the brightlite supplies", \
                      "uitype": "outcommodity", \
                      "uilabel": "Output"}
  std::string out_commod;

  #pragma cyclus var {"tooltip": "reactor libraries to load", \
                      "userlevel": 0, \
                      "doc": "the reactor's burnup & criticality behavior to use"}
  std::vector<std::string> libraries;

  #pragma cyclus var {"tooltip": ["interpolation pairs used for the library", \
                      "Interpolation metric", "Interpolation values"], \
                      "default": {"BURNUP": 42}, \
                      "uitype": ["oneOrMore", "string", "double"], \
                      "userlevel": 2}
  std::map<std::string, double> interpol_pairs;

  #pragma cyclus var {"tooltip": "number of batches", \
                      "userlevel": 1, \
                      "default": 3}
  int batches;

  #pragma cyclus var {"tooltip": "Non-leakage probability", \
                      "doc": "Varies from 0 - 1.", \
                      "userlevel": 1, \
                      "default": 0.98}
  double nonleakage;

  #pragma cyclus var {"default": 1e299, \
                      "userlevel": 2, \
                      "tooltip": "reactor maximum inventory size", \
                      "doc": "total maximum inventory size of the reactor"}
  double max_inv_size;

  #pragma cyclus var {"tooltip": "Target burnup", \
                      "default": 0, \
                      "userlevel": 0, \
                      "units": "MWd/kgIHM"}
  double target_burnup;

  #pragma cyclus var {"default": 1000.0, \
                      "units": "MWe", \
                      "userlevel": 2, \
                      "tooltip": "Electrical production."}
  double generated_power;

  #pragma cyclus var {"units": "kgIHM", \
                      "userlevel": 0, \
                      "tooltip": "Total mass of the core."}
  double core_mass;

  #pragma cyclus var {"default": 0.33, \
                      "userlevel": 2, \
                      "tooltip": "Thermal to electric conversion rate."}
  double efficiency;

  #pragma cyclus var {"default": 4197, \
                      "units": "cm2", \
                      "userlevel": 3, \
                      "tooltip": "Total area of the fuel in core. Used for cylindrical flux calculation."}
  double fuel_area;

  #pragma cyclus var {"default": 5, \
                      "userlevel": 3, \
                      "tooltip": "Delta to be used for cylindrical flux calculation."}
  double cylindrical_delta;

  #pragma cyclus var {"default": 0.0222, \
                      "units": "cm-1", \
                      "userlevel": 3, \
                      "tooltip": "Macroscopic absorption cross section of the moderator."}
  double mod_Sig_a;

  #pragma cyclus var {"default": 3.46, \
                      "units": "cm-1", \
                      "userlevel": 3, \
                      "tooltip": "Macroscopic transport cross section of the moderator."}
  double mod_Sig_tr;

  #pragma cyclus var {"default": 0.0, \
                      "units": "cm-1", \
                      "userlevel": 3, \
                      "tooltip": "Macroscopic fission cross section of the moderator."}
  double mod_Sig_f;

  #pragma cyclus var {"default": 50, \
                      "units": "cm", \
                      "userlevel": 3, \
                      "tooltip": "Radial thickness of the moderator used for cylindrical flux calculation."}
  double mod_thickness;

  #pragma cyclus var {"default": 0.12, \
                      "units": "cm-1", \
                      "userlevel": 3, \
                      "tooltip": "Macroscopic transport cross section of the fuel."}
  double fuel_Sig_tr;

  #pragma cyclus var {"default": 10, \
                      "userlevel": 2, \
                      "tooltip": "Timestep [days] for the burnup calculation."}
  double burnupcalc_timestep;

  #pragma cyclus var {"default": 1, \
                      "userlevel": 2, \
                      "tooltip": "Flux calculation method. 1:Uniform, 2:Inv.Neut.Prod, 3:Cylindrical"}
  int flux_mode;

  #pragma cyclus var {"default": 0, \
                      "userlevel": 2, \
                      "tooltip": "Disadvantage calculation. 0:Off, 1:On"}
  int DA_mode;

  #pragma cyclus var {"capacity": "max_inv_size"}
  cyclus::toolkit::ResourceBuff inventory;

  #pragma cyclus var {"default": 480, \
                      "units": "months", \
                      "userlevel": 1, \
                      "tooltip": "Time before reactor is shutdown after startup."}
  int reactor_life;

  #pragma cyclus var {"units": "NUCID", \
                      "userlevel": 3, \
                      "tooltip": "List of fissile isotopes for conversion ratio calculation."}
  std::vector<std::string> CR_fissile;

  #pragma cyclus var {"userlevel": 3, \
                      "default": 70, \
                      "tooltip": "Lower bound of the fission product mass number for conversion ratio calculation."}
  int CR_lower;

  #pragma cyclus var {"userlevel": 3, \
                      "default": 160, \
                      "tooltip": "Upper bound of the fission product mass number for conversion ratio calculation."}
  int CR_upper;

  #pragma cyclus var {"userlevel": 3, \
                      "default": -1, \
                      "tooltip": "Target for conversion ratio."}
  double CR_target;

  #pragma cyclus var {"userlevel": 2, \
                      "default": 28, \
                      "tooltip": "The amount of time the reactor is offline to load new fuel and shuffle older fuel", \
                      "label": "Outage Period"}
  double outage_time;


 private:
  bool shutdown;
  bool record;
  int cycle_end_;
  int start_time_;
  fuelBundle fuel_library_;
  fuelBundle core_;
  double ss_fluence = 0;
  double p_time = 0;
  double p_frac = 0;
  int outage_shutdown = 0;
  int steady_state = 0;

};

}  // namespace reactor

#endif  // CYCLUS_REACTORS_REACTOR_FACILITY_H_
