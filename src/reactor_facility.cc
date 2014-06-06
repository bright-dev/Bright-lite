#include "reactor_facility.h"

namespace reactor {

ReactorFacility::ReactorFacility(cyclus::Context* ctx)
    : cyclus::Facility(ctx) {};

std::string ReactorFacility::str() {
  return Facility::str();
}

void ReactorFacility::Tick() {}

void ReactorFacility::Tock() {}

std::set<cyclus::RequestPortfolio<cyclus::Material>::Ptr> ReactorFacility::GetMatlRequests() {
  using cyclus::RequestPortfolio;
  using cyclus::Material;
  using cyclus::Composition;
  using cyclus::CompMap;
  using cyclus::CapacityConstraint;

  CompMap cm;
  Material::Ptr target = Material::CreateUntracked(core_mass/batches, 
                          Composition::CreateFromAtom(cm));
  CapacityConstraint<Material> cc(1.0);

  RequestPortfolio<Material>::Ptr port(new RequestPortfolio<Material>());
  port->AddRequest(target, this, in_commod);
  port->AddConstraint(cc);

  std::set<RequestPortfolio<Material>::Ptr> ports;
  ports.insert(port);
  return ports;
}

extern "C" cyclus::Agent* ConstructReactorFacility(cyclus::Context* ctx) {
  return new ReactorFacility(ctx);
}

}
