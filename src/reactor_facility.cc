#include "reactor_facility.h"

namespace reactor {

ReactorFacility::ReactorFacility(cyclus::Context* ctx)
    : cyclus::Facility(ctx) {
  cycle_end_ = 0;
  
};

std::string ReactorFacility::str() {
  return Facility::str();
}

void ReactorFacility::Tick() {}

void ReactorFacility::Tock() {
  cyclus::Context* ctx = context();
  if (ctx->time() != cycle_end_)
    return;

  cycle_end_ = ctx->time() + 18;
}

std::set<cyclus::RequestPortfolio<cyclus::Material>::Ptr> ReactorFacility::GetMatlRequests() {
  using cyclus::RequestPortfolio;
  using cyclus::Material;
  using cyclus::Composition;
  using cyclus::CompMap;
  using cyclus::CapacityConstraint;

  std::set<RequestPortfolio<Material>::Ptr> ports;
  cyclus::Context* ctx = context();
  if (ctx->time() != cycle_end_)
    return ports;

  CompMap cm;
  Material::Ptr target = Material::CreateUntracked(core_mass/batches, 
                          Composition::CreateFromAtom(cm));
  CapacityConstraint<Material> cc(1.0);

  RequestPortfolio<Material>::Ptr port(new RequestPortfolio<Material>());
  port->AddRequest(target, this, in_commod);
  port->AddConstraint(cc);

  ports.insert(port);
  return ports;
}

void ReactorFacility::AcceptMatlTrades(
      const std::vector< std::pair<cyclus::Trade<cyclus::Material>,
      cyclus::Material::Ptr> >& responses) {
  std::vector<std::pair<cyclus::Trade<cyclus::Material>, 
                        cyclus::Material::Ptr> >::const_iterator it;
  for (it = responses.begin(); it != responses.end(); ++it) {
    inventory.Push(it->second);
  }
}

extern "C" cyclus::Agent* ConstructReactorFacility(cyclus::Context* ctx) {
  return new ReactorFacility(ctx);
}

}
