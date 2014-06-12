#include "reactor_facility.h"

namespace reactor {

ReactorFacility::ReactorFacility(cyclus::Context* ctx)
    : cyclus::Facility(ctx) {
  cycle_end_ = ctx->time();
};

std::string ReactorFacility::str() {
  return Facility::str();
}

void ReactorFacility::Tick() {
  std::cout << "TICK" << std::endl;
  if(fuel_library_.name.size() == 0){
    std::ifstream inf(libraries[0] +"/manifest.txt");
    std::string line;
    std::string iso_name;
    fuel_library_.name = libraries[0];
    while(getline(inf, line)){
      std::istringstream iss(line);
      isoInformation iso;
      iss >> iso_name;
      iso.name = atoi(iso_name.c_str());
      iso.region = 0;
      iso.fraction.push_back(0);
      fuel_library_.iso.push_back(iso);
    }
    if (libraries.size() == 1){
      DataReader2(fuel_library_.name, fuel_library_.iso);
    } else {
    ///interpolation stuff
    }
    core_ = std::vector<fuelBundle>(batches);
    for(int i = 0; i < core_.size(); i++){
      core_[i] = fuel_library_;
      core_[i].batch_fluence = 0;
    } 
  }
}

void ReactorFacility::Tock() {
  std::cout << "TOCK" << std::endl;
  cyclus::Context* ctx = context();
  if (ctx->time() != cycle_end_)
    return;
  // Pop materials out of inventory
  std::vector<cyclus::Material::Ptr> manifest;
  std::cout << "inventory:   "<< inventory.count() << std::endl;
  manifest = cyclus::ResCast<cyclus::Material>(inventory.PopN(inventory.count()));
  // convert materials into fuel bundles
  cyclus::CompMap comp;
  cyclus::CompMap::iterator it;
  if( core_.size() < batches){
    fuelBundle bundle;
    core_.push_back(bundle);
  } 
  core_[batches-1].batch_fluence = 0;
  std::cout << "MANIFEST:   " << manifest.size() << std::endl;
  for(int i = 0; i < manifest.size(); i++){
     comp = manifest[i]->comp()->mass();
     int j = 0;
     int fl_iso = core_[i].iso[j].name;
     int comp_iso;
     for (it = comp.begin(); it != comp.end(); ++it){
       comp_iso = pyne::nucname::zzaaam(it->first);
       if(fl_iso < comp_iso) {
         while(fl_iso < comp_iso && j < core_[i].iso.size()-1){
           fl_iso = core_[i].iso[j++].name;        
         }
       }
       if(fl_iso == comp_iso){
         core_[i].iso[j].fraction[0] = it->second;
         j+=1;
         fl_iso = core_[i].iso[j].name;
       } 
    }
  }
  /// pass fuel bundles to burn-up calc
  std::map<double, std::map<int, double> > reactor_return;
  reactor_return = burnupcalc(core_, nonleakage, 0.0001); 
  /// convert fuel bundles into materials
  int i = 0;
  for(std::map<double, std::map<int, double> >::iterator rr = reactor_return.begin(); rr != reactor_return.end(); ++rr){
    cyclus::CompMap out_comp;
    for(std::map<int, double>::iterator c = rr->second.begin(); c != rr->second.end(); ++c){
      out_comp[pyne::nucname::zzaaam_to_id(c->first)] = c->second;
    }
    manifest[i]->Transmute(cyclus::Composition::CreateFromMass(out_comp));
    inventory.Push(manifest[i]);
    ++i;
  }
  // cycle end update
  cycle_end_ = ctx->time() + ceil(reactor_return.rbegin()->first/28);
}

std::set<cyclus::RequestPortfolio<cyclus::Material>::Ptr> ReactorFacility::GetMatlRequests() {
  using cyclus::RequestPortfolio;
  using cyclus::Material;
  using cyclus::Composition;
  using cyclus::CompMap;
  using cyclus::CapacityConstraint;
  std::cout << "Requests" << std::endl;
  std::set<RequestPortfolio<Material>::Ptr> ports;
  cyclus::Context* ctx = context();
  if (ctx->time() != cycle_end_)
    return ports;
  CompMap cm;
  Material::Ptr target = Material::CreateUntracked(core_mass/batches, 
                          Composition::CreateFromAtom(cm));

  RequestPortfolio<Material>::Ptr port(new RequestPortfolio<Material>());
  port->AddRequest(target, this, in_commod);
  double qty;
  if(inventory.count() == 0){
     for(int i = 0; i < batches - 1; i++){
        port->AddRequest(target, this, in_commod);
     }
     qty = core_mass;
  } else {
     qty = core_mass/batches;
  }
  CapacityConstraint<Material> cc(qty);

  port->AddConstraint(cc);
  ports.insert(port);
  return ports;
}

// MatlBids //
std::set<cyclus::BidPortfolio<cyclus::Material>::Ptr>
  ReactorFacility::GetMatlBids(
    cyclus::CommodMap<cyclus::Material>::type& commod_requests) {
  using cyclus::BidPortfolio;
  using cyclus::CapacityConstraint;
  using cyclus::Converter;
  using cyclus::Material;
  using cyclus::Request;

  std::cout << "Bids" << std::endl;
  cyclus::Context* ctx = context();
  std::set<BidPortfolio<Material>::Ptr> ports;
  if (ctx->time() != cycle_end_)
    return ports;
  // respond to all requests of my commodity
  if (inventory.count() == 0){return ports;}

  std::vector<cyclus::Material::Ptr> manifest;
  manifest = cyclus::ResCast<Material>(inventory.PopN(inventory.count()));

  BidPortfolio<Material>::Ptr port(new BidPortfolio<Material>());
  std::vector<Request<Material>*>& requests = commod_requests[out_commod];
  std::vector<Request<Material>*>::iterator it;
  for (it = requests.begin(); it != requests.end(); ++it) {
    Request<Material>* req = *it;
    if (req->commodity() == out_commod) {
      Material::Ptr offer = Material::CreateUntracked(core_mass/batches, manifest[0]->comp());
      port->AddBid(req, offer, this);
    }
  }
  inventory.PushAll(manifest);

  
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

void ReactorFacility::GetMatlTrades(
    const std::vector< cyclus::Trade<cyclus::Material> >& trades,
    std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                          cyclus::Material::Ptr> >& responses) {
  using cyclus::Material;
  using cyclus::Trade;

  std::vector< cyclus::Trade<cyclus::Material> >::const_iterator it;
  for (it = trades.begin(); it != trades.end(); ++it) {
    responses.push_back(std::make_pair(*it, cyclus::ResCast<Material>(inventory.Pop())));
  }
}
extern "C" cyclus::Agent* ConstructReactorFacility(cyclus::Context* ctx) {
  return new ReactorFacility(ctx);
}

}
