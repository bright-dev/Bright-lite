#include "reactor_facility.h"

namespace reactor {

ReactorFacility::ReactorFacility(cyclus::Context* ctx)
    : cyclus::Facility(ctx) {
      cycle_end_ = ctx->time();
      start_time_ = cycle_end_;
};

std::string ReactorFacility::str() {
  return Facility::str();
}

void ReactorFacility::Tick() {
    //std::cout << "tick begin, inventory size: " << inventory.count() << std::endl;

    // if the reactor has just been deployed
    if(fuel_library_.name.size() == 0){
        std::ifstream inf(libraries[0] +"/manifest.txt"); //opens manifest file
        std::string line;
        std::string iso_name;
        fuel_library_.name = libraries[0]; //for now only one entry in here
        while(getline(inf, line)){
            std::istringstream iss(line);
            isoInformation iso; //creates temporary iso
            iss >> iso_name;
            iso.name = atoi(iso_name.c_str()); //adds name to iso
            iso.region = 0; //default region is the zeroth region
            iso.fraction = 0; //zero fraction
            fuel_library_.all_iso.push_back(iso); //adds this temp iso to fuel_library
        }
        if (libraries.size() < 2){
        //adds a built isoInfo for each isotope that can be in fuel, stores it in all_iso
            DataReader2(fuel_library_.name, fuel_library_.all_iso);
        } else {
            for(int i = 0; i < libraries.size(); i++){
                fuel_library_.interpol_libs.push_back(libraries[i]);
            }
            for(std::map<std::string, double>::iterator c = interpol_pairs.begin(); c != interpol_pairs.end(); ++c){
                interpol_pair pair;
                pair.metric = c -> first;
                pair.value = c -> second;
                fuel_library_.interpol_pairs.push_back(pair);
            }
            fuel_library_ = lib_interpol(fuel_library_);
            //mass_check(fuel_library_);
        }
        //adds general info about the fuel in fuel_library_
        ///if theres value, dont update field
        fuel_library_.name = libraries[0];
        fuel_library_.base_flux = flux_finder(fuel_library_.name);
        fuel_library_.base_mass = core_mass;
        fuel_library_.base_power = generated_power;
        fuel_library_.pnl = nonleakage;
        fuel_library_.target_BU = target_burnup;
        fuel_library_.fuel_area = fuel_area;
        fuel_library_.cylindrical_delta = cylindrical_delta;
        fuel_library_.mod_Sig_a = mod_Sig_a;
        fuel_library_.mod_Sig_tr = mod_Sig_tr;
        fuel_library_.mod_Sig_f = mod_Sig_f;
        fuel_library_.mod_thickness = mod_thickness;
        fuel_library_.fuel_Sig_tr = fuel_Sig_tr;
        fuel_library_.disadv_a = disadv_a;
        fuel_library_.disadv_b = disadv_b;
        fuel_library_.disadv_mod_siga = disadv_mod_siga;
        fuel_library_.disadv_mod_sigs = disadv_mod_sigs;
        fuel_library_.disadv_fuel_sigs = disadv_fuel_sigs;


        batch_info empty_batch;
        //creates empty batch entries
        for(int i = 0; i < batches; i++){
            fuel_library_.batch.push_back(empty_batch);
            fuel_library_.batch[i].batch_fluence = 0; //fluence of the batch is set to zero
            fuel_library_.batch[i].fraction.push_back(1);//fraction of the batch is set to one
        }

        /*
        core_ = std::vector<fuelBundle>(batches); //core_ is size of batches, which will be changed
        for(int i = 0; i < core_.size(); i++){
          core_[i] = fuel_library_[i];
        }*/

/************************output file*********************************/
        std::ofstream outfile("../output_cyclus_recent.txt");

        outfile << "Fuel library name: " << fuel_library_.name << "\r\n";
        outfile << "Batches: " << batches << "\r\n";
        outfile << "libcheck: " << fuel_library_.libcheck << "\r\n";
        outfile << "Leakage: " << fuel_library_.pnl << "\r\n";
        outfile << "Burnup calculation timestep: " << burnupcalc_timestep << " [days]\r\n";
        outfile << "Flux calculation method: " << flux_mode << "\r\n";
        if(flux_mode == 3){//if the flux calculation is the cylindrical method
        outfile << "   Fuel area: " << fuel_library_.fuel_area << " [cm2]\r\n";
        outfile << "   Fuel Sig_tr: " << fuel_library_.fuel_Sig_tr << " [cm-1]\r\n";
        outfile << "   Delta: " << fuel_library_.cylindrical_delta << " [cm]\r\n";
        outfile << "   Moderator thickness: " << fuel_library_.mod_thickness << " [cm]\r\n";
        outfile << "   Moderator Sig_a: " << fuel_library_.mod_Sig_a << " [cm-1]\r\n";
        outfile << "   Moderator Sig_tr: " << fuel_library_.mod_Sig_tr << " [cm-1]\r\n";
        outfile << "   Moderator Sig_f: " << fuel_library_.mod_Sig_f << " [cm-1]\r\n";
        }
        outfile << "\r\n|Thermal disadvantage calculation:\r\n";
        outfile << "| Fuel radius: " << fuel_library_.disadv_a << " [cm]\r\n";
        outfile << "| Moderator radius: " << fuel_library_.disadv_b << " [cm]\r\n";
        outfile << "| Moderator Sigma_a: " << fuel_library_.disadv_mod_siga << " [cm-1]\r\n";
        outfile << "| Moderator Sigma_s: " << fuel_library_.disadv_mod_sigs << " [cm-1]\r\n";
        outfile << "| Fuel Sigma_s: " << fuel_library_.disadv_fuel_sigs << " [cm-1]\r\n\r\n";
        outfile << "Base flux: " << fuel_library_.base_flux << "\r\n";
        outfile << "Base power: " << fuel_library_.base_power << "\r\n";
        outfile << "Base mass: " << fuel_library_.base_mass << "\r\n\r\n\r\n";

        outfile.close();
/************************End of output file*********************************/
    }

    //check if its time to decommission the facility
    if(start_time_ + reactor_life <= cycle_end_){
        //std::cout << "here: " << start_time_ + reactor_life << "  " << cycle_end_ << std::endl;

        //take care of the fuel





    }



    //std::cout << "end tick" << std::endl;
}

void ReactorFacility::Tock() {
    //std::cout << "Begin tock\n" << std::endl;
    if(inventory.count() == 0){return;}
    cyclus::Context* ctx = context();
    if (ctx->time() != cycle_end_) {
        //std::cout << "time: "<< ctx->time()<< "  not end of cycle.  End of cycle: " << cycle_end_ << std::endl;/// <--------
        return;
    }

    // Pop materials out of inventory
    std::vector<cyclus::Material::Ptr> manifest;
    manifest = cyclus::ResCast<cyclus::Material>(inventory.PopN(inventory.count()));

    cyclus::CompMap comp;
    cyclus::CompMap::iterator it;
    if(manifest.size() > fuel_library_.batch.size()){
        for(int i = 0; i < manifest.size() - fuel_library_.batch.size(); i++){
            batch_info temp_batch;
            temp_batch.batch_fluence = 0;
            fuel_library_.batch.push_back(temp_batch);
        }
    }
    //each batch
    for(int i = 0; i < manifest.size(); i++){
    //build correct isoinfo and fraction for each batch
    ///put the stuff in comp fraction to fuel_library_.batch.iso fraction using values in fuel_library_.all_iso

        comp = manifest[i]->comp()->mass(); //store the fractions of i'th batch in comp
        int comp_iso;
        //each iso in comp
        for (it = comp.begin(); it != comp.end(); ++it){
            comp_iso = pyne::nucname::zzaaam(it->first);
            //each iso in all_iso
            for(int j = 0; j < fuel_library_.all_iso.size(); j++){
                int fl_iso = fuel_library_.all_iso[j].name;
                if(fl_iso == comp_iso && fuel_library_.batch[i].batch_fluence == 0){
                    //std::cout << "i: " << i << "  " << fl_iso << "  " << comp_iso << std::endl;
                    isoInformation temp_iso;
                    temp_iso = fuel_library_.all_iso[j];
                    temp_iso.fraction = it->second;
                    fuel_library_.batch[i].iso.push_back(temp_iso);
                }
            }

        }
    }
    //collapse iso's, read struct effects, reorder the fuel batches accordingly
    batch_reorder();
  // pass fuel bundles to burn-up calc
  fuel_library_ = burnupcalc(fuel_library_, flux_mode, DA_mode, burnupcalc_timestep);
  //std::cout << "after burnupcalc" << std::endl;
  // convert fuel bundle into materials
  for(int i = 0; i < fuel_library_.batch.size(); i++){
    cyclus::CompMap out_comp;
    for(std::map<int, double>::iterator c = fuel_library_.batch[i].comp.begin(); c != fuel_library_.batch[i].comp.end(); ++c){
      if(c->second < 0){
        out_comp[pyne::nucname::zzaaam_to_id(c->first)] = 0;
      } else {
        out_comp[pyne::nucname::zzaaam_to_id(c->first)] = c->second;
      }
    }
    manifest[i]->Transmute(cyclus::Composition::CreateFromMass(out_comp));
    inventory.Push(manifest[i]);
  }

  // cycle end update
  cycle_end_ = ctx->time() + ceil(fuel_library_.batch[fuel_library_.batch.size()-1].batch_fluence/(86400*fuel_library_.base_flux*28));
  //std::cout << "Cycle length: " << ceil(fuel_library_.batch[fuel_library_.batch.size()-1].batch_fluence/(86400*fuel_library_.base_flux*28)) << std::endl;
  //std::cout << "Time :: " <<cycle_end_ << std::endl;

    /************************output file*********************************/
    std::ofstream outfile;
    outfile.open("../output_cyclus_recent.txt", std::ios::app);

    outfile << " Cycle length: " << ceil(fuel_library_.batch[fuel_library_.batch.size()-1].batch_fluence/(86400*fuel_library_.base_flux*28)) << " [months]";

    outfile << "\r\n\r\n\r\n";

    outfile.close();
    /************************End of output file**************************/

}

std::set<cyclus::RequestPortfolio<cyclus::Material>::Ptr> ReactorFacility::GetMatlRequests() {
  //std::cout << "Getmatlrequests begin" << std::endl;
  using cyclus::RequestPortfolio;
  using cyclus::Material;
  using cyclus::Composition;
  using cyclus::CompMap;
  using cyclus::CapacityConstraint;
  std::set<RequestPortfolio<Material>::Ptr> ports;
  cyclus::Context* ctx = context();
  if (ctx->time() != cycle_end_ && inventory.count() != 0){
    return ports;
  }

  CompMap cm;
  Material::Ptr target = Material::CreateUntracked(core_mass/batches,
                          Composition::CreateFromAtom(cm));

  RequestPortfolio<Material>::Ptr port(new RequestPortfolio<Material>());
  double qty;

  if(inventory.count() == 0){
     for(int i = 0; i < batches; i++){
        //checks to see if there is a next in_commod to request, otherwise puts the first commod request
        if(in_commods.size() > i+1){
            port->AddRequest(target, this, in_commods[i+1]);
        } else {
            port->AddRequest(target, this, in_commods[0]);
        }
     }

     qty = core_mass;
  } else {
     port->AddRequest(target, this, in_commods[0]);
     qty = core_mass/batches;
  }
  CapacityConstraint<Material> cc(qty);

  port->AddConstraint(cc);
  ports.insert(port);
  //std::cout << "end getmatlrequests" << std::endl;
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
  cyclus::Context* ctx = context();
  std::set<BidPortfolio<Material>::Ptr> ports;
  if (ctx->time() != cycle_end_){
    return ports;
  }
  // respond to all requests of my commodity
  if (inventory.count() == 0){
    //std::cout << "YAYa8?" << std::endl; 
    return ports;}
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
  CapacityConstraint<Material> cc(core_mass/batches);

  port->AddConstraint(cc);
  ports.insert(port);
  //std::cout << "end getmatlbids" << std::endl;
  return ports;
}

void ReactorFacility::AcceptMatlTrades(const std::vector< std::pair<cyclus::Trade<cyclus::Material>, cyclus::Material::Ptr> >& responses) {
    //std::cout << "begin accptmatltrades" << std::endl;
    std::vector<std::pair<cyclus::Trade<cyclus::Material>, cyclus::Material::Ptr> >::const_iterator it;
    cyclus::Composition::Ptr compost;
    //std::cout << "YAYa9?" << std::endl;
    for (it = responses.begin(); it != responses.end(); ++it) {

        inventory.Push(it->second);
        compost = it->second->comp();
        cyclus::CompMap cmap = compost->mass();
        cyclus::CompMap::iterator cit;

/************************output file*********************************/
        std::ofstream outfile;
        outfile.open("../output_cyclus_recent.txt", std::ios::app);
        outfile << "Composition of fresh batch:\r\n";

        for(cit = cmap.begin(); cit != cmap.end(); ++cit){
            outfile << "  Isotope: " << cit->first << "  Fraction: " << cit->second << "\r\n";
        }
        outfile << "\r\n\n";

        outfile.close();
/************************End of output file**************************/
  }
  //std::cout << "end acceptmatltrades" << std::endl;
}

void ReactorFacility::GetMatlTrades(const std::vector< cyclus::Trade<cyclus::Material> >& trades,
    std::vector<std::pair<cyclus::Trade<cyclus::Material>,cyclus::Material::Ptr> >& responses) {
    using cyclus::Material;
    using cyclus::Trade;
    std::cout << "YAYa10?" << std::endl;
    //std::cout << "begin getmatltrades" << std::endl;
    std::vector< cyclus::Trade<cyclus::Material> >::const_iterator it;
    cyclus::Material::Ptr discharge = cyclus::ResCast<Material>(inventory.Pop());
    fuel_library_.batch.erase(fuel_library_.batch.begin());
    for (it = trades.begin(); it != trades.end(); ++it) {
        responses.push_back(std::make_pair(*it, discharge));
    }
    //std::cout << "end getmatltrades" << std::endl;
}


double ReactorFacility::burnup_test(cyclus::Material::Ptr new_batch ){
/*
    fuelBundle bundle;
    bundle = fuel_library_;
    cyclus::CompMap comp;
    cyclus::CompMap::iterator it;
    comp = new_batch->comp()->mass();
    int j = 0;
    int fl_iso = bundle.all_iso[j].name;
    int comp_iso;
    for (it = comp.begin(); it != comp.end(); ++it){
        comp_iso = pyne::nucname::zzaaam(it->first);
        if(fl_iso < comp_iso) {
            while(fl_iso < comp_iso && j < bundle.all_iso.size()-1){
                fl_iso = bundle.all_iso[j++].name;
            }
        }
        if(fl_iso == comp_iso){
            if (bundle.batch_fluence == 0){
                bundle.all_iso[j].fraction[0] = it->second;
            }
            j+=1;
            fl_iso = bundle.all_iso[j].name;
        }
    }
    std::vector<fuelBundle> temp_core;
    for(int i = 1; i < core_.size(); i++){
        temp_core.push_back(core_[i]);
    }
    temp_core.push_back(bundle);
    std::vector<fuelInfo> reactor_return;
    reactor_return = burnupcalc(core_, nonleakage, 0.0001);
    return reactor_return[reactor_return.size()].burnup;
*/
}

void ReactorFacility::start_up(){
    std::vector<double> target_burnups;
    for(int i = 0; i < batches; i++){
        target_burnups.push_back(target_burnup*i/batches);
    }
    for(int i = 0; i < batches; i++){
        std::pair<double, std::pair<double, std::map<int, double> > > test;
        //test = blending_calc(fuel_library_, target_BUd, flux_mode, DA_mode, burnupcalc_timestep));
        //std::cout << test.first << std::endl;
    }

}


void ReactorFacility::batch_reorder(){
//collapses each batch first, thgien orders them
    std::cout << "Begin batch_reorder" << std::endl;
    double k0, k1;
    fuel_library_ = StructReader(fuel_library_);
    fuel_library_ = regionCollapse(fuel_library_);
    bool test = false;
    for(int i = 0; i < fuel_library_.batch.size(); i++){
        if(fuel_library_.batch[i].batch_fluence != 0){
            test = true;
        }
    }
    if(test == true){return;}
    std::cout << "BEGIN SS_BURNUPCALC" << std::endl;
    
    double burnup = SS_burnupcalc(fuel_library_.batch[0].collapsed_iso, 3, 30, 0.98, fuel_library_.base_flux);
    
    std::cout << "Result: " << burnup << std::endl;
    std::cout << "END SS_BURNUPCALC" << std::endl;

    fuelBundle temp_fuel = fuel_library_;
    fuel_library_.batch.clear();

    //goes through all the batches and orders them from lowest k to highest
    //uses temp_fuel to store the unordered batches
    while(temp_fuel.batch.size() != 0){
        int lowest = 0;
        for(int i = 1; i < temp_fuel.batch.size(); i++){
            k0 = temp_fuel.batch[lowest].collapsed_iso.neutron_prod[1]/temp_fuel.batch[lowest].collapsed_iso.neutron_dest[1];
            k1 = temp_fuel.batch[i].collapsed_iso.neutron_prod[1]/temp_fuel.batch[i].collapsed_iso.neutron_dest[1];
            if(k1 < k0){
                //if k of i'th batch is smaller than k of current lowest batch
                lowest = i;
            }
        }
        fuel_library_.batch.push_back(temp_fuel.batch[lowest]);
        temp_fuel.batch.erase(temp_fuel.batch.begin() + lowest);
    }
    //std::cout << " End batch_reorder" << std::endl;
}


extern "C" cyclus::Agent* ConstructReactorFacility(cyclus::Context* ctx) {
  return new ReactorFacility(ctx);
}

}
