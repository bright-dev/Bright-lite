#include "fuelfab_facility.h"

namespace fuelfab {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    FuelfabFacility::FuelfabFacility(cyclus::Context* ctx)
        : cyclus::Facility(ctx) {};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    std::string FuelfabFacility::str() {
      return Facility::str();
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void FuelfabFacility::Tick() {
        if(inventory.size() == 0){
            cyclus::toolkit::ResourceBuff resource;
            for(int i = 0; i < in_commods.size(); i++){
                inventory.push_back(resource);
            }
        }//outputs whats inside the facility
        //std::cout << "FuelFab Inventory:" << std::endl;
        for(int i = 0; i < inventory.size(); i++){
            //sets max inventory at startup
            if(inventory[i].count() == 0){
                inventory[i].set_capacity(maximum_storage/inventory.size());
            }

            //std::cout << "  " << i+1 << ": " << inventory[i].quantity() << std::endl;
            if(inventory[i].count() != 0){
                std::vector<cyclus::Material::Ptr> manifest;
                manifest = cyclus::ResCast<cyclus::Material>(inventory[i].PopN(inventory[i].count()));

                cyclus::CompMap comp;
                cyclus::CompMap::iterator it;

                comp = manifest[0]->comp()->mass();

                inventory[i].PushAll(manifest);
            }
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void FuelfabFacility::Tock() {

    }

    std::set<cyclus::RequestPortfolio<cyclus::Material>::Ptr> FuelfabFacility::GetMatlRequests() {
        using cyclus::RequestPortfolio;
        using cyclus::Material;
        using cyclus::Composition;
        using cyclus::CompMap;
        using cyclus::CapacityConstraint;
        std::set<RequestPortfolio<Material>::Ptr> ports;
        cyclus::Context* ctx = context();
        CompMap cm;
        Material::Ptr target = Material::CreateUntracked(maximum_storage/inventory.size(), Composition::CreateFromAtom(cm));
        RequestPortfolio<Material>::Ptr port(new RequestPortfolio<Material>());
        for(int i = 0; i < inventory.size(); i++){
            double qty = inventory[i].space();

            port->AddRequest(target, this, in_commods[i]);

            CapacityConstraint<Material> cc(qty);
            port->AddConstraint(cc);
            ports.insert(port);
        }
        return ports;
    }

    // MatlBids //
    std::set<cyclus::BidPortfolio<cyclus::Material>::Ptr>FuelfabFacility::GetMatlBids(cyclus::CommodMap<cyclus::Material>::type& commod_requests) {
    using cyclus::BidPortfolio;
    using cyclus::CapacityConstraint;
    using cyclus::Converter;
    using cyclus::Material;
    using cyclus::Request;
    using reactor::ReactorFacility;
    cyclus::Context* ctx = context();
    std::set<BidPortfolio<Material>::Ptr> ports;

    // respond to all requests of my commodity
    int inventory_test = 0;
    for(int i = 0; i < inventory.size(); i++){
        if(inventory[i].count() > 0){
            inventory_test += 1;
        }
    }
    CapacityConstraint<Material> cc(1);
    if (inventory_test == 0){return ports;}
    std::vector<Request<Material>*>& requests = commod_requests[out_commod];
    std::vector<Request<Material>*>::iterator it;
    std::map<cyclus::Trader*, int> facility_request;
    Request<Material>* req;
    //Set up a map for facility requestor to number of requests
    for (it = requests.begin(); it != requests.end(); ++it) {
        Request<Material>* req = *it;
        facility_request[req->requester()] += 1;
    }
    //Iterate through requesters
    std::map<cyclus::Trader*, int>::iterator id;
    for(id = facility_request.begin(); id != facility_request.end(); ++id){
        double limit, nlimit;
        //Cast requester as a reactor if possible
        ReactorFacility* reactor = dynamic_cast<ReactorFacility*>(id->first);
        if (!reactor){
           throw cyclus::CastError("No reactor for fuelfab facility.");
        } else {
            // Check for reactor start up
            if (reactor->inventory.count() == 0){
                limit = reactor->start_up(inventory);
                double temp_limit = limit;
                int k = 1;
                //For each request from reactor create a trade portfolio
                for(it = requests.begin(); it!=requests.end(); ++it){
                    Request<Material>* req = *it;
                    if(req->requester() == id->first){
                        BidPortfolio<Material>::Ptr port(new BidPortfolio<Material>());
                        limit = temp_limit*k/(reactor->batches);
                        double qty = limit + nlimit;
                        qty = qty <= 0 ? 1 : qty;
                        CapacityConstraint<Material> cc(qty);
                        cyclus::Material::Ptr manifest;
                        manifest = cyclus::ResCast<Material>(inventory[0].Pop());
                        Material::Ptr offer = Material::CreateUntracked(limit, manifest->comp());
                        inventory[0].Push(manifest);
                        manifest = cyclus::ResCast<Material>(inventory[1].Pop());
                        offer->Absorb(Material::CreateUntracked(nlimit, manifest->comp()));
                        inventory[1].Push(manifest);
                        port->AddBid(req, offer, this);
                        port->AddConstraint(cc);
                        ports.insert(port);
                        k++;
                    }
                }
                //Same for refueling
            } else{
                for(it = requests.begin(); it!=requests.end(); ++it){
                    Request<Material>* req = *it;
                    if(req->requester() == id->first){
                        BidPortfolio<Material>::Ptr port(new BidPortfolio<Material>());
                        limit = reactor->blend_next(inventory);
                        nlimit = reactor->core_mass/reactor->batches-limit;
                        double qty = limit + nlimit;
                        qty = qty <= 0 ? 1 : qty;
                        CapacityConstraint<Material> cc(qty);
                        cyclus::Material::Ptr manifest;
                        manifest = cyclus::ResCast<Material>(inventory[0].Pop());
                        Material::Ptr offer = Material::CreateUntracked(limit, manifest->comp());
                        inventory[0].Push(manifest);
                        manifest = cyclus::ResCast<Material>(inventory[1].Pop());
                        offer->Absorb(Material::CreateUntracked(nlimit, manifest->comp()));
                        inventory[1].Push(manifest);
                        port->AddBid(req, offer, this);
                        port->AddConstraint(cc);
                        ports.insert(port);
                        }
                    }
                }
            }
        }
    }

    void FuelfabFacility::AcceptMatlTrades(const std::vector< std::pair<cyclus::Trade<cyclus::Material>, cyclus::Material::Ptr> >& responses) {

        std::vector<std::pair<cyclus::Trade<cyclus::Material>, cyclus::Material::Ptr> >::const_iterator it;
        for (it = responses.begin(); it != responses.end(); ++it) {
            for(int i = 0; i < in_commods.size(); i++){
                if(it->first.request->commodity() == in_commods[i]){
                    inventory[i].Push(it->second);
                }
            }
        }
    }

    void FuelfabFacility::GetMatlTrades(const std::vector< cyclus::Trade<cyclus::Material> >& trades,
                                        std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                                        cyclus::Material::Ptr> >& responses) {
        using cyclus::Material;
        using cyclus::Trade;
        std::map<cyclus::Trader*, int> facility_request;
        std::vector< cyclus::Trade<cyclus::Material> >::const_iterator it;
        cyclus::Material::Ptr manifest;
        //Setting up comp map of requesters to number of requests
        for(it = trades.begin(); it != trades.end(); ++it){
            cyclus::Request<Material> req = *it->request;
            facility_request[req.requester()] += 1;
        }
        std::map<cyclus::Trader*, int>::iterator id;
        //iterate through requests
        for(id = facility_request.begin(); id != facility_request.end(); ++id){
            double limit, nlimit;
            //cast requester as reactor
            reactor::ReactorFacility* reactor = dynamic_cast<reactor::ReactorFacility*>(id->first);
            if (!reactor){
               throw cyclus::CastError("No reactor for fuelfab facility.");
            } else {
                //start up
                if (reactor->inventory.count() == 0){
                    limit = reactor->start_up(inventory);
                    double temp_limit = limit;
                    int k = 1;
                    for(it = trades.begin(); it!=trades.end(); ++it){
                        cyclus::Request<Material> req = *it->request;
                        if(req.requester() == id->first){
                            limit = temp_limit*k/(reactor->batches);
                            nlimit = reactor->core_mass/reactor->batches-limit;
                            manifest = cyclus::ResCast<Material>(inventory[0].Pop());
                            Material::Ptr offer = manifest->ExtractComp(0., manifest->comp());
                            //catch if over extracting from current fuel fab material
                            if(limit > manifest->quantity()){
                                double bonus = manifest->quantity();
                                offer->Absorb(manifest->ExtractComp(manifest->quantity(), manifest->comp()));
                                manifest = cyclus::ResCast<Material>(inventory[0].Pop());
                                offer->Absorb(manifest->ExtractComp(limit-bonus, manifest->comp()));
                            } else {
                                offer = manifest->ExtractComp(limit, manifest->comp());
                            }
                            std::cout << "limit " << limit << " : quant" << manifest->quantity() << std::endl;
                            inventory[0].Push(manifest);
                            manifest = cyclus::ResCast<Material>(inventory[1].Pop());
                            //catch for second resource buffer
                            if(nlimit > manifest->quantity()){
                                double bonus = manifest->quantity();
                                std::cout << "TEST" << std::endl;
                                offer->Absorb(manifest->ExtractComp(bonus, manifest->comp()));
                                manifest = cyclus::ResCast<Material>(inventory[1].Pop());
                                offer->Absorb(manifest->ExtractComp(nlimit-bonus, manifest->comp()));
                            } else {
                                std::cout << "ELSE TEST" << std::endl;
                                offer->Absorb(manifest->ExtractComp(nlimit, manifest->comp()));
                            }
                            std::cout << "nlimit " << nlimit << " : quant" << manifest->quantity() << std::endl;
                            inventory[1].Push(manifest);
                            responses.push_back(std::make_pair(*it, offer));
                            k++;
                        }
                    }
                } else{
                    //non-start up refueling
                    for(it = trades.begin(); it!=trades.end(); ++it){
                        cyclus::Request<Material> req = *it->request;
                        if(req.requester() == id->first){
                            limit = reactor->blend_next(inventory);
                            nlimit = reactor->core_mass/reactor->batches-limit;
                            manifest = cyclus::ResCast<Material>(inventory[0].Pop());
                            Material::Ptr offer = manifest->ExtractComp(0., manifest->comp());
                            if(limit > manifest->quantity()){
                                double bonus = manifest->quantity();
                                offer->Absorb(manifest->ExtractComp(manifest->quantity(), manifest->comp()));
                                manifest = cyclus::ResCast<Material>(inventory[0].Pop());
                                offer->Absorb(manifest->ExtractComp(limit-bonus, manifest->comp()));
                            } else {
                                offer = manifest->ExtractComp(limit, manifest->comp());
                            }
                            std::cout << "limit " << limit << " : quant" << manifest->quantity() << std::endl;
                            inventory[0].Push(manifest);
                            manifest = cyclus::ResCast<Material>(inventory[1].Pop());
                            if(nlimit > manifest->quantity()){
                                double bonus = manifest->quantity();
                                std::cout << "TEST_Base" << std::endl;
                                offer->Absorb(manifest->ExtractComp(bonus, manifest->comp()));
                                manifest = cyclus::ResCast<Material>(inventory[1].Pop());
                                offer->Absorb(manifest->ExtractComp(nlimit-bonus, manifest->comp()));
                            } else {
                                std::cout << "ELSE TEST_base" << std::endl;
                                offer->Absorb(manifest->ExtractComp(nlimit, manifest->comp()));
                            }
                            std::cout << "nlimit " << nlimit << " : quant" << manifest->quantity() << std::endl;
                            inventory[1].Push(manifest);
                            responses.push_back(std::make_pair(*it, offer));
                        }
                    }
                }
            }
        }
    }
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    extern "C" cyclus::Agent* ConstructFuelfabFacility(cyclus::Context* ctx) {
        return new FuelfabFacility(ctx);
    }

}
