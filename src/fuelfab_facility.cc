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
        }
        for(int i = 0; i < inventory.size(); i++){
            if(inventory[i].count() == 0){
                inventory[i].set_capacity(maximum_storage/inventory.size());
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
        Material::Ptr target = Material::CreateUntracked(1,
                              Composition::CreateFromAtom(cm));
        RequestPortfolio<Material>::Ptr port(new RequestPortfolio<Material>());
        for(int i = 0; i < inventory.size(); i++){
            double qty = inventory[i].space();
            for(std::map<std::string, double>::iterator it = in_commods.begin(); it!= in_commods.end(); it++){
                port->AddRequest(target, this, it->first);
            }
            CapacityConstraint<Material> cc(qty);
            port->AddConstraint(cc);
            ports.insert(port);
        }
        return ports;
    }

    // MatlBids //
    std::set<cyclus::BidPortfolio<cyclus::Material>::Ptr>FuelfabFacility::GetMatlBids(
        cyclus::CommodMap<cyclus::Material>::type& commod_requests) {
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
        /** Quick Hack */
        CapacityConstraint<Material> cc(1);
        if (inventory_test == 0){std::cout << "YAY5?" << std::endl; return ports;}
        for (int i = 0; i < inventory.size(); i++){
            std::cout << "YAY5.5?" << std::endl;
            std::vector<cyclus::Material::Ptr> manifest;
            manifest = cyclus::ResCast<Material>(inventory[i].PopN(inventory[i].count()));
            CapacityConstraint<Material> cc(1);
            BidPortfolio<Material>::Ptr port(new BidPortfolio<Material>());
            std::vector<Request<Material>*>& requests = commod_requests[out_commod];
            std::vector<Request<Material>*>::iterator it;
            for (it = requests.begin(); it != requests.end(); ++it) {
                Request<Material>* req = *it;
                ReactorFacility* reactor = dynamic_cast<ReactorFacility*>(req->requester());
                if (reactor == NULL){
                   throw cyclus::CastError("Nope!");
                } else {
                    if (req->commodity() == out_commod) {
                        Material::Ptr offer = Material::CreateUntracked(20, manifest[0]->comp());
                        if (reactor->inventory.count() == 0){
                            std::cout << "yay" << std::endl;
                        } else if(reactor->burnup_test(offer) == reactor->target_burnup){
                            port->AddBid(req, offer, this);
                            port->AddConstraint(cc);
                            ports.insert(port);
                        };
                    }
                }
            }
            inventory[i].PushAll(manifest);
        }
        return ports;
    }

    void FuelfabFacility::AcceptMatlTrades(const std::vector< std::pair<cyclus::Trade<cyclus::Material>,
                                            cyclus::Material::Ptr> >& responses) {
        std::cout << "YAY8?" << std::endl;
        std::vector<std::pair<cyclus::Trade<cyclus::Material>, cyclus::Material::Ptr> >::const_iterator it;
        for (it = responses.begin(); it != responses.end(); ++it) {
            int j = 0;
            for (std::map<std::string, double>::iterator c = in_commods.begin(); c!= in_commods.end(); c++){
                if(it->first.request->commodity() == c->first){
                    inventory[j].Push(it->second);
                }
                j++;
            }
        }
    }

    void FuelfabFacility::GetMatlTrades(const std::vector< cyclus::Trade<cyclus::Material> >& trades,
                                        std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                                        cyclus::Material::Ptr> >& responses) {
        using cyclus::Material;
        using cyclus::Trade;
        std::cout << "YAY7?" << std::endl;
        std::vector< cyclus::Trade<cyclus::Material> >::const_iterator it;
        for(int i = 0; i < inventory.size(); i++){
            cyclus::Material::Ptr discharge = cyclus::ResCast<Material>(inventory[i].Pop());
            for (it = trades.begin(); it != trades.end(); ++it) {
                responses.push_back(std::make_pair(*it, discharge));
            }
        }
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    extern "C" cyclus::Agent* ConstructFuelfabFacility(cyclus::Context* ctx) {
        return new FuelfabFacility(ctx);
    }

}
