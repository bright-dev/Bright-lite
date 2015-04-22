#include "reprocess_facility.h"



namespace reprocess {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ReprocessFacility::ReprocessFacility(cyclus::Context* ctx)
        : cyclus::Facility(ctx) {};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    std::string ReprocessFacility::str() {
      return Facility::str();
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void ReprocessFacility::Tick() {
    //inventory check is done in this phase, inventories initialized in the beginning of simulation
    	//std::cout << std::endl << std::endl << "~~tick~~" << std::endl;
    	//std::cout << "Input inventory: " << input_inventory.quantity() << std::endl;
    	//std::cout << "Waste inventory: " << waste_inventory.quantity() << std::endl;
    	/*for(int i = 0; i < out_inventory.size(); i++){
    	  std::cout << "Out inventory " << i << ": " << out_inventory[i].quantity() << std::endl;
    	}*/

    	//check to see beginning of simulation
    	if(pi != 3.141592){
    	  pi = 3.141592;
    	  int nucid = 0; //for erroneous inputs this could be a long int
    	  double eff;
    	  std::string line;
    	  //populate out_commods, out_eff from input, initiate each associated output_inventory
    	  std::ifstream fin(repro_input_path);
           if(!fin.good()){
             std::cout << "Error! Failed reading reprocessing plant input file." << '\n';
             std::cout << "  Given path: " << repro_input_path << '\n';
             std::cout << "  Please correct the path in cyclus input file for repro_input_path variable." << '\n';
           } else {
             std::cout << "Reprocess facility starting (ID:" << id() << ") reading input " << repro_input_path << std::endl;
           }
           //reads the input file and populates out_eff with the read values
          while(getline(fin, line)){
            if(line.find("BEGIN") == 0){
              cyclus::toolkit::ResourceBuff temp_pass;
              out_inventory.push_back(temp_pass); //initializes the corresponding inventory
              std::map<int, double> tempmap;
              while(getline(fin, line)){
                std::istringstream iss(line);
                if(line.find("END") == 0){break;}
                iss >> nucid >> eff;
                //std::cout << " read input: " << nucid << "  " << eff << std::endl;
                if(nucid < 10000000 || nucid > 2000000000){
                  std::cout << "Error reading isotope identifier in reprocess facility." << '\n';
                  std::cout << "Isotope: " << nucid << " isn't in the form 'zzaaammmm'." << '\n';
                }
                if(eff < 0 || eff > 1){
                  std::cout << "Error in isotope removal efficiency in reprocess facility." << '\n';
                  std::cout << "  Efficiency is set to " << eff << " for " << nucid << "." << '\n';
                }
                tempmap[nucid] = eff;
              }
            out_eff.push_back(tempmap);
            }

          }
         if(out_eff.size() == 0){
           std::cout << "Error populating removal efficiencies in reprocess facility. Efficiencies not read in cyclus." << '\n';
         }
    	}

    	//std::cout << "-//tick//-" << std::endl;
    }

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    void ReprocessFacility::Tock() {
    //the reprocessing is done in this phase
    	//std::cout << "~~REPOtock~~" << std::endl;

    	//Creates a vector Material (manifest) that has been popped out.
    	std::vector<cyclus::Material::Ptr> manifest;
    	 double extract = output_capacity >= input_inventory.quantity() ? input_inventory.quantity() : output_capacity;
         manifest = cyclus::ResCast<cyclus::Material>(input_inventory.PopQty(extract));
         double tot_mass = 0; //total mass of the reprocessed batch (all isotopes)
         //std::cout << "Amount of inventory batches being reprocessed: " << manifest.size() << '\n';

         for(int o = 0; o < out_eff.size(); o++){
             for(int m = 0; m < manifest.size(); m++){
               cyclus::CompMap temp_comp; //stores the masses of extracted isotopes
               //cyclus::CompMap out_comp = out_eff[o].mass();
               cyclus::CompMap mani_comp = manifest[m]->comp()->mass();
               double  size_tad = manifest[m]->quantity();
               //std::cout << "Size M "<< m << " " << size_tad << std::endl;
               std::map<int, double>::iterator out_it;
               cyclus::CompMap::iterator mani_it;
               for(out_it = out_eff[o].begin(); out_it != out_eff[o].end(); ++out_it){
                 for(mani_it = mani_comp.begin(); mani_it != mani_comp.end(); ++mani_it){
                   if(out_it->first == mani_it->first){
                     //std::cout << mani_it->first << "   " << mani_it->second << "   " << out_it->second << std::endl;
                     temp_comp[mani_it->first] = mani_it->second * out_it->second;
                     tot_mass += size_tad * mani_it->second * out_it->second;
                     //std::cout << "TOT+AMSS " << tot_mass <<std::endl;
                   }
                 }
               }
               /*for(mani_it = temp_comp.begin(); mani_it != temp_comp.end(); ++mani_it){
                    std::cout << mani_it->first << " amount " << mani_it->second << std::endl;
               }

               for(mani_it = mani_comp.begin(); mani_it != mani_comp.end(); ++mani_it){
                    std::cout << mani_it->first << " amount " << mani_it->second << std::endl;
               }*/
               //Puts the extracted material in the corresponding out_inventory
               //std::cout << "TOTAL MASS "<<tot_mass << std::endl;
               //std::cout << manifest[m]->quantity() << std::endl;
               cyclus::Resource::Ptr resource = cyclus::ResCast<cyclus::Resource>(manifest[m]->ExtractComp(tot_mass, cyclus::Composition::CreateFromMass(temp_comp)));
               //std::cout << tot_mass << std::endl;
               out_inventory[o].Push(resource);
               //std::cout << tot_mass << std::endl;
               tot_mass = 0;
             }
         }
         //adds the remaining materials (the waste) in the waste_inventory
         waste_inventory.PushAll(manifest);

    	//std::cout << "-//REPOtock//-" << std::endl;
    }

    std::set<cyclus::RequestPortfolio<cyclus::Material>::Ptr> ReprocessFacility::GetMatlRequests() {
        //std::cout << "~~REPOGetReq~~" << std::endl;
        using cyclus::RequestPortfolio;
        using cyclus::Material;
        using cyclus::Composition;
        using cyclus::CompMap;
        using cyclus::CapacityConstraint;
        std::set<RequestPortfolio<Material>::Ptr> ports;
        cyclus::Context* ctx = context();
        CompMap cm;

        Material::Ptr target = Material::CreateUntracked(input_capacity, Composition::CreateFromAtom(cm));

        RequestPortfolio<Material>::Ptr port(new RequestPortfolio<Material>());

        double qty = input_inventory.space();
        //std::cout << "REPO " << qty << std::endl;
        port->AddRequest(target, this, in_commod[0]);

        CapacityConstraint<Material> cc(qty);
        port->AddConstraint(cc);
        ports.insert(port);



        //std::cout << "-//REPOGetReqEND//-" << std::endl;
        return ports;
    }

    // MatlBids //
    std::set<cyclus::BidPortfolio<cyclus::Material>::Ptr>ReprocessFacility::GetMatlBids(
            cyclus::CommodMap<cyclus::Material>::type& commod_requests) {
        //std::cout << "~~REPOGetBid~~" << std::endl;

        using cyclus::Bid;
        using cyclus::BidPortfolio;
        using cyclus::CapacityConstraint;
        using cyclus::Material;
        using cyclus::Request;
        using cyclus::Converter;

        std::set<BidPortfolio<Material>::Ptr> ports;

        if(waste_inventory.quantity() != 0){
            std::vector<cyclus::Material::Ptr> manifest = cyclus::ResCast<Material>(waste_inventory.PopN(waste_inventory.count()));
            std::vector<Request<Material>*>& requests = commod_requests[commod_out[commod_out.size()-1]];
            std::vector<Request<Material>*>::iterator it;

            for (it = requests.begin(); it != requests.end(); ++it) {
                BidPortfolio<Material>::Ptr port(new BidPortfolio<Material>());
                Request<Material>* req = *it;
                Material::Ptr offer = Material::CreateUntracked(output_capacity, manifest[0]->comp());
                port->AddBid(req, offer, this);
                ports.insert(port);
            }
            waste_inventory.PushAll(manifest);
        }
        for(int i = 0; i < out_inventory.size(); i++){
            if(out_inventory[i].count() != 0){
                std::vector<cyclus::Material::Ptr> manifest = cyclus::ResCast<Material>(out_inventory[i].PopN(out_inventory[i].count()));
                std::vector<Request<Material>*>& requests = commod_requests[commod_out[i]];
                std::vector<Request<Material>*>::iterator it;

                for (it = requests.begin(); it != requests.end(); ++it) {
                    BidPortfolio<Material>::Ptr port(new BidPortfolio<Material>());
                    Request<Material>* req = *it;
                    Material::Ptr offer = Material::CreateUntracked(output_capacity, manifest[0]->comp());
                    port->AddBid(req, offer, this);
                    ports.insert(port);
                }
                out_inventory[i].PushAll(manifest);
            }
        }
        //std::cout << "~~REPOGetBidEND~~" << std::endl;
        return ports;

    }

void ReprocessFacility::AcceptMatlTrades(const std::vector< std::pair<cyclus::Trade<cyclus::Material>, cyclus::Material::Ptr> >& responses) {
    //std::cout << "~~REPOAcptTrades~~" << std::endl;

    std::vector< std::pair<cyclus::Trade<cyclus::Material>, cyclus::Material::Ptr> >::const_iterator it;
    for (it = responses.begin(); it != responses.end(); ++it) {
        //std::cout << it->first.request->commodity() << std::endl;
        for(int i = 0; i < in_commod.size(); i++){
            if(it->first.request->commodity() == in_commod[i]){
                input_inventory.Push(it->second);
            }
        }
    }

    //std::cout << "-//REPO ENDAcptTrades//-" << std::endl;
}

void ReprocessFacility::GetMatlTrades(const std::vector< cyclus::Trade<cyclus::Material> >& trades,
    std::vector<std::pair<cyclus::Trade<cyclus::Material>,cyclus::Material::Ptr> >& responses) {
    //std::cout << "REPOMAT TRADES BEGIN" << std::endl;
    using cyclus::Material;
    using cyclus::Trade;
    std::vector< cyclus::Trade<cyclus::Material> >::const_iterator it;

    std::vector<cyclus::Material::Ptr> waste = cyclus::ResCast<Material>(waste_inventory.PopN(waste_inventory.count()));
    for(int i = 1; i < waste.size(); i++){
        waste[0]->Absorb(waste[i]);
    }
    for (it = trades.begin(); it != trades.end(); ++it) {
        //std::cout << "Trade Commodity: "<<it->request->commodity() << std::endl;
        //std::cout << "Trade Commodity: "<<it->request->commodity() << std::endl;
        if(it->request->commodity() == commod_out[commod_out.size()-1]){
            responses.push_back(std::make_pair(*it, waste[0]));
        }
    }
    for(int i = 0; i < out_inventory.size(); i++){
        std::vector<cyclus::Material::Ptr> discharge = cyclus::ResCast<Material>(out_inventory[i].PopN(out_inventory[i].count()));
        for(int j = 1; j < discharge.size(); j++){
            discharge[0]->Absorb(discharge[j]);
        }
        for (it = trades.begin(); it != trades.end(); ++it){
            if(it->request->commodity() == commod_out[i]){
                responses.push_back(std::make_pair(*it, discharge[0]));
            }
        }
    //std::cout << "REPOMAT TRADES END" << std::endl;
    }

}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    extern "C" cyclus::Agent* ConstructReprocessFacility(cyclus::Context* ctx) {
        return new ReprocessFacility(ctx);
    }


}
