#include "blmod.h"


using namespace cycamore {
    blmod::blmod(cyclus::Context* ctx) : cyclus::FacilityModel(ctx),
    commod_price_(0),
    capcity(std::numeric_limits<double>::max()) {}

    blmod::~blmod() {}

    /**
    FOR READING!
    */
    void InitModuleMembers(cyclus::QueryEngine* qe){
        cyclus::QueryEngine* input = qe->QueryElement("input");

        cyclus::QueryEngine* commodities = input->QueryElement("commodities");
        in_commods_.push_back(commodities->GetElementContent("incommodity"));

        batches = atoi(input->GetElementContent("batches").c_str());
        burnup = atof(input->GetElementContent("burnup").c_str());
    }

    /**
    FOR NAMES!
    */
    std::string blmod::str() {
        return "I AM BLMOD!";
    }

    /**
    FOR CLONES!
    */
    cyclus::Model* blmod::Clone() {
        blmod* b = new blmod(*this);
        b->InitFrom(this);

        b.batches = batches;
        b.burnup = burnup;
        b->in_commods_ = InputCommodities();

        return b;
    }
    /**
    FOR TICKING!
    */
    void blmod::HandleTick(int time){
        using std::string;

        double requestAmt = 1.;

        cyclus::MarketModel* market = cyclus::MarketModel::MarketForCommod(in_commods_[0]);
        cyclus::Communicator* recipient = dynamic_cast<cyclus::Communicator*>(market);

        cyclus::GenericResource::Ptr request_res =
            cyclus:: GenericResource::Create(Model::context(),
                                             requestAmt,
                                             "kg",
                                             in_commods_[0]);

        cyclus::Transaction trans(this, cyclus::REQUEST);
        trans.SetCommod(in_commods_[0]);
        trans.SetMinFrac(0.);
        trans.SetPrice(commod_price_);
        trans.SetResource(request_res);

        cyclus::Message::Ptr request(new cyclus::Message(this, recipient, trans));
        request.SendOn();
    }

    /**
    FOR TOCKIN'
    */
    void blmod::HandleTock(int time){
        std::cout << FacName() << "is tockin'\n";
    }

    /**
    FOR RESOURCES!
    */

}
