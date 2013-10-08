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

        cyclus::QueryEngine* output = qe->QueryElement("output");
        out_commods_.push_back(output->GetElementContent("outcommodity"));

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
        b->out_commods_ = OutputCommodities();

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
            cyclus::GenericResource::Create(Model::context(),
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

        market = cyclus::MarketModel::MarketForCommod(out_commods_[0]);
        recipient = dynamic_cast<cyclus::Communicator*>(market);

        trans = BuildTransaction();
        cyclus::GenericResource::Ptr trade_res =
            cyclus::GenericResource::Create(Model::context(),
                                            trade_res,
                                            "kg",
                                            out_commods_[0]);
        cyclus::Message::Ptr offer(new cyclus::Message(this, recipient, trans));
        offer.SendOn();
    }

    /**
    FOR TOCKIN'
    */
    void blmod::HandleTock(int time){
        std::cout << FacName() << "is tockin'\n";
    }

    /**
    FOR TRANSACTING!
    */
    cyclus::Transaction blmod::BuildTransaction() {
        using cyclus::Model;
        using cyclus::Material;
        std::map<int, double> isomap = std::map<int, double>();
        cyclus::Composition::Ptr out = cyclus::Composition::CreateFromMass(isomap);


        cyclus::Context* ctx = model::context();
        Material::Ptr trade_res = Material::Create(ctx, offer_amt, out);

        cyclus::Transaction trans(this, cyclus::OFFER);
        trans.SetCommod(out_commod_);
        trans.SetMinFrac(0.0);
        trans.SetPrice(commod_price_);
        trans.SetResource(trade_res);

        return trans;
    }
}
