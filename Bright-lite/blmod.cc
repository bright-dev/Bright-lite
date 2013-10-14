#include "blmod.h"


blmod::blmod(cyclus::Context* ctx) : cyclus::FacilityModel(ctx),
    commod_price_(0){};

blmod::~blmod() {}

/**
FOR READING!
*/
void blmod::InitModuleMembers(cyclus::QueryEngine* qe){
    cyclus::QueryEngine* input = qe->QueryElement("input");

    cyclus::QueryEngine* commodities = input->QueryElement("commodities");
    in_commods_.push_back(commodities->GetElementContent("incommodity"));

    cyclus::QueryEngine* output = qe->QueryElement("output");
    out_commods_.push_back(output->GetElementContent("outcommodity"));

    batches = atoi(input->GetElementContent("batches").c_str());
    burnup = atof(input->GetElementContent("burnup").c_str());

    isoInformation iso_info;
    isoInformation iso_info1;
    std::ifstream inf("/home/robert/Bright-lite/u235data.txt");
    std::ifstream inf1("/home/robert/Bright-lite/u238data.txt");
    std::ofstream outf("/home/robert/Bright-lite/test.txt");
    if (!inf){
        cerr << "Could not read file yo for U-235\n";
    }
    if (!inf1){
        cerr << "Could not readdf file for U-238\n";
    }
    iso_info = BuildIsotope(inf);
    mass_stream.push_back(iso_info);
    iso_info1 = BuildIsotope(inf1);
    mass_stream.push_back(iso_info1);
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

    (*b).batches = batches;
    (*b).burnup = burnup;
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
    request->SendOn();

    market = cyclus::MarketModel::MarketForCommod(out_commod_);
    recipient = dynamic_cast<cyclus::Communicator*>(market);

    trans = BuildTransaction();
    cyclus::GenericResource::Ptr trade_res =
        cyclus::GenericResource::Create(Model::context(),
                                        1.0,
                                        "kg",
                                        out_commod_);
    cyclus::Message::Ptr offer(new cyclus::Message(this, recipient, trans));
    offer->SendOn();
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
    Material::Ptr mat = inventory_.PopOne();

    map<int, double> mass_frac = mat->comp()->mass();
    enrichment = mass_frac[92235]/mat->quantity();
    std::map<int, double> isomap = burnupcalc( FuelBuilder(mass_stream, enrichment), batches, 0.00001).second;
    cyclus::Composition::Ptr out = cyclus::Composition::CreateFromMass(isomap);

    cyclus::Context* ctx = Model::context();
    Material::Ptr trade_res = Material::Create(ctx, 1.0, out);

    cyclus::Transaction trans(this, cyclus::OFFER);
    trans.SetCommod(out_commod_);
    trans.SetMinFrac(0.0);
    trans.SetPrice(commod_price_);
    trans.SetResource(trade_res);

    return trans;
}
/**
FOR RESOURCES!
*/
void blmod::AddResource(cyclus::Transaction trans, std::vector<cyclus::Resource::Ptr> manifest) {
    inventory_.PushAll(cyclus::MatBuff::ToMat(manifest));
}





















