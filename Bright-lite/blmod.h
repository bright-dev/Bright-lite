#ifndef BLMOD_H_INCLUDED
#define BLMOD_H_INCLUDED

#include <string>
#include <stdlib.h>

#include "facility_model.h"
#include "query_engine.h"
#include "mat_buff.h"
#include "logger.h"

#include "structures.h"
#include "parser.h"

#include <sstream>
#include <limits>

#include <boost/lexical_cast.hpp>

#include "origenBuilder.h"
#include "burnupcalc.h"

#include "context.h"
#include "logger.h"
#include "generic_resource.h"
#include "error.h"
#include "cyc_limits.h"
#include "market_model.h"


class blmod : public cyclus::FacilityModel {
    public:
    blmod(cyclus::Context* ctx);
    virtual ~blmod();

    /**
    */
    vector<isoInformation> mass_stream;

    virtual cyclus::Model* Clone();
    /**

    */
    virtual void InitModuleMembers(cyclus::QueryEngine* qe);

    /**
    */
    std::string schema();

    /**
    */
    virtual std::string str();

    /* Agent Methods */
    /**
    */
    virtual void HandleTick(int time);

    /**
    */
    virtual void HandleTock(int time);

    /**
    */
    virtual void AddResource(cyclus::Transaction trans,
                           std::vector<cyclus::Resource::Ptr> manifest);
    /**
    */
    void SendOffer(cyclus::Transaction trans);
    /**
    */
    virtual void ReceiveMessage(cyclus::Message::Ptr msg){};

    /**
    */
    protected:
    cyclus::Transaction BuildTransaction();

    /* blmod methods! */


    /* blmod attributes */
    int batches;
    double burnup;
    double enrichment;
    double commod_price_;
    cyclus::MatBuff inventory_;
    std::vector<std::string> in_commods_;
    std::string out_commod_;
};

#endif // BLMOD_H_INCLUDED
