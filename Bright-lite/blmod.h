#ifndef BLMOD_H_INCLUDED
#define BLMOD_H_INCLUDED

#include <string>
#include <stdlib.h>

#include "facility_model.h"
#include "query_engine.h"
#include "mat_buf.h"
#include "logger.h"

#include "structures.h"
#include "parser.h"

#include <sstream>
#include <limits>

#include <boost/lexical_cast.hpp>

#include "sink_facility.h"

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

    virtual cyclus::Model* Clone();
    /**

    */
    virtual void InitModuleMembers(cyclus::QueryEngine* qe);

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
    void SendOffer(cyclus::Transaction trans);
    /**
    */
    virtual void RecieveMessage(cyclus::Message::Ptr msg) {};

    /**
    */
    protected:
    cyclus::Transaction BuildTransaction();

    /* blmod methods! */


    /* blmod attributes */
    int batches;
    double burnup;
    cyclus::MatBuff inventory_;
    std::string out_commod_;
}

#endif // BLMOD_H_INCLUDED
