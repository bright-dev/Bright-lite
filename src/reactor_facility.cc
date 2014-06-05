#include "reactor_facility.h"

namespace reactor {

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ReactorFacility::ReactorFacility(cyclus::Context* ctx)
    : cyclus::Facility(ctx) {};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
std::string ReactorFacility::str() {
  return Facility::str();
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void ReactorFacility::Tick() {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void ReactorFacility::Tock() {}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
extern "C" cyclus::Agent* ConstructReactorFacility(cyclus::Context* ctx) {
  return new ReactorFacility(ctx);
}

}
