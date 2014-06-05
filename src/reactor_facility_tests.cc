#include <gtest/gtest.h>

#include "reactor_facility.h"

#include "context.h"
#include "facility_tests.h"
#include "agent_tests.h"

using reactor::ReactorFacility;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
class ReactorFacilityTest : public ::testing::Test {
 protected:
  cyclus::TestContext tc_;
  ReactorFacility* src_facility_;

  virtual void SetUp() {
    src_facility_ = new ReactorFacility(tc_.get());
  }

  virtual void TearDown() {}
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(ReactorFacilityTest, clone) {
  ReactorFacility* cloned_fac =
      dynamic_cast<ReactorFacility*> (src_facility_->Clone());
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(ReactorFacilityTest, InitialState) {
  // Test things about the initial state of the facility here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(ReactorFacilityTest, Print) {
  EXPECT_NO_THROW(std::string s = src_facility_->str());
  // Test ReactorFacility specific aspects of the print method here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(ReactorFacilityTest, Tick) {
  ASSERT_NO_THROW(src_facility_->Tick());
  // Test ReactorFacility specific behaviors of the Tick function here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(ReactorFacilityTest, Tock) {
  EXPECT_NO_THROW(src_facility_->Tock());
  // Test ReactorFacility specific behaviors of the Tock function here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cyclus::Agent* ReactorFacilityConstructor(cyclus::Context* ctx) {
  return new ReactorFacility(ctx);
}

// required to get functionality in cyclus agent unit tests library
#ifndef CYCLUS_AGENT_TESTS_CONNECTED
int ConnectAgentTests();
static int cyclus_agent_tests_connected = ConnectAgentTests();
#define CYCLUS_AGENT_TESTS_CONNECTED cyclus_agent_tests_connected
#endif // CYCLUS_AGENT_TESTS_CONNECTED

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INSTANTIATE_TEST_CASE_P(ReactorFac, FacilityTests,
                        ::testing::Values(&ReactorFacilityConstructor));

INSTANTIATE_TEST_CASE_P(ReactorFac, AgentTests,
                        ::testing::Values(&ReactorFacilityConstructor));
