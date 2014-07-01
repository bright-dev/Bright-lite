#include <gtest/gtest.h>

#include "fuelfab_facility.h"

#include "context.h"
#include "facility_tests.h"
#include "agent_tests.h"

using fuelfab::FuelfabFacility;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
class FuelfabFacilityTest : public ::testing::Test {
 protected:
  cyclus::TestContext tc_;
  FuelfabFacility* src_facility_;

  virtual void SetUp() {
    src_facility_ = new FuelfabFacility(tc_.get());
  }

  virtual void TearDown() {}
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(FuelfabFacilityTest, clone) {
  FuelfabFacility* cloned_fac =
      dynamic_cast<FuelfabFacility*> (src_facility_->Clone());
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(FuelfabFacilityTest, InitialState) {
  // Test things about the initial state of the facility here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(FuelfabFacilityTest, Print) {
  EXPECT_NO_THROW(std::string s = src_facility_->str());
  // Test FuelfabFacility specific aspects of the print method here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(FuelfabFacilityTest, Tick) {
  ASSERT_NO_THROW(src_facility_->Tick());
  // Test FuelfabFacility specific behaviors of the Tick function here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(FuelfabFacilityTest, Tock) {
  EXPECT_NO_THROW(src_facility_->Tock());
  // Test FuelfabFacility specific behaviors of the Tock function here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cyclus::Agent* FuelfabFacilityConstructor(cyclus::Context* ctx) {
  return new FuelfabFacility(ctx);
}

// required to get functionality in cyclus agent unit tests library
#ifndef CYCLUS_AGENT_TESTS_CONNECTED
int ConnectAgentTests();
static int cyclus_agent_tests_connected = ConnectAgentTests();
#define CYCLUS_AGENT_TESTS_CONNECTED cyclus_agent_tests_connected
#endif // CYCLUS_AGENT_TESTS_CONNECTED

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
INSTANTIATE_TEST_CASE_P(FuelfabFac, FacilityTests,
                        ::testing::Values(&FuelfabFacilityConstructor));

INSTANTIATE_TEST_CASE_P(FuelfabFac, AgentTests,
                        ::testing::Values(&FuelfabFacilityConstructor));
