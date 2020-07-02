#include "gtest/gtest.h"
#include "RNGModel.h"

TEST(NegativeBinomialModel, Fit_r_between_0_and_1_rounded_to_1)
{
    NegativeBinomialModel actual_model(0.25, 0.01);
    NegativeBinomialModel expected_model(1.0, 0.01);
    EXPECT_EQ(actual_model, expected_model);
}


TEST(NegativeBinomialModel, Fit_r_rounded_down)
{
    NegativeBinomialModel actual_model(2.49, 0.01);
    NegativeBinomialModel expected_model(2.0, 0.01);
    EXPECT_EQ(actual_model, expected_model);
}


TEST(NegativeBinomialModel, Fit_r_rounded_up)
{
    NegativeBinomialModel actual_model(3.51, 0.01);
    NegativeBinomialModel expected_model(4.0, 0.01);
    EXPECT_EQ(actual_model, expected_model);
}
