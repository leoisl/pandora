#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "test_helpers.h"

using namespace std;
using ::testing::Return;


TEST(GCPSampleInfoModelTest, produce_data___random_value_is_even)
{
    GCPSampleInfoModel sample_info_model(3, &default_genotyping_options, std::make_shared<ConstantModel>(10));

    SampleInfo actual_sample_info = sample_info_model.produce_data();
    SampleInfo expected_sample_info = SampleInfo(3, 2, &default_genotyping_options);
    expected_sample_info.set_coverage_information({{5}, {5}}, {{5}, {5}});

    EXPECT_EQ(actual_sample_info, expected_sample_info);
}

TEST(GCPSampleInfoModelTest, produce_data___random_value_is_odd)
{
    GCPSampleInfoModel sample_info_model(3, &default_genotyping_options, std::make_shared<ConstantModel>(11));

    SampleInfo actual_sample_info = sample_info_model.produce_data();
    SampleInfo expected_sample_info = SampleInfo(3, 2, &default_genotyping_options);
    expected_sample_info.set_coverage_information({{5}, {5}}, {{6}, {6}});

    EXPECT_EQ(actual_sample_info, expected_sample_info);
}


class GCPGenotyperAdapterTest___Fixture : public ::testing::Test {
public:
    class SampleInfoMock : public SampleInfo {
    public:
        using SampleInfo::SampleInfo;
        MOCK_METHOD(boost::optional<IndexAndConfidenceAndMaxLikelihood>, get_confidence,
                    (), (const override));
    };

    GCPGenotyperAdapterTest___Fixture()
        : sample_info_mock(0, 2, &default_genotyping_options),
          gcp_genotyper_adapter(sample_info_mock) {}

    void SetUp() override {}

    void TearDown() override {}

    SampleInfoMock sample_info_mock;
    GCPGenotyperAdapter gcp_genotyper_adapter;
};

TEST_F(GCPGenotyperAdapterTest___Fixture, get_genotype_confidence___invalid_confidence)
{
    EXPECT_CALL(sample_info_mock, get_confidence).Times(1).WillOnce(Return(boost::none));
    EXPECT_EQ(0, gcp_genotyper_adapter.get_genotype_confidence());
}


TEST_F(GCPGenotyperAdapterTest___Fixture, get_genotype_confidence___valid_confidence)
{
    EXPECT_CALL(sample_info_mock, get_confidence).Times(1).WillOnce(Return(std::make_tuple((size_t)3, 8.5, -10.5)));
    EXPECT_NEAR(8.5, gcp_genotyper_adapter.get_genotype_confidence(), 0.000001);
}
