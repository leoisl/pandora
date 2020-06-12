#include "GCPWrapper.h"

SampleInfo GCPSampleInfoModel::produce_data()
{
    // we always simulate 2-allele SNPs only
    SampleInfo sample_info(sample_index, 2, genotyping_options);

    uint32_t coverage_for_correct_allele = kmer_coverage_model->get_random_coverage_for_correct_allele();
    uint32_t coverage_for_incorrect_allele = kmer_coverage_model->get_random_coverage_for_incorrect_allele();

    // split the simulated coverages into forward and reverse coverages
    std::vector<std::vector<uint32_t>> allele_to_forward_coverages {
        { uint32_t(floor(coverage_for_correct_allele/2.0)) }, { uint32_t(floor(coverage_for_incorrect_allele/2.0)) }
    };
    std::vector<std::vector<uint32_t>> allele_to_reverse_coverages {
        { uint32_t(ceil(coverage_for_correct_allele/2.0)) }, { uint32_t(ceil(coverage_for_incorrect_allele/2.0)) }
    };

    sample_info.set_coverage_information(
        allele_to_forward_coverages, allele_to_reverse_coverages);

    return sample_info;
}


// NB : trivial method, not tested
GCPSampleInfoModels::GCPSampleInfoModels(const KmerCoverageModels &kmer_coverage_models, GenotypingOptions const* genotyping_options) {
    for (uint32_t sample_index = 0; sample_index < kmer_coverage_models.size(); ++sample_index) {
        std::shared_ptr<GCPSampleInfoModel> gcp_sample_info_model =
            std::make_shared<GCPSampleInfoModel>(sample_index, genotyping_options, kmer_coverage_models[sample_index]);
        push_back(gcp_sample_info_model);
    }
}


double GCPGenotyperAdapter::get_genotype_confidence() const
{
    double confidence { 0.0 };

    auto index_and_confidence_and_max_likelihood_optional
        = sample_info.get_confidence();
    if (index_and_confidence_and_max_likelihood_optional) {
        confidence = std::get<1>(*index_and_confidence_and_max_likelihood_optional);
    }

    return confidence;
}


// NB : trivial method, not tested
void GCPWrapper::train() {
    // simulate 10000000 confidences according to our Model and our Genotyper
    GCP::Simulator<SampleInfo, GCPGenotyperAdapter> genotype_confidence_simulator(model.get());
    std::vector<GCP::GenotypeConfidence> simulated_confidences = genotype_confidence_simulator.simulate(10000000);

    // create the Percentiler
    genotype_confidence_percentiler = std::make_shared<GCP::Percentiler>(simulated_confidences);
}


// NB : trivial method, not tested
GCPWrappers::GCPWrappers (const GCPSampleInfoModels &sample_info_models) {
    for (const std::shared_ptr<GCPSampleInfoModel> &sample_info_model : sample_info_models) {
        GCPWrapper gcp_wrapper(sample_info_model);
        push_back(gcp_wrapper);
    }
}