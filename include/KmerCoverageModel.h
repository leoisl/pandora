#ifndef PANDORA_KMERCOVERAGEMODEL_H
#define PANDORA_KMERCOVERAGEMODEL_H

#include "RNGModel.h"

class NotImplementedModel : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

/**
 * Represents a KmerCoverage Model, capable of generating kmer coverages of correct and incorrect alleles.
 * Kmer coverages are generated using RNGModels, and this model is parameterized according to estimate_parameters() function.
 */
class KmerCoverageModel {
public:
    KmerCoverageModel(ModelType model_type_for_correct_alleles,
                      ModelType model_type_for_incorrect_alleles,
                      uint32_t exp_depth_covg,
                      double error_rate,
                      double mean_kmer_coverage,
                      double binomial_parameter_p,
                      double negative_binomial_parameter_r,
                      double negative_binomial_parameter_p);

    uint32_t get_random_coverage_for_correct_allele() {
        return model_for_correct_alleles->get_random_value();
    }

    uint32_t get_random_coverage_for_incorrect_allele() {
        return model_for_incorrect_alleles->get_random_value();
    }

    bool operator==(const KmerCoverageModel& rhs) const
    {
        return std::tie(*model_for_correct_alleles, *model_for_incorrect_alleles)
            == std::tie(*rhs.model_for_correct_alleles, *rhs.model_for_incorrect_alleles);
    }
    bool operator!=(const KmerCoverageModel& rhs) const { return !(rhs == *this); }

private:
    std::shared_ptr<RNGModel> model_for_correct_alleles;
    std::shared_ptr<RNGModel> model_for_incorrect_alleles;
};


/**
 * Container that stores and manages a vector of KmerCoverageModel, one per sample.
 */
class KmerCoverageModels : public std::vector<std::shared_ptr<KmerCoverageModel>>{
public:
    using std::vector<std::shared_ptr<KmerCoverageModel>>::vector;
};

#endif // PANDORA_KMERCOVERAGEMODEL_H
