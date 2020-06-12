#include "KmerCoverageModel.h"

KmerCoverageModel::KmerCoverageModel(   ModelType model_type_for_correct_alleles,
                                        ModelType model_type_for_incorrect_alleles,
                                        uint32_t exp_depth_covg,
                                        double error_rate,
                                        double mean_kmer_coverage,
                                        double binomial_parameter_p,
                                        double negative_binomial_parameter_r,
                                        double negative_binomial_parameter_p) {
    switch (model_type_for_correct_alleles) {
    case ModelType::BinomialModelType:
        model_for_correct_alleles = std::make_shared<BinomialModel>(mean_kmer_coverage, binomial_parameter_p);
        break;
    case ModelType::NegativeBinomialModelType:
        model_for_correct_alleles = std::make_shared<NegativeBinomialModel>(negative_binomial_parameter_r, negative_binomial_parameter_p);
        break;
    case ModelType::ConstantModelType:
        model_for_correct_alleles = std::make_shared<ConstantModel>(exp_depth_covg);
        break;
    default:
        throw NotImplementedModel("To model coverage of correct alleles, use the Binomial, "
                                  "Negative Binomial or Constant models.");
    }

    switch (model_type_for_incorrect_alleles) {
    case ModelType::BinomialModelType:
        model_for_incorrect_alleles = std::make_shared<BinomialModel>(mean_kmer_coverage, error_rate);
        break;
    case ModelType::ConstantModelType:
        model_for_incorrect_alleles = std::make_shared<ConstantModel>((uint32_t)(exp_depth_covg * error_rate));
        break;
    // TODO: add NegativeBinomialModel here?
    default:
        throw NotImplementedModel("To model coverage of incorrect alleles, use the Binomial, "
                                  "or Constant models.");
    }
}