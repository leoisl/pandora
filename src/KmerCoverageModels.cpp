#include "KmerCoverageModels.h"


std::shared_ptr<RNGModel> RNGModel::create_new_RNG_model(ModelType model_type,
                                                         uint32_t exp_depth_covg,
                                                         double mean_kmer_coverage,
                                                         double binomial_parameter_p,
                                                         double negative_binomial_parameter_r,
                                                         double negative_binomial_parameter_p) {
    switch (model_type) {
    case ModelType::BinomialModelType:
        return std::make_shared<BinomialModel>(mean_kmer_coverage, binomial_parameter_p);
    case ModelType::NegativeBinomialModelType:
        return std::make_shared<NegativeBinomialModel>(negative_binomial_parameter_r, negative_binomial_parameter_p);
    case ModelType::ConstantModelType:
        return std::make_shared<ConstantModel>(exp_depth_covg);
    }
}