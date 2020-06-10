#ifndef PANDORA_KMERCOVERAGEMODELS_H
#define PANDORA_KMERCOVERAGEMODELS_H

#include <memory>
#include <random>
#include <boost/random/negative_binomial_distribution.hpp>
#include <boost/random/binomial_distribution.hpp>



enum ModelType { NegativeBinomialModelType, BinomialModelType, ConstantModelType };
/**
 * Represents a Random Number Generator generic model, which can produce random values.
 */
class RNGModel {
public:
    virtual uint32_t get_random_value() = 0;

    // NB: see function estimate_parameters() for details.
    static std::shared_ptr<RNGModel> create_new_RNG_model(ModelType model_type,
                                                          uint32_t exp_depth_covg,
                                                          double mean_kmer_coverage,
                                                          double binomial_parameter_p,
                                                          double negative_binomial_parameter_r,
                                                          double negative_binomial_parameter_p);

    bool operator==(const RNGModel& rhs) const
    {
        return random_number_generator == rhs.random_number_generator;
    }
    bool operator!=(const RNGModel& rhs) const { return !(rhs == *this); }

protected:
    std::default_random_engine random_number_generator;
    RNGModel() : random_number_generator(){}
};

/**
 * Represents a Negative Binomial Random Number Generator model, which can produce random values.
 */
class NegativeBinomialModel : public RNGModel {
public:
    NegativeBinomialModel(double negative_binomial_parameter_r, double negative_binomial_parameter_p) :
        RNGModel(),
        nb_distribution(negative_binomial_parameter_r, negative_binomial_parameter_p) {}

    virtual uint32_t get_random_value() {
        return uint32_t(std::max(nb_distribution(random_number_generator), 0));
    }

    bool operator==(const NegativeBinomialModel& rhs) const
    {
        return std::tie(static_cast<const RNGModel&>(*this), nb_distribution)
            == std::tie(static_cast<const RNGModel&>(rhs), rhs.nb_distribution);
    }
    bool operator!=(const NegativeBinomialModel& rhs) const { return !(rhs == *this); }

private:
    boost::random::negative_binomial_distribution<> nb_distribution;
};

/**
 * Represents a Binomial Random Number Generator model, which can produce random values.
 */
class BinomialModel : public RNGModel {
public:
    BinomialModel(double binomial_parameter_trials, double binomial_parameter_p) :
        RNGModel(),
        binomial_distribution(binomial_parameter_trials, binomial_parameter_p) {}

    virtual uint32_t get_random_value() {
        return uint32_t(std::max(binomial_distribution(random_number_generator), 0));
    }

    bool operator==(const BinomialModel& rhs) const
    {
        return std::tie(static_cast<const RNGModel&>(*this), binomial_distribution)
            == std::tie(static_cast<const RNGModel&>(rhs), rhs.binomial_distribution);
    }
    bool operator!=(const BinomialModel& rhs) const { return !(rhs == *this); }

private:
    boost::random::binomial_distribution<> binomial_distribution;
};

/**
 * Represents a Constant Model, which just produces always the same value. This is used to
 * cover the trivial case where no reads of a sample map, so it is a very specific and
 * not very useful case.
 */
class ConstantModel : public RNGModel {
public:
    ConstantModel(uint32_t value) :
        RNGModel(),
        value(value) {}

    virtual uint32_t get_random_value() {
        return value;
    }

    bool operator==(const ConstantModel& rhs) const
    {
        return std::tie(static_cast<const RNGModel&>(*this), value)
            == std::tie(static_cast<const RNGModel&>(rhs), rhs.value);
    }
    bool operator!=(const ConstantModel& rhs) const { return !(rhs == *this); }

private:
    uint32_t value;
};


/**
 * Container that stores and manages a vector of RNGModels, one per sample.
 * Holds a pointer to RNGModel so that we can apply polymorphism.
 */
class RNGModels : public std::vector<std::shared_ptr<RNGModel>>{
public:
using std::vector<std::shared_ptr<RNGModel>>::vector;
};

#endif // PANDORA_KMERCOVERAGEMODELS_H
