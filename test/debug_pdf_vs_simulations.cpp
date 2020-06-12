#include "gtest/gtest.h"
#include "RNGModel.h"
#include <boost/math/distributions/negative_binomial.hpp>
#include <fstream>
#include <iostream>
#include <chrono>

/**
 * Code used to debug if the simulated coverages model well Negative Binomial distribution.
 * See https://github.com/leoisl/pandora_debug/blob/master/debug_percentiles/debug_percentiles.ipynb
 * See https://github.com/leoisl/GCP/issues/11
 */


using namespace std;

vector<double> get_pdf(double negative_binomial_parameter_r,
    double negative_binomial_parameter_p, uint32_t limit=300) {
    boost::math::negative_binomial_distribution<> nb_distribution(negative_binomial_parameter_r, negative_binomial_parameter_p);

    std::vector<double> probs;
    for (uint32_t i=0; i<limit; ++i) {
        probs.push_back(boost::math::pdf(nb_distribution, i));
    }
    return probs;
}

vector<double> get_pdf_based_on_simulations(double negative_binomial_parameter_r,
                                            double negative_binomial_parameter_p,
                                            uint32_t nb_of_simulations,
                                            uint32_t limit=300) {
    NegativeBinomialModel negative_binomial_model(
        negative_binomial_parameter_r, negative_binomial_parameter_p);
    map<uint32_t , uint32_t> counter;
    for (uint32_t i=0; i<nb_of_simulations; ++i) {
        uint32_t value = negative_binomial_model.get_random_value();
        counter[value]++;
    }

    std::vector<double> probs;
    for (uint32_t i=0; i<limit; ++i) {
        probs.push_back(((double)counter[i])/nb_of_simulations);
    }
    return probs;
}

void print_vector_to_file(const std::vector<double> &pdf, const std::string &filename) {
    ofstream output(filename);
    output << "trial,prob" << endl;

    uint32_t trial = 0;
    for (double value : pdf) {
        output << trial << "," << value << endl;
        ++trial;
    }
    output.close();
}

void debug_pdf_vs_simulations() {
//    double negative_binomial_parameter_r = 2.41921; // the value for 063_STEC
    double negative_binomial_parameter_r = 2; // rounded value
    double negative_binomial_parameter_p = 0.0335631;

    std::vector<uint32_t> nb_of_simulations{1000, 10000, 100000, 1000000, 10000000};

    std::vector<double> pdf = get_pdf(negative_binomial_parameter_r, negative_binomial_parameter_p);
    print_vector_to_file(pdf, "pdf");

    for (uint32_t nb_of_simulation : nb_of_simulations) {
        std::vector<double> simulated_pdf = get_pdf_based_on_simulations(negative_binomial_parameter_r,
            negative_binomial_parameter_p, nb_of_simulation);
        std::string filename = std::string("simulated_pdf_") + to_string(nb_of_simulation);
        print_vector_to_file(simulated_pdf, filename);
    }
}


// Uncomment this to run debug_pdf_vs_simulations()
//TEST(debug_pdf_vs_simulations_dummy, debug_pdf_vs_simulations_dummydummy) {
//    debug_pdf_vs_simulations();
//}