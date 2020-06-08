#ifndef __ESTIMATEPARAMETERS_H_INCLUDED__ // if estimate_parameters.h hasn't been
                                          // included yet...
#define __ESTIMATEPARAMETERS_H_INCLUDED__

#include <cstring>
#include <cstdint>
#include <boost/filesystem.hpp>
#include "GCPWrapper.h"

namespace fs = boost::filesystem;

double fit_mean_covg(const std::vector<uint32_t>&, const uint8_t);

double fit_variance_covg(const std::vector<uint32_t>&, double&, const uint8_t);

void fit_negative_binomial(double&, double&, float&, float&);

uint32_t find_mean_covg(std::vector<uint32_t>&);

int find_prob_thresh(std::vector<uint32_t>&);

/**
 * Estimate Kmer coverage models parameters for a sample.
 * This has the side effect of changing kmer_prg_with_coverage objects.
 * @return : the expected depth coverage of the sample, and a RNGModel object that models random coverages of this sample.
 */
using ExpDepthCovg = uint32_t;
std::pair<ExpDepthCovg, std::shared_ptr<RNGModel>> estimate_parameters(const std::shared_ptr<pangenome::Graph> &pangraph,
                                                                       const std::string& outdir, const uint32_t k, float& e_rate, const uint32_t covg,
                                                                       bool& bin, const uint32_t sample_id);

#endif
