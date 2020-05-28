#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>
#include <numeric>
#include <algorithm>
#include "utils.h"
#include "pangenome/pangraph.h"
#include "pangenome/pannode.h"
#include "estimate_parameters.h"

#define assert_msg(x) !(std::cerr << "Assertion failed: " << x << std::endl)

double fit_mean_covg(
    const std::vector<uint32_t>& kmer_covg_dist, const uint8_t zero_thresh)
{

    double sum = 0;
    double total = 0;
    for (uint32_t i = zero_thresh; i < kmer_covg_dist.size(); ++i) {
        sum += kmer_covg_dist[i] * i;
        total += kmer_covg_dist[i];
    }

    if (total == 0)
        return 0;

    BOOST_LOG_TRIVIAL(info) << "Found mean " << sum / total;
    return sum / total;
}

double fit_variance_covg(const std::vector<uint32_t>& kmer_covg_dist, double& mean,
    const uint8_t zero_thresh)
{
    double acc = 0;
    double total = 0;
    for (uint32_t i = zero_thresh; i < kmer_covg_dist.size(); ++i) {
        acc += (i - mean) * (i - mean) * kmer_covg_dist[i];
        total += kmer_covg_dist[i];
    }

    if (total == 0)
        return 0;

    BOOST_LOG_TRIVIAL(info) << "Found variance " << acc / total;
    return acc / total;
}

void fit_negative_binomial(double& mean, double& variance, float& p, float& r)
{
    assert(mean > 0 and variance > 0);
    assert(mean < variance);
    p = mean / variance;
    r = (mean * p / (1 - p) + variance * p * p / (1 - p)) / 2;
    BOOST_LOG_TRIVIAL(info) << "Negative binomial parameters p: " << p
                             << " and r: " << r;
}

uint32_t find_mean_covg(std::vector<uint32_t>& kmer_covg_dist)
{
    // tries to return the position in vector at which the maximum of the second peak
    // occurs expects at least 3 increases of covg to decide we are out of the first
    // peak Note that if there are localPRGs which occur twice (or more) we expect the
    // kmer covgs to have an additional smaller peak(s)
    bool first_peak = true;
    uint32_t max_covg = 0;
    uint32_t noise_buffer = 0;
    uint32_t run_of_zeroes = 0;

    for (uint32_t i = 1; i != kmer_covg_dist.size(); ++i) {
        if (kmer_covg_dist[i] == 0 and kmer_covg_dist[i - 1] == 0) {
            run_of_zeroes += 1;
        }
        if (kmer_covg_dist[i] <= kmer_covg_dist[i - 1]) {
            // only interested in where we stop being in a decreasing section
            continue;
        } else if (first_peak and run_of_zeroes < 10 and noise_buffer < 3) {
            // increase buffer -> have to see several increases to believe not noise
            noise_buffer += 1;
            run_of_zeroes = 0;
            continue;
        } else if (first_peak) {
            // have seen several increases now, so probably out of first peak
            first_peak = false;
            run_of_zeroes = 0;
            max_covg = i;
        } else if (kmer_covg_dist[i] > kmer_covg_dist[max_covg]) {
            // have a new max
            max_covg = i;
        }
    }

    // if still have first_peak true, was a mistake to try with this covg
    if (first_peak) {
        BOOST_LOG_TRIVIAL(info)
            << "Did not find 2 distinct peaks - use default error rate";
        max_covg = 0;
    }

    return max_covg;
}

int find_prob_thresh(std::vector<uint32_t>& kmer_prob_dist)
{
    // finds position at which minimum occurs between two peaks
    // naive way is to pick window with minimal positive covg
    // sligtly less naive way is to find this minimal value between two peaks more than
    // 10 apart Note we are working with log probs so expect all negative, lots near 0
    // and some larger negatives. We expect the threshold to be not very far from 0 Also
    // note that we have forced kmer_prob_dist to have size 200 so we can cast all the
    // long things as ints since we know they are between 0 and 200
    if (kmer_prob_dist.empty()) {
        return 0;
    }

    std::vector<uint32_t> kmer_prob_dist_temp = kmer_prob_dist; //makes a copy so that we can manipulate it

    // first find the first peak
    int first_peak = (int) distance(kmer_prob_dist_temp.begin(),
        max_element(kmer_prob_dist_temp.begin(), kmer_prob_dist_temp.end()));

    // then remove the first peak and values surrounding it, to not get any close value
    int surrounding_length = 3;
    for (int i=first_peak-surrounding_length; i<=first_peak+surrounding_length; ++i) {
        if (i>=0 and i<kmer_prob_dist_temp.size()) {
            kmer_prob_dist_temp[i] = 0;
        }
    }

    // find second peak
    int second_peak = (int) distance(kmer_prob_dist_temp.begin(),
                                    max_element(
                                        kmer_prob_dist_temp.begin(), kmer_prob_dist_temp.end()));

    // ensure first peak < second peak
    if (first_peak > second_peak) {
        // swap
        int temp_peak = first_peak;
        first_peak = second_peak;
        second_peak = temp_peak;
    }

    // peak offset to skip the artificial minimum
    int peak_offset = 2;

    int minimum_peak_between_the_two_highest_peaks =
        (int)distance(kmer_prob_dist.begin(),
            min_element(kmer_prob_dist.begin() + first_peak + peak_offset, kmer_prob_dist.begin() + second_peak));


    BOOST_LOG_TRIVIAL(info) << "Minimum found between " << (first_peak - 200) << " and "
                             << (second_peak - 200) << " at " << (minimum_peak_between_the_two_highest_peaks - 200);
    return minimum_peak_between_the_two_highest_peaks - 200;
}

uint32_t estimate_parameters(std::shared_ptr<pangenome::Graph> pangraph,
    const std::string& outdir, const uint32_t k, float& e_rate, const uint32_t covg,
    bool& bin, const uint32_t& sample_id)
{
    BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] Start";
    BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] k: " << k;
    BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] error rate: " << e_rate;
    BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] Global coverage: " << covg;
    BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] Binomial? " << bin;
    BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] sample_id: " << bin;


    uint32_t exp_depth_covg = covg;

    // ignore trivial case
    if (pangraph->nodes.empty()) {
        BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] Trivial case, expected coverage: " << exp_depth_covg;
        BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] End";
        return exp_depth_covg;
    }

    std::vector<uint32_t> kmer_covg_dist(1000,
        0); // empty vector of zeroes to represent kmer coverages between 0 and 1000
    std::vector<uint32_t> kmer_prob_dist(
        200, 0); // empty vector of zeroes to represent the
    // distribution of log probs between -200 and 0
    uint32_t c, mean_covg;
    unsigned long num_reads = 0;
    int thresh;
    float binomial_parameter_p;
    float negative_binomial_parameter_p = 0, negative_binomial_parameter_r = 0;

    // first we estimate error rate
    for (const auto& node : pangraph->nodes) {
        num_reads += node.second->covg;
        for (uint32_t i = 1;
             i != node.second->kmer_prg_with_coverage.kmer_prg->nodes.size() - 1;
             ++i) // NB first and last kmer in kmergraph are null
        {
            c = node.second->kmer_prg_with_coverage.get_covg(i, 0, sample_id)
                + node.second->kmer_prg_with_coverage.get_covg(i, 1, sample_id);
            if (c < 1000) {
                kmer_covg_dist[c] += 1;
            }
        }
    }

    // estimate actual coverage on these prgs
    num_reads = num_reads / pangraph->nodes.size();
    BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] Coverage on PRGs: " << num_reads;

    // save coverage distribution
    fs::create_directories(outdir);
    std::ofstream handle;
    handle.open(outdir + "/kmer_covgs.txt");
    assert(!handle.fail()
        or assert_msg("Could not open file " << outdir + "/kmer_covgs.txt"));
    for (uint32_t j = 0; j != kmer_covg_dist.size(); ++j) {
        handle << j << "\t" << kmer_covg_dist[j] << std::endl;
    }
    handle.close();

    // evaluate error rate
    auto mean = fit_mean_covg(kmer_covg_dist, covg / 10);
    BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] Kmer coverage mean: " << mean;

    auto var = fit_variance_covg(kmer_covg_dist, mean, covg / 10);
    BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] Kmer coverage variance: " << var;

    if (mean > var) {
        BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] Mean > var, reestimating...";

        auto zero_thresh = 2;
        mean = fit_mean_covg(kmer_covg_dist, zero_thresh);
        BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] Reestimated kmer coverage mean: " << mean;

        var = fit_variance_covg(kmer_covg_dist, mean, zero_thresh);
        BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] Reestimated kmer coverage variance: " << var;
    }

    if ((bin and num_reads > 30 and covg > 30)
        or (not bin and abs(var - mean) < 2 and mean > 10 and num_reads > 30
            and covg > 2)) {
        BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] Fitting binomial... ";
        bin = true;
        mean_covg = find_mean_covg(kmer_covg_dist);
        BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] Fitting binomial: Mean kmer coverage: " << mean_covg;

        if (exp_depth_covg < 1)
            exp_depth_covg = mean;
        BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] Fitting binomial: Expected coverage: " << exp_depth_covg;

        auto old_e_rate = e_rate;
        if (mean_covg > 0 and mean_covg < covg) {
            e_rate = -log((float)mean_covg / covg) / k;
        }
        BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] Fitting binomial: Estimated error rate updated from " << old_e_rate << " to " << e_rate;


    } else if (not bin and num_reads > 30 and covg > 2 and mean < var) {
        BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] Fitting negative binomial... ";
        fit_negative_binomial(
            mean, var, negative_binomial_parameter_p, negative_binomial_parameter_r);
        exp_depth_covg = mean;
        BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] Fitting negative binomial: Negative binomial parameter p: " << negative_binomial_parameter_p;
        BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] Fitting negative binomial: Negative binomial parameter r: " << negative_binomial_parameter_r;
        BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] Fitting negative binomial: Expected coverage: " << exp_depth_covg;
    } else {
        BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] Insufficient coverage to update error rate";
        exp_depth_covg = fit_mean_covg(kmer_covg_dist, covg / 10);
        exp_depth_covg = std::max(exp_depth_covg, (uint)1);
        BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] Insufficient coverage to update error rate: Expected coverage: " << exp_depth_covg;
    }

    // find probability threshold
    for (const auto& node : pangraph->nodes) {
        node.second->kmer_prg_with_coverage.set_exp_depth_covg(exp_depth_covg);
        if (bin)
            node.second->kmer_prg_with_coverage.set_binomial_parameter_p(e_rate);
        else
            node.second->kmer_prg_with_coverage.set_negative_binomial_parameters(
                negative_binomial_parameter_p, negative_binomial_parameter_r);

        for (uint32_t i = 1;
             i < node.second->kmer_prg_with_coverage.kmer_prg->nodes.size() - 1;
             ++i) // NB first and last kmer in kmergraph are null
        {
            if (bin)
                binomial_parameter_p
                    = node.second->kmer_prg_with_coverage.bin_prob(i, sample_id);
            else
                binomial_parameter_p
                    = node.second->kmer_prg_with_coverage.nbin_prob(i, sample_id);
            for (int j = 0; j < 200; ++j) {
                if ((float)j - 200 <= binomial_parameter_p
                    and (float) j + 1 - 200 > binomial_parameter_p) {
                    kmer_prob_dist[j] += 1;
                    break;
                }
            }
        }
    }

    // save probability distribution
    handle.open(outdir + "/kmer_probs.txt");
    assert(!handle.fail()
        or assert_msg("Could not open file " << outdir + "/kmer_probs.txt"));
    for (int j = 0; (uint32_t)j != kmer_prob_dist.size(); ++j) {
        handle << j - 200 << "\t" << kmer_prob_dist[j] << std::endl;
    }
    handle.close();

    // evaluate threshold
    auto it = kmer_prob_dist.begin();
    while (*it == 0 and it != kmer_prob_dist.end() - 1) {
        ++it;
    }
    ++it;
    // it now represents most negative prob bin with non-zero coverage.
    // if there are enough remaining kmers, estimate thresh, otherwise use a default
    if (std::accumulate(it, kmer_prob_dist.end(), (uint32_t)0) > 1000) {
        thresh = find_prob_thresh(kmer_prob_dist);
        BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] Estimated threshold for true kmers: " << thresh;
    } else {
        thresh = (int)std::distance(kmer_prob_dist.begin(), it) - 200;
        BOOST_LOG_TRIVIAL(info)
            << "[Estimate Parameters] Did not find enough non-zero coverage kmers to estimated threshold for "
               "true kmers. Using the naive threshold: "
            << thresh;
    }

    // set threshold in each kmer graph
    for (auto& node : pangraph->nodes) {
        node.second->kmer_prg_with_coverage.set_thresh(thresh);
    }

    BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] Expected coverage: " << exp_depth_covg;
    BOOST_LOG_TRIVIAL(info) << "[Estimate Parameters] End";
    return exp_depth_covg;
}
