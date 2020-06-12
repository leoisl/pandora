#ifndef PANDORA_OPTIONSAGGREGATOR_H
#define PANDORA_OPTIONSAGGREGATOR_H

#include <vector>
#include <cstdint>

class GenotypingOptions {
private:
    std::vector<uint32_t> sample_index_to_exp_depth_covg;
    double error_rate;
    uint16_t confidence_threshold;
    uint32_t min_allele_covg;
    double min_fraction_allele_covg;
    uint32_t min_site_total_covg;
    uint32_t min_site_diff_covg;
    bool snps_only;

public:
    GenotypingOptions(const std::vector<uint32_t>& sampleIndexToExpDepthCovg,
        double errorRate, uint16_t confidenceThreshold, uint32_t minAlleleCovg,
        double minFractionAlleleCovg, uint32_t minSiteTotalCovg,
        uint32_t minSiteDiffCovg, bool snpsOnly)
        : sample_index_to_exp_depth_covg(sampleIndexToExpDepthCovg)
        , error_rate(errorRate)
        , confidence_threshold(confidenceThreshold)
        , min_allele_covg(minAlleleCovg)
        , min_fraction_allele_covg(minFractionAlleleCovg)
        , min_site_total_covg(minSiteTotalCovg)
        , min_site_diff_covg(minSiteDiffCovg)
        , snps_only(snpsOnly)
    {
    }

    const std::vector<uint32_t>& get_sample_index_to_exp_depth_covg() const
    {
        return sample_index_to_exp_depth_covg;
    }

    double get_error_rate() const { return error_rate; }

    uint16_t get_confidence_threshold() const { return confidence_threshold; }

    uint32_t get_min_allele_covg() const { return min_allele_covg; }

    double get_min_fraction_allele_covg() const { return min_fraction_allele_covg; }

    uint32_t get_min_site_total_covg() const { return min_site_total_covg; }

    uint32_t get_min_site_diff_covg() const { return min_site_diff_covg; }

    uint32_t get_min_kmer_covg(uint32_t sample_index) const {
        return sample_index_to_exp_depth_covg[sample_index] / 10;
    }

    bool is_snps_only() const { return snps_only; }

    void add_exp_depth_covg(uint32_t exp_depth_covg)
    {
        sample_index_to_exp_depth_covg.push_back(exp_depth_covg);
    }
};

#endif // PANDORA_OPTIONSAGGREGATOR_H
