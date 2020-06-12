#ifndef PANDORA_GCPWRAPPER_H
#define PANDORA_GCPWRAPPER_H

#include "sampleinfo.h"
#include "GCP.h"
#include "OptionsAggregator.h"
#include "KmerCoverageModel.h"

// forward declarations
class SampleInfo;


/**
 * Class that teaches the library how to produce simulated data for each sample.
 */
class GCPSampleInfoModel : public GCP::Model<SampleInfo> {
public:
    /**
     * Builds a GCPSampleInfoModel
     * @param sample_index : Required to build a SampleInfo later.
     * @param genotyping_options :  Required to build a SampleInfo later.
     * @param kmer_coverage_model : Kmer coverage model used to draw kmer coverages from.
     */
    GCPSampleInfoModel(uint32_t sample_index,
                       GenotypingOptions const* genotyping_options,
                       std::shared_ptr<KmerCoverageModel> kmer_coverage_model) :
        Model<SampleInfo>(), sample_index(sample_index), genotyping_options(genotyping_options),
        kmer_coverage_model(kmer_coverage_model) {}

    SampleInfo produce_data() override;

private:
    uint32_t sample_index;
    GenotypingOptions const* genotyping_options;
    std::shared_ptr<KmerCoverageModel> kmer_coverage_model;
};

/**
 * Container that stores and manages a vector of GCPSampleInfoModels, one per sample.
 */
class GCPSampleInfoModels : public std::vector<std::shared_ptr<GCPSampleInfoModel>>{
public:
    GCPSampleInfoModels(const KmerCoverageModels &kmer_coverage_models, GenotypingOptions const* genotyping_options);
};


/**
 * Adapts Pandora's Genotyper to the genotyper needed by the GCP library
 */
class GCPGenotyperAdapter {
public:
    GCPGenotyperAdapter (const SampleInfo &sample_info) : sample_info(sample_info) {}
    double get_genotype_confidence() const;

private:
    const SampleInfo &sample_info;
};

/**
 * Abstracts the GCP library functionality into this class.
 * Allows to get the confidence percentile for genotype confidences for a given GCPSampleInfoModel
 */
class GCPWrapper {
public:
    GCPWrapper(const std::shared_ptr<GCPSampleInfoModel> &model) : model(model){
        if (model != nullptr) {
            train();
        }

    }

    virtual GCP::GenotypePercentile get_confidence_percentile(GCP::GenotypeConfidence queried_confidence) const {
        return genotype_confidence_percentiler->get_confidence_percentile(queried_confidence);
    }

private:
    std::shared_ptr<GCPSampleInfoModel> model;
    std::shared_ptr<GCP::Percentiler> genotype_confidence_percentiler;

    void train();
};

/**
 * Container that stores and manages a vector of GCPWrappers, one per sample.
 */
class GCPWrappers : public std::vector<GCPWrapper>{
public:
    GCPWrappers (const GCPSampleInfoModels &sample_info_models);
};


#endif // PANDORA_GCPWRAPPER_H
