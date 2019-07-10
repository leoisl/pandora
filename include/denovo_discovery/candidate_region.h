#ifndef PANDORA_CANDIDATE_REGION_H
#define PANDORA_CANDIDATE_REGION_H

#include <vector>
#include <algorithm>
#include <set>
#include "interval.h"
#include "denovo_utils.h"
#include "fastaq.h"
#include "fastaq_handler.h"
#include "localPRG.h"
#include "denovo_discovery.h"
#include <omp.h>


using CandidateRegionIdentifier = std::tuple<Interval, std::string>;


template<>
struct std::hash<CandidateRegionIdentifier> {
    size_t operator()(const CandidateRegionIdentifier &id) const;
};


struct TmpPanNode {
    const PanNodePtr &pangraph_node;
    const std::shared_ptr<LocalPRG> &local_prg;
    const std::vector<KmerNodePtr> &kmer_node_max_likelihood_path;
    const std::vector<LocalNodePtr> &local_node_max_likelihood_path;

    TmpPanNode(const PanNodePtr &pangraph_node, const std::shared_ptr<LocalPRG> &local_prg,
               const std::vector<KmerNodePtr> &kmer_node_max_likelihood_path,
               const std::vector<LocalNodePtr> &local_node_max_likelihood_path);
};


using DenovoPaths = std::vector<std::string>;
using ReadPileup = std::vector<std::string>;
using ReadCoordinates = std::set<ReadCoordinate>;

class DenovoDiscovery;

class CandidateRegion {
public:
    ReadCoordinates read_coordinates;
    std::string max_likelihood_sequence;
    std::string left_flanking_sequence;
    std::string right_flanking_sequence;
    ReadPileup pileup;
    fs::path filename;
    DenovoPaths denovo_paths;

    //variables concerning the GATB graph
    bool graph_files_are_created;
    bool graph_was_correctly_built;
    fs::path gatb_graph_filepath;
    LocalAssemblyGraph graph;
    uint32_t max_path_length;
    double expected_kmer_covg;

    CandidateRegion(const Interval &interval, std::string name);

    CandidateRegion(const Interval &interval, std::string name, const uint_least16_t &interval_padding);

    virtual ~CandidateRegion();

    bool operator==(const CandidateRegion &rhs) const;

    bool operator!=(const CandidateRegion &rhs) const;

    Interval get_interval() const;

    const std::string &get_name() const;

    CandidateRegionIdentifier get_id() const;

    std::string get_max_likelihood_sequence_with_flanks() const;

    void add_pileup_entry(const std::string &read, const ReadCoordinate &read_coordinate);

    void write_denovo_paths_to_file(const fs::path &output_directory) const;

    void create_local_assembly_graph (const DenovoDiscovery &denovo, uint32_t threads = 1);

    //this is a RAM efficiency-related method that can be called **after** all candidate regions pileups were loaded into memory
    //after this happens, all the read pileups are themselves stored,
    //so there is no need to store the read coordinates for these pileups, saving some RAM
    void clear_read_coordinates() { read_coordinates.clear(); }

    //this is a RAM efficiency-related method that can be called **after** CandidateRegion::create_local_assembly_graph() method
    //after the local assembly graph is created for a candidate region, there is no need to store the read pileups for it,
    //so we can remove them and save some RAM
    void clear_pileup() { pileup.clear(); }

    //this is a RAM efficiency-related method that can be called **after** CandidateRegion::write_denovo_paths_to_file() method
    //after writing the de-novo paths, we don't need to store any information about this candidate region, freeing up memory
    //for the others
    void clear();

private:
    const Interval interval;
    const std::string name;
    const uint_least16_t interval_padding;
    //TODO: check if we might have issues here with copy constructor - I think not: OpenMP locks work on the address of the lock I think
    omp_lock_t add_pileup_entry_lock; //synchronizes multithreaded access to the add_pileup_entry() method

    void initialise_filename();

    Fastaq generate_fasta_for_denovo_paths() const;
};


std::vector<Interval> identify_low_coverage_intervals(const std::vector<uint32_t> &covg_at_each_position,
                                                      const uint32_t &min_required_covg = 2,
                                                      const uint32_t &min_length = 1);

using CandidateRegions = std::unordered_map<CandidateRegionIdentifier, CandidateRegion>;

CandidateRegions find_candidate_regions_for_pan_node(const TmpPanNode &pangraph_node_components,
                                                     const uint_least16_t &candidate_region_interval_padding = 0);


using ReadId = uint32_t;
using PileupConstructionMap = std::map<ReadId, std::vector<std::pair<CandidateRegion *, const ReadCoordinate *>>>;

PileupConstructionMap construct_pileup_construction_map(CandidateRegions &candidate_regions);

void
load_all_candidate_regions_pileups_from_fastq(const fs::path &reads_filepath, CandidateRegions &candidate_regions,
                                              const PileupConstructionMap &pileup_construction_map, const uint32_t threads,
                                              const bool clear_read_coordinates /*set this to true if you want to be RAM efficient*/);

void create_all_local_assembly_graphs(CandidateRegions &candidate_regions, const DenovoDiscovery &denovo,
                                      const uint32_t threads, const bool clear_pileup /*set this to true if you want to be RAM efficient*/);

#endif //PANDORA_CANDIDATE_REGION_H
