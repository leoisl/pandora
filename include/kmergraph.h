#ifndef __KMERGRAPH_H_INCLUDED__   // if kmergraph.h hasn't been included yet...
#define __KMERGRAPH_H_INCLUDED__

class LocalPRG;

#include <cstdint>
#include <vector>
#include <iostream>
#include <set>
#include "prg/path.h"
#include "kmernode.h"
#include "pangenome/ns.cpp"


struct condition {
    prg::Path q;

    condition(const prg::Path &);

    bool operator()(const KmerNodePtr) const;
};

struct pCompKmerNode {
    bool operator()(KmerNodePtr, KmerNodePtr);
};


class KmerGraph {
private:
    uint32_t reserved_size;
    uint32_t k;

public:
    uint32_t shortest_path_length;
    std::vector<KmerNodePtr> nodes;
    std::set<KmerNodePtr, pCompKmerNode> sorted_nodes; // representing ordering of the nodes compatible with dp

    KmerGraph();

    KmerGraph(const KmerGraph &);

    KmerGraph &operator=(const KmerGraph &);

    virtual ~KmerGraph() = default;

    void clear();

    KmerNodePtr add_node(const prg::Path &);

    KmerNodePtr add_node_with_kh(const prg::Path &, const uint64_t &, const uint8_t &num = 0);

    void add_edge(KmerNodePtr, KmerNodePtr);

    void remove_shortcut_edges();

    void check() const;

    void discover_k();

    uint32_t min_path_length();

    void save(const std::string &, const std::shared_ptr<LocalPRG> = nullptr);
    void load(const std::string &);

    bool operator==(const KmerGraph &y) const;

    friend std::ostream &operator<<(std::ostream &out, KmerGraph const &data);

    friend uint32_t
    estimate_parameters(std::shared_ptr<pangenome::Graph>, const std::string &, const uint32_t, float &, const uint32_t,
                        bool &, const uint32_t &sample_id);


    //friends
    friend struct condition;
    friend class KmerGraphWithCoverage;
};

#endif