#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <set>
#include <stdint.h>
#include "inthash.h"
#include "minimizer.h"
#include "seq.h"
#include "utils.h"

using std::vector;
using namespace std;

Seq::Seq (uint32_t i, string n, string p, uint32_t w, uint32_t k): id(i), name(n), seq(p) {
    minimizer_sketch (w, k);
}

Seq::~Seq()
{
    for (auto c : sketch)
    {
        delete c;
    }
}

void Seq::initialize(uint32_t i, string n, string p, uint32_t w, uint32_t k)
{
    id = i;
    name = n;
    seq = p;
    sketch.clear();
    minimizer_sketch (w, k);
}

void Seq::minimizer_sketch (const uint32_t w, const uint32_t k)
{
    if (seq.length()+1 < w+k) {//cout << "Sequence too short to sketch" << endl; 
        return;}

    // initializations
    uint c, i;
    uint64_t shift1 = 2 * (k - 1), mask = (1ULL<<2*k) - 1, smallest, kmer[2] = {0,0}, kh[2] = {0,0};
    vector<Minimizer*> vm;
    vm.reserve(w);
    Minimizer* m;

    char myArray[seq.size()+1];//as 1 char space for null is also required
    strcpy(myArray, seq.c_str());

    for(uint32_t buff=0; buff < seq.length() ; ++buff)
    {
	c = nt4(seq[buff]);
        if (c < 4) { // not an ambiguous base
            kmer[0] = (kmer[0] << 2 | c) & mask;           // forward k-mer
            kmer[1] = (kmer[1] >> 2) | (3ULL^c) << shift1; // reverse k-mer	
	    kh[0] = hash64(kmer[0], mask);
            kh[1] = hash64(kmer[1], mask);
	} else {
	    cout << "bad letter - not sure how to handle this so skipping read" << endl;
        }
	
	if (buff >=k-1)
	{
	    //make a mini
	    m = new Minimizer(min(kh[0], kh[1]), buff-k+1, buff+1, (kh[0]<=kh[1]));
	    vm.push_back(m);
	}

	if (vm.size()==w)
	{
	    smallest = std::numeric_limits<uint64_t>::max();
	    // find smallest khash value for the w kmers
	    for (i=0; i!=vm.size(); ++i)
	    {
		smallest = min(vm[i]->kmer, smallest);
	    }
	    // add these minimizers to sketch, and delete others
	    for (i=0; i!=vm.size(); ++i)
            {
		if (vm[i]->kmer == smallest)
		{
		    sketch.insert(vm[i]);
		} else {
		    delete vm[i];
		}
	    }
	    vm.clear();		
	} else if ( buff >= w+k-1 and min(kh[0], kh[1]) <= smallest)
	{
	    sketch.insert(vm.back());
	    // delete all but the last kmer
	    for (i=0; i!=vm.size()-1; ++i)
	    {
		delete vm[i];
	    }
	    vm.clear();
	}

	if (buff == seq.length()-1)
        {
	    // delete remaining elements of vm
	    for (i=0; i!=vm.size(); ++i)
            {
                delete vm[i];
            }
            vm.clear();
	}
    }
    cout << now() << "Sketch size " << sketch.size() << " for read " << name << endl;
    return;
}

/*void Seq::minimizer_sketch (const uint32_t w, const uint32_t k)
{
    //cout << "Start sketching" << endl;
    // If sequence too short, just return
    if (seq.length()+1 < w+k) {//cout << "Sequence too short to sketch" << endl; 
	return;}

    // initializations
    string kmer;
    pair<uint64_t, uint64_t> kh;
    uint64_t smallest;
    Minimizer *m;
    Minimizer *m_previous;
    KmerHash hash;

    // force inclusion of first kmer
    //kmer = seq.substr(0, k);
    //kh = hash.kmerhash(kmer, k);
    //m = new Minimizer(min(kh.first, kh.second), 0, k, 0);
    //sketch.insert(m);
    //m_previous = m;

    // for each window position
    for(uint32_t wpos=0; wpos <= seq.length()-w-k+1 ; ++wpos)
    {
	//cout << "wpos: " << wpos << endl;
    	// if wpos==0 or the previous minimizer starts outside current window, calculate new mini from scratch
        if((wpos == 0) or (m_previous->pos.start < wpos))
	{
            smallest = std::numeric_limits<uint64_t>::max();
	    // find the lex smallest kmer in the window
	    for (uint32_t i = 0; i < w; i++)
	    {
	    	kmer = seq.substr(wpos+i, k);
	        kh = hash.kmerhash(kmer, k);
	        smallest = min(smallest, min(kh.first, kh.second));
	    }
	    for (uint32_t i = 0; i < w; i++)
	    {
	        kmer = seq.substr(wpos+i, k);
	        kh = hash.kmerhash(kmer, k);
	        if (kh.first == smallest or kh.second == smallest)
                {
		    m = new Minimizer(min(kh.first, kh.second), wpos+i, wpos+i+k, (kh.first<=kh.second));
		    sketch.insert(m);
		    //cout << "found from scratch: " << *m << " so sketch size is now " << sketch.size() << endl;
		    m_previous = m;
                }
	    }
        } else {
        // otherwise only need to do something if the kmer from the newest position at end of new window is smaller or equal to the previous smallest
	    kmer = seq.substr(wpos+w-1, k);
            kh = hash.kmerhash(kmer, k);
	    //cout << "End kh for wpos: " << min(kh.first, kh.second) << " compared to previous smallest: " << smallest << endl;
	    if(kh.first <= smallest or kh.second <= smallest)
	    {
	        m = new Minimizer(min(kh.first, kh.second), wpos+w-1, wpos+w-1+k, (kh.first<=kh.second));
		sketch.insert(m);
		m_previous = m;
            }

            smallest = min(smallest, min(kh.first, kh.second));
	}
    }

    // force inclusion of last
    //kmer = seq.substr(seq.length()-k, k);
    //kh = hash.kmerhash(kmer, k);
    //m = new Minimizer(min(kh.first, kh.second), seq.length()-k, seq.length(), 0);
    //if (!(*m_previous==*m))
    //{
    //    sketch.insert(m);
    //}
    
    cout << now() << "Sketch size " << sketch.size() << " for read " << name << endl;
    return;
}*/

std::ostream& operator<< (std::ostream & out, Seq const& data) {
    out << data.name;
    return out ;
}

