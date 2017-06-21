#include "gtest/gtest.h"
#include "minihit.h"
#include "minimizer.h"
#include "minirecord.h"
#include "interval.h"
#include "path.h"
#include "inthash.h"
#include <stdint.h>
#include <iostream>
#include <algorithm>

using namespace std;

class MinimizerHitTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test
    // (right before the destructor).
  }
};

TEST_F(MinimizerHitTest,create){
    Minimizer* m;
    KmerHash hash;
    pair<uint64_t, uint64_t> kh = hash.kmerhash("ACGTA", 5);
    m = new Minimizer(min(kh.first,kh.second), 0,5,0);
    deque<Interval> d = {Interval(7,8), Interval(10, 14)};
    Path p;
    p.initialize(d);
    MiniRecord* mr;
    mr = new MiniRecord(0,p,0);
    MinimizerHit mh(1, m, mr);
    uint32_t j(1);
    EXPECT_EQ(j, mh.read_id);
    EXPECT_EQ(Interval(0,5), mh.read_interval);
    j=0;
    EXPECT_EQ(j, mh.prg_id);
    EXPECT_EQ(p, mh.prg_path);
    bool b = true;
    EXPECT_EQ(b, mh.strand);

    kh = hash.kmerhash("hell", 4);
    m = new Minimizer(min(kh.first,kh.second),1,5,0);
    EXPECT_DEATH(MinimizerHit(1, m, mr), "");

    delete m;
    delete mr;
    //TEST SECOND CONSTRUCTOR!!
}

TEST_F(MinimizerHitTest,checkStrand){
    Minimizer* m;
    KmerHash hash;
    pair<uint64_t, uint64_t> kh = hash.kmerhash("ACGTA", 5);
    m = new Minimizer(min(kh.first,kh.second), 0,5,0);
    deque<Interval> d = {Interval(7,8), Interval(10, 14)};
    Path p;
    p.initialize(d);
    MiniRecord* mr;
    mr = new MiniRecord(0,p,0);
    MinimizerHit mh(1, m, mr);
    EXPECT_EQ(mh.strand, true);

    delete m;
    delete mr;
    m = new Minimizer(min(kh.first,kh.second), 0,5,1);
    mr = new MiniRecord(0,p,1);
    MinimizerHit mh1(1, m, mr);
    EXPECT_EQ(mh1.strand, true);

    delete m;
    delete mr;
    m = new Minimizer(min(kh.first,kh.second), 0,5,1);
    mr = new MiniRecord(0,p,0);
    MinimizerHit mh2(1, m, mr);
    EXPECT_EQ(mh2.strand, false);
 
    delete m;
    delete mr;
    m = new Minimizer(min(kh.first,kh.second), 0,5,0);
    mr = new MiniRecord(0,p,1);
    MinimizerHit mh3(1, m, mr);
    EXPECT_EQ(mh3.strand, false);

    delete m;
    delete mr;
}

TEST_F(MinimizerHitTest,equals){
    Minimizer* m;
    KmerHash hash;
    pair<uint64_t,uint64_t> kh = hash.kmerhash("ACGTA", 5);
    m = new Minimizer(min(kh.first,kh.second), 0,5,0);
    deque<Interval> d = {Interval(7,8), Interval(10, 14)};
    Path p;
    p.initialize(d);
    MiniRecord* mr;
    mr = new MiniRecord(0,p,0);
    MinimizerHit mh1(1, m, mr);

    delete m;
    m = new Minimizer(min(kh.first,kh.second), 0,5,0);
    d = {Interval(7,9), Interval(11, 14)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0,p,0);
    MinimizerHit mh2(1, m, mr);

    EXPECT_EQ(mh1, mh1);
    EXPECT_EQ(mh2, mh2);
    EXPECT_EQ((mh1==mh2), false);
    delete m;
    delete mr;
}

TEST_F(MinimizerHitTest,compare){
    set<MinimizerHit> hits;
    KmerHash hash;

    Minimizer* m;
    pair<uint64_t,uint64_t> kh = hash.kmerhash("ACGTA", 5);
    m = new Minimizer(min(kh.first,kh.second), 1,6,0);
    deque<Interval> d = {Interval(7,8), Interval(10, 14)};
    Path p;
    p.initialize(d);
    MiniRecord* mr;
    mr = new MiniRecord(0,p,0);
    MinimizerHit mh1(1, m, mr);

    MinimizerHit mh2(0, m, mr);

    delete m;
    m = new Minimizer(min(kh.first,kh.second), 0,5,0);

    d = {Interval(6,10), Interval(11, 12)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0,p,0);
    MinimizerHit mh4(1, m, mr);

    d = {Interval(6,10), Interval(12, 13)};
    p.initialize(d);
    delete mr;
    mr = new MiniRecord(0,p,0);
    MinimizerHit mh5(1, m, mr);

    hits.insert(mh1);
    hits.insert(mh2);
    hits.insert(mh4);
    hits.insert(mh5);
    
    vector<MinimizerHit> expected;
    expected.push_back(mh1);
    expected.push_back(mh2);
    expected.push_back(mh4);
    expected.push_back(mh5);

    uint32_t j(1);
    for (set<MinimizerHit>::iterator it=hits.begin(); it!=--hits.end(); ++it)
    {
        EXPECT_EQ(expected[j], *it);
        j++;
    }
    EXPECT_EQ(expected[0], *(--hits.end()));
    delete m;
    delete mr;
}

