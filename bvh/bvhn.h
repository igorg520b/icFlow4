#ifndef BVHN_H
#define BVHN_H

#include <vector>
#include "kdop8.h"
#include "ConcurrentPool.h"

namespace icy { class BVHN; }

class icy::BVHN
{
public:
    static ConcurrentPool<std::vector<BVHN*>> VectorFactory;
    static ConcurrentPool<BVHN> BVHNFactory;

    kDOP8 box;
    bool isLeaf;
    bool test_self_collision;   // can disable self-collision tests on fragments
    int level;
    std::pair<unsigned,unsigned> feature; // if leaf, refers to the feature that is enveloped by this kDOP

    BVHN();
    void Build(std::vector<BVHN*> *bvs, int level_);
    void Update();
    void SelfCollide(std::vector<std::pair<unsigned,unsigned>> &broad_list);
    void Collide(BVHN *b, std::vector<std::pair<unsigned,unsigned>> &broad_list);

private:
    BVHN *child1, *child2;
};

#endif // BVHN_H
