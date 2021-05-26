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
    int level;
    unsigned featureIdx; // if leaf, refers to the feature that is enveloped by this kDOP

    BVHN();
    void Initialize(std::vector<BVHN*> *bvs, int level_);
    void Update();
    void UpdateLeaf();  // use tentative coordinates from elem
    void SelfCollide(std::vector<unsigned> &broad_list);
    void Collide(BVHN *b, std::vector<unsigned> &broad_list);

private:
    BVHN *child1, *child2;
};

#endif // BVHN_H
