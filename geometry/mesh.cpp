#include "mesh.h"
#include "model.h"
#include <numeric>
#include <algorithm>
#include <iterator>

// code from icy::Geometry class that changes less often

#include <bits/stdc++.h>

icy::Mesh::Mesh(double CharacteristicLengthMax)
{

}

void icy::Mesh::Reset(double CharacteristicLengthMax)
{
    elems.clear();
    nodes.clear();

    // invoke Gmsh
}

