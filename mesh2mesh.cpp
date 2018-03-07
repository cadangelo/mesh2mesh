#include <iostream>
#include <cstdlib>
#include <cstdio>

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/ElemEvaluator.hpp"
#include "moab/CN.hpp"
#include "moab/SpatialLocator.hpp"


int main(int argc, char **argv){

moab::Core mbi;
moab::ErrorCode rval;

// Load file
rval = mbi.load_file(argv[1]);
MB_CHK_SET_ERR(rval, "Error loading file");

// Get all voxels
moab:: Range elems;
rval = mbi.get_entities_by_dimension(0, 3, elems);
MB_CHK_SET_ERR(rval, "Error getting 3d elements");

// Create a tree to use for the location service
moab::AdaptiveKDTree tree(&mbi);

// Specify an evaluator based on linear hexes
moab::ElemEvaluator el_eval(&mbi);

// Build the SpatialLocator
moab::SpatialLocator sl(&mbi, elems, &tree);


// Get the box extents
  moab::CartVect box_extents, pos;
  moab::BoundBox box = sl.local_box();
  box_extents = box.bMax - box.bMin;

// Query at random places in the tree
  moab::CartVect params;
  int is_inside = 0;
  moab::EntityHandle elem;
  pos = {-119, -119, 1};
  std::cout << "pos " << pos.array()[0] << pos.array()[1] << pos.array()[2] << std::endl;
  rval = sl.locate_point(pos.array(), elem, params.array(), &is_inside, 0.0, 0.0);MB_CHK_ERR(rval);
  std::cout << "elem " << elem << std::endl;
  if (is_inside)
    std::cout << "is inside " << std::endl;

  std::cout << "Mesh contains " << elems.size() << " elements of type "
            << moab::CN::EntityTypeName(mbi.type_from_handle(*elems.begin())) << std::endl;
  std::cout << "Bounding box min-max = (" << box.bMin[0] << "," << box.bMin[1] << "," << box.bMin[2] << ")-("
            << box.bMax[0] << "," << box.bMax[1] << "," << box.bMax[2] << ")" << std::endl;

return 0;

}
