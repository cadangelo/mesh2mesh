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


//int num_groups(moab::Tag tag){
//  moab::ErrorCode rval;
//  int tag_size;
//  rval = mbi.tag_get_bytes(tag, *(&tag_size));
//  if (rval != moab::MB_SUCCESS)
//      throw std::runtime_error("Problem getting tag size.");
//  return tag_size/sizeof(double);
//}

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

// Get flux tag
std::string flux_tag_name ("flux");
moab::Tag flux_tag;
rval = mbi.tag_get_handle(flux_tag_name.c_str(),
                           moab::MB_TAG_VARLEN,
                           moab::MB_TYPE_DOUBLE,
                           flux_tag);
std::cout << "rval " << rval << std::endl;
//num_e_groups = num_groups(src_tag);
//  std::vector<double> pdf(num_ves*num_e_groups); 
//  rval = mesh->tag_get_data(src_tag, ves, &pdf[0]);

// Query point
  moab::CartVect params;
  int is_inside = 0;
  moab::EntityHandle elem;
  double data;
  pos = {-119, -119, 1};
  std::cout << "pos " << pos.array()[0] << pos.array()[1] << pos.array()[2] << std::endl;
  rval = sl.locate_point(pos.array(), elem, params.array(), &is_inside, 0.0, 0.0);MB_CHK_ERR(rval);
  std::cout << "elem " << elem << std::endl;
  rval = mbi.tag_get_data(flux_tag, &elem, 1, &data);
  std::cout << "tag data " << data << std::endl;

  if (is_inside)
    std::cout << "is inside " << std::endl;

  std::cout << "Mesh contains " << elems.size() << " elements of type "
            << moab::CN::EntityTypeName(mbi.type_from_handle(*elems.begin())) << std::endl;
  std::cout << "Bounding box min-max = (" << box.bMin[0] << "," << box.bMin[1] << "," << box.bMin[2] << ")-("
            << box.bMax[0] << "," << box.bMax[1] << "," << box.bMax[2] << ")" << std::endl;



return 0;

}
