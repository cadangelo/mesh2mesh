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
#include "moab/Util.hpp"
#include "moab/GeomUtil.hpp"
#include "MBTagConventions.hpp"

moab::Core mbi;
moab::Tag flux_tag;

moab::ErrorCode get_mesh_elements(std::string filename,
                                  ){

 moab::ErrorCode rval;

 // Load mesh from file 1 into fileset 1
 rval = mbi.create_meshset(moab::MESHSET_SET, fileset); MB_CHK_ERR(rval);
 rval = mbi.load_file(filename.c_str(), &fileset);
 MB_CHK_SET_ERR(rval, "Error loading file.");
 
 // Get all 3D elements in fileset 1
 rval = mbi.get_entities_by_dimension(fileset, 3, ves);MB_CHK_SET_ERR(rval, "Error getting 3d elements");

 return MB_SUCCESS;
}

int main(int argc, char **argv){

moab::ErrorCode rval;

moab::Range ves;
moab::EntityHandle fileset;

// Get flux tag
std::string flux_tag_name ("flux");
rval = mbi.tag_get_handle(flux_tag_name.c_str(),
//                           moab::MB_TAG_VARLEN,
//                           moab::MB_TYPE_DOUBLE,
                          1,
                          moab::MB_TYPE_INTEGER,
                          flux_tag,
                          moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
MB_CHK_SET_ERR(rval, "Error getting flux tag.");

for (int i = 1; i < argc; ++i){
  std::string filename = argv[i];
  std::cout << filename << std::endl;
  

}


}
