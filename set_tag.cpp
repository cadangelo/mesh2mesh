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

int main(int argc, char **argv){

moab::ErrorCode rval;

moab::Range ves;
moab::EntityHandle fileset;

std::string filename = argv[1];
std::cout << filename << std::endl;
// Load mesh from file 1 into fileset 1
rval = mbi.create_meshset(moab::MESHSET_SET, fileset); MB_CHK_ERR(rval);
rval = mbi.load_file(filename.c_str(), &fileset);
MB_CHK_SET_ERR(rval, "Error loading file.");

// Get all 3D elements in fileset 1
rval = mbi.get_entities_by_dimension(fileset, 3, ves);MB_CHK_SET_ERR(rval, "Error getting 3d elements");

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
int flux = std::stoi(argv[2]);
//int flux =1 ;

moab::Range::iterator it;
for( it = ves.begin(); it != ves.end(); ++ it){
     // set the flux tag on the ve from set 2 we are mapping to
     rval = mbi.tag_set_data(flux_tag, &(*it), 1, &flux);MB_CHK_ERR(rval);
}

// Write out mesh 2 w/ mapped data
moab::EntityHandle output_list[] = {fileset};
std::string filenum = std::to_string(flux);
std::string basename = "taggeddata.h5m";
rval = mbi.write_mesh((filenum+basename).c_str(), output_list, 1);



}
