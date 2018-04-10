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
                                  std::map<moab::EntityHandle, int> &tet_flux_map){

moab::ErrorCode rval;
moab::EntityHandle fileset;
moab::Range ves;
moab::Range::iterator it;
int flux;
//std::map<moab::EntityHandle, int> tet_flux_map;

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

 // Load mesh from file 1 into fileset 1
 rval = mbi.create_meshset(moab::MESHSET_SET, fileset); MB_CHK_ERR(rval);
 rval = mbi.load_file(filename.c_str(), &fileset);
 MB_CHK_SET_ERR(rval, "Error loading file.");
 
 // Get all 3D elements in fileset 1
 ves.clear();
 rval = mbi.get_entities_by_dimension(fileset, 3, ves);MB_CHK_SET_ERR(rval, "Error getting 3d elements");
 MB_CHK_SET_ERR(rval, "Error getting tets.");


 // for each tet, get flux tag
 for (it = ves.begin(); it != ves.end(); ++it){
   // get the flux tag on the ve 
   rval = mbi.tag_get_data(flux_tag, &(*it), 1, &flux);//MB_CHK_ERR(rval);
   MB_CHK_SET_ERR(rval, "Error getting flux tag.");
   
   // add flux to map
   tet_flux_map[*it] = flux;
 }
 std::cout << "flux tag " << flux << std::endl;

 return moab::MB_SUCCESS;
}

int main(int argc, char **argv){

moab::ErrorCode rval;

moab::Range ves;
moab::EntityHandle fileset;

std::map<moab::EntityHandle, int> tet_flux_map;
std::map<moab::EntityHandle, int> tot_flux_map;

int num_steps = 0;
// for each mesh file on run line
for (int i = 1; i < argc; ++i){
  std::string filename = argv[i];
  std::cout << filename << std::endl;


  // get blank mesh vals
  if (i == 1){
    rval = get_mesh_elements(filename, tot_flux_map);
    std::cout << "blank mesh " << tot_flux_map.begin()->second << std::endl;
  }
  else{
  ++num_steps;
  // get flux vals from tet mesh
  rval = get_mesh_elements(filename, tet_flux_map);
    std::cout << "1 mesh " << tet_flux_map.begin()->second << std::endl;
  

  // keep running avg of flux mesh vals
  std::map<moab::EntityHandle, int>::iterator mit;
  moab::EntityHandle tet_id;
  for(mit = tet_flux_map.begin(); mit!=tet_flux_map.end(); ++mit){
    tet_id = mit->first;
    tot_flux_map[tet_id] = (tot_flux_map[tet_id]+tet_flux_map[tet_id])/num_steps;
  }
  std::cout << "tot num steps " << num_steps << std::endl;
  }
}


}
