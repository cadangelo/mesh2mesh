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

moab::Core mbi;

struct Coord{

  double x;
  double y;
  double z;
};

int num_groups(moab::Tag tag) {
  moab::ErrorCode rval;
  int tag_size;
  rval = mbi.tag_get_bytes(tag, *(&tag_size));
  if (rval != moab::MB_SUCCESS)
      throw std::runtime_error("Problem getting tag size.");
  return tag_size/sizeof(double);
}

moab::ErrorCode setup(std::string file1, std::string file2, 
                      moab::Range &ves1,
                      moab::Range &ves2,
                      moab::EntityHandle &fileset2,
                      moab::EntityHandle &set1){
moab::ErrorCode rval;

// create sets for mesh elements 
moab::EntityHandle fileset1;
rval = mbi.create_meshset(moab::MESHSET_SET, fileset1); MB_CHK_ERR(rval);
rval = mbi.create_meshset(moab::MESHSET_SET, fileset2); MB_CHK_ERR(rval);
rval = mbi.create_meshset(moab::MESHSET_SET, set1); MB_CHK_ERR(rval);

// Load mesh from file 1 into fileset 1
rval = mbi.load_file(file1.c_str(), &fileset1);
MB_CHK_SET_ERR(rval, "Error loading file 1");
// Load mesh from file 2 into fileset 2
rval = mbi.load_file(file2.c_str(), &fileset2);
MB_CHK_SET_ERR(rval, "Error loading file 2");

// Get all 3D elements in fileset 1
rval = mbi.get_entities_by_dimension(fileset1, 3, ves1);MB_CHK_SET_ERR(rval, "Error getting 3d elements");
rval = mbi.add_entities(set1, ves1);
// Get all 3D elements in fileset 2
rval = mbi.get_entities_by_dimension(fileset2, 3, ves2);MB_CHK_SET_ERR(rval, "Error getting 3d elements");

return moab::MB_SUCCESS;

}

//Coord find_centroid(moab::EntityHandle ve){
//
//moab::ErrorCode rval;
//const moab::EntityHandle *connectivity;
//int number_nodes = 0;
//rval = mbi.get_connectivity(ve, connectivity, number_nodes,true);
//Coord coord;
//
//coord.x=0.0;
//coord.y=0.0;
//coord.z=0.0;
//
//for(int i = 0; i< number_nodes; i++)
//{
//   double node_coords[3];
//   rval = mbi.get_coords(&(connectivity[i]), 1, node_coords);
//  
//   coord.x+=node_coords[0];
//   coord.y+=node_coords[1];
//   coord.z+=node_coords[2];
//}
//
//coord.x/=(double)number_nodes;
//coord.y/=(double)number_nodes;
//coord.z/=(double)number_nodes;
//
//return coord;
//}

moab::CartVect find_centroid(moab::EntityHandle ve){

moab::ErrorCode rval;
const moab::EntityHandle *connectivity;
int number_nodes = 0;
rval = mbi.get_connectivity(ve, connectivity, number_nodes,true);
moab::CartVect coord;

coord.array()[0]=0.0;
coord.array()[1]=0.0;
coord.array()[2]=0.0;

for(int i = 0; i< number_nodes; i++)
{
   double node_coords[3];
   rval = mbi.get_coords(&(connectivity[i]), 1, node_coords);
  
   coord.array()[0]+=node_coords[0];
   coord.array()[1]+=node_coords[1];
   coord.array()[2]+=node_coords[2];
}

coord.array()[0]/=(double)number_nodes;
coord.array()[1]/=(double)number_nodes;
coord.array()[2]/=(double)number_nodes;

return coord;
}


int main(int argc, char **argv){

moab::ErrorCode rval;

// Setup loads files and populates volume element ranges
moab::Range ves1, ves2;
moab::EntityHandle fileset2, set1;
rval = setup(argv[1], argv[2], ves1, ves2, fileset2,set1);
std::cout << "num ves2 " << ves2.size() << std::endl;
std::cout << "num ves1 " << ves1.size() << std::endl;


// Create a tree to use for the location service
moab::AdaptiveKDTree tree(&mbi);

// Specify an evaluator based on linear hexes
moab::ElemEvaluator el_eval(&mbi);

// Build the SpatialLocator
moab::SpatialLocator sl(&mbi, ves1, &tree);


// Get flux tag
std::string flux_tag_name ("flux");
moab::Tag flux_tag;
rval = mbi.tag_get_handle(flux_tag_name.c_str(),
                           moab::MB_TAG_VARLEN,
                           moab::MB_TYPE_DOUBLE,
                           flux_tag);
int num_e_groups = num_groups(flux_tag);
std::vector<double> groupwise_flux(num_e_groups);


// for each ve in set 2, find centroid
//moab::CartVect point; // centroid from ve in set 2 to find in set 1
moab::EntityHandle leaf; // set of elements that contains point
moab::CartVect params;
int is_inside = 0;
moab::EntityHandle ve1;
unsigned int num;
moab::Range::iterator it;
for( it = ves2.begin(); it != ves2.end(); ++ it){
   //Coord centroid = find_centroid(*it);
   moab::CartVect centroid = find_centroid(*it);
 //  point.array()[0] = centroid.x;
 //  point.array()[1] = centroid.y;
 //  point.array()[2] = centroid.z;
   // find ve in set 1 that encloses centroid point
   rval = sl.locate_point(centroid.array(), leaf, params.array(), &is_inside);MB_CHK_ERR(rval);
   if (is_inside){
     rval = el_eval.find_containing_entity(leaf, centroid.array(), 1e-4, 1e-6, ve1, params.array(), &num);MB_CHK_ERR(rval);

     // get the flux tag on the ve from set 1 enclosing the point
     rval = mbi.tag_get_data(flux_tag, &ve1, 1, &groupwise_flux[0]);//MB_CHK_ERR(rval);

     // set the flux tag on the ve from set 2 we are mapping to
     rval = mbi.tag_set_data(flux_tag, &(*it), 1, &groupwise_flux[0]);//MB_CHK_ERR(rval);
  }// is_inside
}

// Delete idx tag
//  (needed for expand_tags.py to work properly)
std::string idx_tag_name ("idx");
moab::Tag idx_tag;
rval = mbi.tag_get_handle(idx_tag_name.c_str(),
                          idx_tag);
rval = mbi.tag_delete(idx_tag);

// Write out mesh 2 w/ mapped data
moab::EntityHandle output_list[] = {fileset2};
rval = mbi.write_mesh("mappeddata.h5m", output_list, 1);

return 0;

}
