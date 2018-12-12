#include <iostream>
#include <boost/program_options.hpp>
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

namespace po = boost::program_options;

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

//moab::ErrorCode setup(std::string file1, std::string file2, 
moab::ErrorCode setup(int argc, char* argv[], 
                      moab::Range &ves1,
                      moab::Range &ves2,
                      moab::EntityHandle &fileset2,
                      moab::EntityHandle &set1, 
                      bool &expand,
                      std::string &map_tag_name){
  // process command line flags
  po::options_description desc("Allowed options");
  desc.add_options()
      ("help", "Help message")
      ("expand", po::value<bool>(), "Expand vector tags?")
      ("mesh1", po::value<std::string>(), "Mesh with values that will be mapped onto second mesh file")
      ("mesh2", po::value<std::string>(), "Second mesh file that will have values mapped onto it")
      ("tag", po::value<std::string>(), "Name of tag to map from mesh 1 to mesh 2")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);   

  std::string file1, file2; 
  
  if (vm.count("help")) {
      std::cout << desc << std::endl;
      return moab::MB_FAILURE;
  }
  
  if (vm.count("expand")) {
      std::cout << "Expand tags? " << vm["expand"].as<bool>() << std::endl;
      expand = vm["expand"].as<bool>();
  } else {
      std::cout << "Expand flag not set." << std::endl;;
  }
  if (vm.count("mesh1")) {
      std::cout << "Mesh 1 file name " << vm["mesh1"].as<std::string>() << std::endl;
      file1 =  vm["mesh1"].as<std::string>();
  } else {
      std::cout << "Mesh 1 file not set." << std::endl;;
  }
  if (vm.count("mesh2")) {
      std::cout << "Mesh 2 file name " << vm["mesh2"].as<std::string>() << std::endl;
      file2 =  vm["mesh2"].as<std::string>();
  } else {
      std::cout << "Mesh 2 file not set." << std::endl;;
  }
  if (vm.count("tag")) {
      std::cout << "Tag name " << vm["tag"].as<std::string>() << std::endl;
      map_tag_name =  vm["tag"].as<std::string>();
  } else {
      std::cout << "Tag name not set." << std::endl;;
  }

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
bool expand; //expand tags
std::string map_tag_name;
//rval = setup(argv[1], argv[2], ves1, ves2, fileset2, set1, expand);
rval = setup(argc, argv, ves1, ves2, fileset2, set1, expand, map_tag_name);
std::cout << "num vol elements in file 1 " << ves1.size() << std::endl;
std::cout << "num vol elements in file 2 " << ves2.size() << std::endl;
std::cout << "map tag name " << map_tag_name << std::endl;


// Create a tree to use for the location service
moab::AdaptiveKDTree tree(&mbi);

// Specify an evaluator based on linear hexes
moab::ElemEvaluator el_eval(&mbi);

// Build the SpatialLocator
moab::SpatialLocator sl(&mbi, ves1, &tree);


// Get flux tag
//std::string flux_tag_name ("flux");
//std::string flux_tag_name ("source_density");
//std::string flux_tag_name ("photon_result");
//std::string source_tag_name ("source_density");
moab::Tag flux_tag;
moab::Tag source_tag;
rval = mbi.tag_get_handle(map_tag_name.c_str(),
                           moab::MB_TAG_VARLEN,
                           moab::MB_TYPE_DOUBLE,
                           flux_tag,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
//rval = mbi.tag_get_handle(source_tag_name.c_str(),
//                           moab::MB_TAG_VARLEN,
//                           moab::MB_TYPE_DOUBLE,
//                           source_tag,
//                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
int num_e_groups = num_groups(flux_tag);
std::vector<double> groupwise_flux(num_e_groups);
std::vector<double> zero_flux(num_e_groups);
std::fill (zero_flux.begin(), zero_flux.end(), 0);

// for each ve in set 2, find centroid
moab::EntityHandle leaf; // set of elements that contains point
moab::CartVect params;
int is_inside = 0;
moab::EntityHandle ve1;
unsigned int num;
moab::Range::iterator it;
std::string scalar_map_tag_name;
moab::Tag one_group_flux;
std::vector<moab::Tag> single_group_flux_tags;
for( it = ves2.begin(); it != ves2.end(); ++ it){
     // first, set the flux tag on the ve from set 2 to zero 
     rval = mbi.tag_set_data(flux_tag, &(*it), 1, &zero_flux[0]);
     MB_CHK_SET_ERR(rval, "Could not set flux tag to zero on 'map-to' mesh");

   //Coord centroid = find_centroid(*it);
   moab::CartVect centroid = find_centroid(*it);
   // find ve in set 1 that encloses centroid point
   rval = sl.locate_point(centroid.array(), leaf, params.array(), &is_inside);MB_CHK_ERR(rval);
   if (is_inside){
     rval = el_eval.find_containing_entity(leaf, centroid.array(), 1e-1, 1e-1, ve1, params.array(), &num);MB_CHK_ERR(rval);

     if ( ve1 == 0 ){
       std::cout << "Warning, can't find underlying mesh element"<< std::endl; 
     }
     else{
       // get the flux tag on the ve from set 1 enclosing the point
       rval = mbi.tag_get_data(flux_tag, &ve1, 1, &groupwise_flux[0]);//MB_CHK_ERR(rval);
       MB_CHK_SET_ERR(rval, "Could not get flux tag");
//       std::cout << "ve1 el, setting data " << *it << ", " << groupwise_flux[0] << std::endl;
   
       if(expand){
       //expand the vector tag to a scalar tag
       for(int grp=0; grp <= num_e_groups-1; grp++){
         //scalar_map_tag_name = "flux"+std::to_string(grp);
         scalar_map_tag_name = map_tag_name+std::to_string(grp);
         rval = mbi.tag_get_handle(scalar_map_tag_name.c_str(), 1, 
                                   moab::MB_TYPE_DOUBLE, one_group_flux,
                                   moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
         MB_CHK_SET_ERR(rval, "Could not get flux tag handle");
//         single_group_flux_tags[grp] = groupwise_flux[grp];
         rval = mbi.tag_set_data(one_group_flux, &(*it), 1, &(groupwise_flux[grp]));//MB_CHK_ERR(rval);
         MB_CHK_SET_ERR(rval, "Could not set flux tag");
       }
       }

       // set the flux tag on the ve from set 2 we are mapping to
       rval = mbi.tag_set_data(flux_tag, &(*it), 1, &groupwise_flux[0]);//MB_CHK_ERR(rval);
       MB_CHK_SET_ERR(rval, "Could not set flux tag");
//       std::cout << "ve2 el, setting data " << *it << ", " << groupwise_flux[0] << std::endl;
     }
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
std::cout << "writing mapped data" << std::endl;
rval = mbi.write_mesh("mappeddata.h5m", output_list, 1);

return 0;

}
