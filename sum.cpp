#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <stdio.h>
#include <map>

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
#include "DagMC.hpp"
#include "dagmcmetadata.hpp"

moab::Core mbi;
moab::Tag flux_tag;
moab::Tag error_tag;
moab::Tag id_tag;
moab::Tag move_tag;

moab::DagMC *DAG;
dagmcMetaData* DMD;

struct Tet_info{
       //Tet_info(moab::EntityHandle e=0, std::vector<double> f)
       //:eh(e), flux(f)
       //{}
       moab::EntityHandle eh;
       std::vector<double> flux;
       std::vector<double> error;
       std::vector<double> sqrd_error;
};

int num_groups(moab::Tag tag) {
  moab::ErrorCode rval;
  int tag_size;
  rval = mbi.tag_get_bytes(tag, *(&tag_size));
  if (rval != moab::MB_SUCCESS)
      throw std::runtime_error("Problem getting tag size.");
  return tag_size/sizeof(double);
}
moab::ErrorCode get_mesh_elements(std::string filename,
                                  std::map<int, Tet_info> &tet_map,
                                  //std::map<moab::EntityHandle, int> &tet_flux_map,
                                  //std::map<int, moab::EntityHandle> &id_eh_map,
           //                       int &num_e_groups,
                                  moab::EntityHandle &fileset){

moab::ErrorCode rval;
moab::Range ves;
moab::Range::iterator it;


// Load mesh from file into fileset
rval = mbi.create_meshset(moab::MESHSET_SET, fileset); MB_CHK_ERR(rval);
MB_CHK_SET_ERR(rval, "Error creating meshset.");
rval = mbi.load_file(filename.c_str(), &fileset);
MB_CHK_SET_ERR(rval, "Error loading file.");

// Get ID tag
rval = mbi.tag_get_handle( GLOBAL_ID_TAG_NAME,
                              1, 
                              moab::MB_TYPE_INTEGER,
                              id_tag,
                              moab::MB_TAG_DENSE );
// Get flux tag
std::string flux_tag_name ("photon_result");
rval = mbi.tag_get_handle(flux_tag_name.c_str(),
                           moab::MB_TAG_VARLEN,
                           moab::MB_TYPE_DOUBLE,
                           flux_tag,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
int num_e_groups = num_groups(flux_tag);
std::cout << "num e groups " << num_e_groups << std::endl;
//if (id_only == false){
std::vector<double> flux(num_e_groups, 0);
//}
// Get error tag
std::string error_tag_name ("photon_result_rel_error");
rval = mbi.tag_get_handle(error_tag_name.c_str(),
                           moab::MB_TAG_VARLEN,
                           moab::MB_TYPE_DOUBLE,
                           error_tag,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
std::vector<double> error(num_e_groups, 0);
std::vector<double> sqrd_error(num_e_groups, 0);


// Get all 3D elements in fileset 1
ves.clear();
rval = mbi.get_entities_by_dimension(fileset, 3, ves);MB_CHK_SET_ERR(rval, "Error getting 3d elements");
MB_CHK_SET_ERR(rval, "Error getting tets.");
std::cout << " num elements " << ves.size() << std::endl;

int tet_id;
for (it = ves.begin(); it != ves.end(); ++it){
//  if (id_only == false){
    // get the flux tag on the ve 
    rval = mbi.tag_get_data(flux_tag, &(*it), 1, &flux[0]);//MB_CHK_ERR(rval);
    MB_CHK_SET_ERR(rval, "Error getting flux tag.");
//  }
  // get the error tag on the ve
  rval = mbi.tag_get_data(error_tag, &(*it), 1, &error[0] ); MB_CHK_ERR(rval);
  MB_CHK_SET_ERR(rval, "Error getting error tag.");

  // get the id tag on the ve
  rval = mbi.tag_get_data(id_tag, &(*it), 1, &tet_id ); MB_CHK_ERR(rval);
  MB_CHK_SET_ERR(rval, "Error getting id tag.");

  // Add entity handle and flux to tet id map 
//  tet_id = *it;
  tet_map[tet_id].eh = *it;
  tet_map[tet_id].flux = flux;
  tet_map[tet_id].error = error;
  tet_map[tet_id].sqrd_error = sqrd_error;
}
  std::cout << "tet id gme " << tet_id << std::endl;
// std::cout << "flux in last tet id: " << tet_map[tet_id].flux[42] << std::endl;
 std::cout << "tet map size " << tet_map.size() << std::endl;

 return moab::MB_SUCCESS;
}

int main(int argc, char **argv){

moab::ErrorCode rval;

moab::EntityHandle fileset;
moab::EntityHandle blankmeshset;

std::map<int, Tet_info> tet_flux_map;
std::map<int, Tet_info> tot_flux_map;
//std::map<int, Tet_info> avg_flux_map;
std::map<int, Tet_info>::iterator mit;

int num_steps = 0;
//int num_e_groups = 217;
int num_e_groups = 24;
// for each mesh file on run line
for (int i = 1; i < argc; ++i){
  std::cout << "argc, i " << argc << ", " << i << std::endl; 
  std::string filename = argv[i];
  std::cout << filename << std::endl;

  // create tet_id map for blank/final config mesh
  // change this to string comparison for "blankmesh" or "final_config"
  // id_only == true bc no flux tag on blank meshset
  if (i == 1){
    rval = get_mesh_elements(filename, tot_flux_map, blankmeshset);
    MB_CHK_SET_ERR(rval, "Error getting blank mesh file");
  }
  else{
    // if no time weighting info, will average over number of steps
    ++num_steps;
    // get flux vals from tet mesh
    tet_flux_map.clear();
    rval = get_mesh_elements(filename, tet_flux_map, fileset);
    MB_CHK_SET_ERR(rval, "Error getting flux mesh file");
    std::cout << "tet flux " << tet_flux_map[2].eh << std::endl;;
    
    // for each tet ID, keep running total of the flux scored in each configuration 
    moab::EntityHandle tet_id;
    std::cout << "num tets " << tot_flux_map.size() << std::endl;
    for(mit = tot_flux_map.begin(); mit!=tot_flux_map.end(); ++mit){
      tet_id = mit->first;
      for(int j=0; j <= num_e_groups-1; j++){
        tot_flux_map[tet_id].flux[j] += tet_flux_map[tet_id].flux[j];
//        std::cout << "flux " << tot_flux_map[tet_id].flux[j] << std::endl;;
//        std::cout << "sqrd error tet grp " << tot_flux_map[tet_id].sqrd_error[j] << tet_id << j << std::endl;
        tot_flux_map[tet_id].sqrd_error[j] += pow(tet_flux_map[tet_id].error[j], 2);
        tot_flux_map[tet_id].error[j] = sqrt(tot_flux_map[tet_id].sqrd_error[j]);
      }
      // set summed flux data tag
      rval = mbi.tag_set_data(flux_tag, &(tot_flux_map[tet_id].eh), 1, &tot_flux_map[tet_id].flux[0]);
//      std::cout << "flux " <<  tot_flux_map[tet_id].flux[23] << std::endl;;
      MB_CHK_SET_ERR(rval, "Error getting flux tag val.");
      rval = mbi.tag_set_data(error_tag, &(tot_flux_map[tet_id].eh), 1, &tot_flux_map[tet_id].error[0]);
//      std::cout << "sq error " <<  &tot_flux_map[tet_id].sqrd_error[0] << std::endl;;
//      std::cout << "error " <<  &tot_flux_map[tet_id].error[0] << std::endl;;
      MB_CHK_SET_ERR(rval, "Error getting flux tag val.");
    }
    // Write out summed flux mesh 
    moab::EntityHandle output_list[] = {blankmeshset};
    std::string basename = "_cumulative_flux.h5m";
    std::string filenum = std::to_string(num_steps);
  rval = mbi.write_mesh((filenum+basename).c_str(), output_list, 1);
  MB_CHK_SET_ERR(rval, "Error writing out mesh.");
  }
}


}
