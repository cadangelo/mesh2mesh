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
#include "DagMC.hpp"
#include "dagmcmetadata.hpp"

moab::Core mbi;
moab::Tag flux_tag;
moab::Tag id_tag;
moab::Tag move_tag;

moab::DagMC *DAG;
dagmcMetaData* DMD;

struct Tet_info{
       //Tet_info(moab::EntityHandle e=0, std::vector<double> f)
       //:eh(e), flux(f)
       //{}
       moab::EntityHandle eh;
       //int flux;
       std::vector<double> flux;
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
                                  moab::EntityHandle &fileset){

moab::ErrorCode rval;
//moab::EntityHandle fileset;
moab::Range ves;
moab::Range::iterator it;
//int flux;
//std::map<moab::EntityHandle, int> tet_flux_map;

// Get flux tag
std::string flux_tag_name ("flux");
//rval = mbi.tag_get_handle(flux_tag_name.c_str(),
////                           moab::MB_TAG_VARLEN,
////                           moab::MB_TYPE_DOUBLE,
//                          1,
//                          moab::MB_TYPE_INTEGER,
//                          flux_tag,
//                          moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
//MB_CHK_SET_ERR(rval, "Error getting flux tag.");
rval = mbi.tag_get_handle(flux_tag_name.c_str(),
                           moab::MB_TAG_VARLEN,
                           moab::MB_TYPE_DOUBLE,
                           flux_tag);
int num_e_groups = num_groups(flux_tag);
std::vector<double> flux(num_e_groups);

rval = mbi.tag_get_handle( GLOBAL_ID_TAG_NAME,
                              1, 
                              moab::MB_TYPE_INTEGER,
                              id_tag,
                              moab::MB_TAG_DENSE );

 // Load mesh from file 1 into fileset 1
 rval = mbi.create_meshset(moab::MESHSET_SET, fileset); MB_CHK_ERR(rval);
 rval = mbi.load_file(filename.c_str(), &fileset);
 MB_CHK_SET_ERR(rval, "Error loading file.");


 // Get all 3D elements in fileset 1
 ves.clear();
 rval = mbi.get_entities_by_dimension(fileset, 3, ves);MB_CHK_SET_ERR(rval, "Error getting 3d elements");
 MB_CHK_SET_ERR(rval, "Error getting tets.");

int tet_id;
//std::map<int, Tet_info> tet_map;
 // for each tet, get flux tag
 for (it = ves.begin(); it != ves.end(); ++it){
   // get the flux tag on the ve 
   rval = mbi.tag_get_data(flux_tag, &(*it), 1, &flux[0]);//MB_CHK_ERR(rval);
   MB_CHK_SET_ERR(rval, "Error getting flux tag.");
   // get the id tag on the ve
   rval = mbi.tag_get_data(id_tag, &(*it), 1, &tet_id ); MB_CHK_ERR(rval);
   MB_CHK_SET_ERR(rval, "Error getting id tag.");

   // add flux to map
   //tet_flux_map[*it] = flux;
   //tet_flux_map[tet_id] = flux;

   // map id to eh
   //id_eh_map[tet_id] = *it;

   // test struct value
   //Tet_info tet; 
   tet_map[tet_id].eh = *it;
   tet_map[tet_id].flux = flux;
 }
 std::cout << "flux tag " << tet_id <<  flux << std::endl;
 std::cout << "flux tag from struct map " << tet_id <<  tet_map[tet_id].flux << std::endl;


 return moab::MB_SUCCESS;
}

//moab::ErrorCode add_motion_vec_to_geom(){
//  
//  moab::ErrorCode rval;
////  rval = mbi.tag_get_handle( "MOVE_TAG", 32, moab::MB_TYPE_OPAQUE,
////			      move_tag, moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
////  CHECK_ERR(rval);
//
//}

//moab::ErrorCode new_get_tagged_entities( 
//                                        std::string tag_name, 
//                                        std::map<int, moab::Range> &tr_vols_map)
//{
//  moab::ErrorCode rval;
// 
//  DAG = new moab::DagMC(mbi);
//  rval = DAG->load_existing_contents();
//  CHECK_ERR(rval);
//std::map<int, moab::Range> tr_vols_map;
//int num_cells = DAG->num_entities( 3 );
//
//
//  // get moving tag from parse_properties
//  std::string tag_delims = ":";
//  std::vector<std::string> group_name;
//  std::map<std::string, std::string> group_name_synonyms;
//
//  group_name.push_back(tag_name);
//
//  rval = DAG->parse_properties(group_name, group_name_synonyms, tag_delims.c_str());
//  if (moab::MB_SUCCESS != rval) 
//    {
//      std::cerr << "DAGMC failed to parse metadata properties" <<  std::endl;
//      exit(EXIT_FAILURE);
//    }
//
//  // get desired tagged entities
//  moab::EntityHandle tagged_meshset;
//  moab::Range surf_set, vert_set;
//  int num_verts;
//  moab::Range::iterator it, itr;
//
//  rval = mbi.create_meshset(moab::MESHSET_SET, tagged_meshset);
//  MB_CHK_SET_ERR(rval, "Error creating meshset.");
//  //CHECK_ERR(rval);
//
//  std::string val;
//  for( int i = 1; i <= num_cells; ++i ) 
//    { 
//      std::vector<std::string> properties;
// 
//      moab::EntityHandle vol = DAG->entity_by_index( 3, i );
//
//      if( DAG->has_prop( vol, tag_name))
//        {
//          rval = DAG->prop_values(vol, tag_name, properties);
//
//          std::cout << "num trs on vol " << properties.size() << std::endl;
//          // get the value of vol's TR number (val) 
//          rval = DAG->prop_value(vol, tag_name, val);
//          MB_CHK_SET_ERR(rval, "Error getting prop val.");
//          std::cout << "val " << vol << " " << val << std::endl;
//
//          std::map<int, moab::Range>::iterator itt;
//          bool added = false;
//          int tr_num;
//          //if the map is not empty, look for a key that matches the tr # of the volume
//          if(tr_vols_map.size() > 0){
//            for(itt = tr_vols_map.begin(); itt != tr_vols_map.end(); ++itt){
//                
//               //get value of tr #
//               tr_num = itt->first;
//               std::cout << "tr num already in map" << tr_num << std::endl;
//               std::cout << "val num " << val << std::endl; 
//               //if tr# of current vol matches one already in map, add vol to range 
//               // added = true as soon as vol added to a range
//               if(stoi(val) == tr_num){
//                 tr_vols_map[tr_num].insert(vol); 
//                 std::cout << "size of range for that tr #" << tr_vols_map[tr_num].size() << std::endl;
//                 std::cout << "first elemet" << *tr_vols_map[tr_num].begin() << std::endl;
//                 added = true;
//               }
//               // if added not true, keep looping, if true exit loop
//               if(added == true)
//                 break;
//            }
//          }
//          //if map empty or no matching key found, create key
//          //and add vol to range
//          //if(tr_vols_map.size() == 0 | added == false) {
//          if(added == false) {
//            tr_num = stoi(val);
//            tr_vols_map[tr_num].insert(vol); 
//            std::cout << "size of range for that tr #" << tr_vols_map[tr_num].size() << std::endl;
//            std::cout << "first elemet" << *tr_vols_map[tr_num].begin() << std::endl;
//            added = true;
//    
//          }
//          
//        }
//          
//    }
//  delete DAG;
//  return moab::MB_SUCCESS;
//}

int main(int argc, char **argv){

moab::ErrorCode rval;

moab::Range ves;
moab::EntityHandle fileset;

//std::map<moab::EntityHandle, int> tet_flux_map;
//std::map<moab::EntityHandle, std::vector<double>> tet_flux_map;
//std::map<moab::EntityHandle, std::vector<double>> tot_flux_map;
//std::map<moab::EntityHandle, std::vector<double>> avg_flux_map;
//std::map<int, moab::EntityHandle> blankid_eh_map;
//std::map<int, moab::EntityHandle> fluxid_eh_map;

// integer flux tag for quick testing
std::string flux_tag_name ("flux");
//rval = mbi.tag_get_handle(flux_tag_name.c_str(),
////                           moab::MB_TAG_VARLEN,
////                           moab::MB_TYPE_DOUBLE,
//                          1,
//                          moab::MB_TYPE_INTEGER,
//                          flux_tag,
//                          moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
//MB_CHK_SET_ERR(rval, "Error getting flux tag.");
rval = mbi.tag_get_handle(flux_tag_name.c_str(),
                           moab::MB_TAG_VARLEN,
                           moab::MB_TYPE_DOUBLE,
                           flux_tag);
int num_e_groups = num_groups(flux_tag);
std::vector<double> groupwise_flux(num_e_groups);


int num_steps = 0;
  moab::EntityHandle tet_id;

moab::EntityHandle blankmeshset;

//rval = new_get_tagged_entities("tr", tr_vols_map);
//MB_CHK_SET_ERR(rval, "Error getting tr tag.");

std::map<int, Tet_info> tet_flux_map;
std::map<int, Tet_info> blank_tet_map;
std::map<int, Tet_info> tot_flux_map;
std::map<int, Tet_info>::iterator mit;

// for each mesh file on run line
for (int i = 1; i < argc; ++i){
  std::string filename = argv[i];
  std::cout << filename << std::endl;
//  std::cout << "new tot " << tet_id << " " << tot_flux_map[tet_id] << std::endl;


  // get blank mesh vals
  // change this to string comparison for "blankmesh" or "final_config"
  if (i == 1){
    //rval = get_mesh_elements(filename, tot_flux_map, blankid_eh_map, blankmeshset);
    rval = get_mesh_elements(filename, tot_flux_map, blankmeshset);
    MB_CHK_SET_ERR(rval, "Error getting blank mesh file");
//    std::cout << "blank mesh " << tot_flux_map.begin()->second << std::endl;
  }
  else{
  // if no time weighting info, will average over number of steps
  ++num_steps;
  // get flux vals from tet mesh
  tet_flux_map.clear();
  //rval = get_mesh_elements(filename, tet_flux_map, fluxid_eh_map, fileset);
  rval = get_mesh_elements(filename, tet_flux_map, fileset);
  MB_CHK_SET_ERR(rval, "Error getting flux mesh file");
  
  // for each tet ID, keep running total of the flux scored in each configuration 
  for(mit = tot_flux_map.begin(); mit!=tot_flux_map.end(); ++mit){
    tet_id = mit->first;
    for(int i=0; i=num_e_groups-1; i++){
      tot_flux_map[tet_id].flux[i] += tet_flux_map[tet_id].flux[i];
    }
    // if this is the last mesh file, find the avg (simplest is just dividing by number of configs)
    if(i == argc-1){
      //int avg;
      std::vector<double> avg(num_e_groups);
      for(int i=0; i=num_e_groups-1; i++){
        avg[i] = (tot_flux_map[tet_id].flux[i])/num_steps;
      }
//      std::cout << "this is the avg " << avg << std::endl;
      
      //set the flux tag val to be the avg 
      //rval = mbi.tag_set_data(flux_tag, &(*mit), 1, &(avg_flux_map[tet_id]));
      //rval = mbi.tag_set_data(flux_tag, &(blankid_eh_map[tet_id]), 1, &avg);
      //rval = mbi.tag_set_data(flux_tag, &(blankid_eh_map[tet_id]), 1, &avg[0]);
      rval = mbi.tag_set_data(flux_tag, &(blank_tet_map[tet_id].eh), 1, &avg[0]);
      MB_CHK_SET_ERR(rval, "Error setting flux tag.");
      //rval = mbi.tag_get_data(flux_tag, &(blankid_eh_map[tet_id]), 1, &avg);
      //rval = mbi.tag_get_data(flux_tag, &(blankid_eh_map[tet_id]), 1, &avg[0]);//MB_CHK_ERR(rval);
      rval = mbi.tag_get_data(flux_tag, &(blank_tet_map[tet_id].eh), 1, &avg[0]);//MB_CHK_ERR(rval);
      MB_CHK_SET_ERR(rval, "Error getting flux tag val.");
//      std::cout << "from tag:this is the avg " << avg << std::endl;
    }
  }
//  std::cout << "tot num steps " << num_steps << std::endl;
  }
//  std::cout << "= new tot " << tet_id << " " << tot_flux_map[tet_id] << std::endl;
}

// Write out mesh 2 w/ mapped data
moab::EntityHandle output_list[] = {blankmeshset};
rval = mbi.write_mesh("avgflux.h5m", output_list, 1);

}
