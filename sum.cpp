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

// break string into vector of strings based on delimiter
void tokenize(const std::string& str, 
              std::vector<std::string>& tokens,
              const char* delimiters){
  tokens.clear();

  std::string::size_type next_token_end, next_token_start =
                         str.find_first_not_of( delimiters, 0);

  while ( std::string::npos != next_token_start )
    {
      next_token_end = str.find_first_of( delimiters, next_token_start );
      if ( std::string::npos == next_token_end )
        {
	  tokens.push_back(str.substr(next_token_start));
          next_token_start = std::string::npos;
        }
      else
        {
          tokens.push_back( str.substr( next_token_start, next_token_end -
                                        next_token_start ) );
          next_token_start = str.find_first_not_of( delimiters, next_token_end );
        }
    }
}


moab::ErrorCode get_mesh_elements(std::string filename,
                                  bool blank,
                                  std::map<int, Tet_info> &tet_map,
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
std::vector<double> flux(num_e_groups, 0);

//// Get error tag
//std::string error_tag_name ("photon_result_rel_error");
//rval = mbi.tag_get_handle(error_tag_name.c_str(),
//                           moab::MB_TAG_VARLEN,
//                           moab::MB_TYPE_DOUBLE,
//                           error_tag,
//                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
//std::vector<double> error(num_e_groups, 0);
//std::vector<double> sqrd_error(num_e_groups, 0);


// Get all 3D elements in fileset 1
ves.clear();
rval = mbi.get_entities_by_dimension(fileset, 3, ves);MB_CHK_SET_ERR(rval, "Error getting 3d elements");
MB_CHK_SET_ERR(rval, "Error getting tets.");
int tet_id;
for (it = ves.begin(); it != ves.end(); ++it){
  if (blank == false){
    // get the flux tag on the ve 
    rval = mbi.tag_get_data(flux_tag, &(*it), 1, &flux[0]);//MB_CHK_ERR(rval);
    MB_CHK_SET_ERR(rval, "Error getting flux tag.");
  }
//  // get the error tag on the ve
//  rval = mbi.tag_get_data(error_tag, &(*it), 1, &error[0] ); MB_CHK_ERR(rval);
//  MB_CHK_SET_ERR(rval, "Error getting error tag.");

  // get the id tag on the ve
  rval = mbi.tag_get_data(id_tag, &(*it), 1, &tet_id ); MB_CHK_ERR(rval);
  MB_CHK_SET_ERR(rval, "Error getting id tag.");

  // Add entity handle and flux to tet id map 
  tet_map[tet_id].eh = *it;
  tet_map[tet_id].flux = flux;
//  tet_map[tet_id].error = error;
//  tet_map[tet_id].sqrd_error = sqrd_error;
}

 return moab::MB_SUCCESS;
}

moab::ErrorCode sort_filenames(int num_files, char **files, 
                              std::map<int, std::string> &map_step_to_filename){
  std::string step_num_str;
  int step_num;
  std::string filename;
  std::vector<std::string> tokens;
  const char* delimiter = "_";
  const char* delimiter2 = ".";

  for (int i = 1; i < num_files; ++i){
    filename = files[i];
    tokenize(filename, tokens, delimiter);
    step_num_str = tokens[4];
    tokenize(step_num_str, tokens, delimiter2);
    step_num = std::atoi(tokens[0].c_str());
    map_step_to_filename[step_num] = filename;
  }
}

int main(int argc, char **argv){

std::map<int, std::string> map_step_to_filename;
sort_filenames(argc, argv, map_step_to_filename);

moab::ErrorCode rval;

moab::EntityHandle fileset;
moab::EntityHandle blankmeshset;

std::map<int, Tet_info> tet_flux_map;
std::map<int, Tet_info> tot_flux_map;

// Info about problem that should be read from file or command line
double tot_time = 27.0;
int tot_num_steps = argc - 1;
std::cout << "num steps " << tot_num_steps << std::endl;
double time_per_step = tot_time/tot_num_steps;
int step = 0;
int num_e_groups = 24;

std::map<int, std::string>::iterator step_it;
for (step_it = map_step_to_filename.begin(); step_it != map_step_to_filename.end(); ++step_it){

  // create tet_id map for final mesh
  std::string filename = step_it->second;
  bool blank;
  if (step_it->first == 0){
    blank = true;
    rval = get_mesh_elements(filename, blank, tot_flux_map, blankmeshset);
    MB_CHK_SET_ERR(rval, "Error getting blank mesh file");
  }
  // get flux vals from tet mesh
  tet_flux_map.clear();
  blank = false;
  rval = get_mesh_elements(filename, blank, tet_flux_map, fileset);
  MB_CHK_SET_ERR(rval, "Error getting flux mesh file");
  
  // for each tet ID, keep running total of the flux scored in each configuration 
  moab::EntityHandle tet_id;
  std::map<int, Tet_info>::iterator mit;
  for(mit = tot_flux_map.begin(); mit!=tot_flux_map.end(); ++mit){
    tet_id = mit->first;
    for(int j=0; j <= num_e_groups-1; j++){
      if (step_it->first == 0){
        tot_flux_map[tet_id].flux[j] = time_per_step*tet_flux_map[tet_id].flux[j];
      }
      else{
        tot_flux_map[tet_id].flux[j] += time_per_step*tet_flux_map[tet_id].flux[j];
  //        tot_flux_map[tet_id].sqrd_error[j] += pow(tet_flux_map[tet_id].error[j], 2);
  //        tot_flux_map[tet_id].error[j] = sqrt(tot_flux_map[tet_id].sqrd_error[j]);
      }
    }
    // set summed flux data tag
    rval = mbi.tag_set_data(flux_tag, &(tot_flux_map[tet_id].eh), 1, &tot_flux_map[tet_id].flux[0]);
    MB_CHK_SET_ERR(rval, "Error getting flux tag val.");
//    rval = mbi.tag_set_data(error_tag, &(tot_flux_map[tet_id].eh), 1, &tot_flux_map[tet_id].error[0]);
//    MB_CHK_SET_ERR(rval, "Error getting flux tag val.");
  }
  // Write out summed flux mesh 
  moab::EntityHandle output_list[] = {blankmeshset};
  std::string basename = "_cumulative_flux.h5m";
  std::string filenum = std::to_string(step);
  rval = mbi.write_mesh((filenum+basename).c_str(), output_list, 1);
  MB_CHK_SET_ERR(rval, "Error writing out mesh.");

  //Advance to next time step
  ++step;
}


}
