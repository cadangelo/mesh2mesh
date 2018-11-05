#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"

moab::Core mbi;

moab::ErrorCode setup(std::string file1, 
                      moab::Range &ves1,
                      moab::EntityHandle &set1){
moab::ErrorCode rval;

// create sets for mesh elements 
moab::EntityHandle fileset1;
rval = mbi.create_meshset(moab::MESHSET_SET, fileset1); MB_CHK_ERR(rval);
rval = mbi.create_meshset(moab::MESHSET_SET, set1); MB_CHK_ERR(rval);

// Load mesh from file 1 into fileset 1
rval = mbi.load_file(file1.c_str(), &fileset1);
MB_CHK_SET_ERR(rval, "Error loading file 1");

// Get all 3D elements in fileset 1
rval = mbi.get_entities_by_dimension(fileset1, 3, ves1);MB_CHK_SET_ERR(rval, "Error getting 3d elements");
rval = mbi.add_entities(set1, ves1);

return moab::MB_SUCCESS;
}

double sciToDub(const std::string& str){
  std::istringstream ss(str);
  double d = 0;
  ss >> d;
  if (ss.fail()) {
     std::string s = "Unable to format ";
     s += str;
     s += " as a number!";
     throw (s);
  }
  return (d);
}

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

//std::map<double, int> ebound_to_group_map

int main(int argc, char **argv){

moab::ErrorCode rval;

// Setup loads files and populates volume element ranges
moab::Range ves1;
moab::EntityHandle fileset2, set1;
rval = setup(argv[1], ves1, set1);
std::cout << "num ves1 " << ves1.size() << std::endl;

moab::Tag flux_tag;
std::string flux_tag_name ("n_flux");
rval = mbi.tag_get_handle(flux_tag_name.c_str(),
                           175,//moab::MB_TAG_VARLEN,
                           moab::MB_TYPE_DOUBLE,
                           flux_tag,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
const char* delimiter = "_";
const char* delimiter2 = "-";

std::vector<double> collected_ebounds;
std::map<int, double> group_to_ebound;

moab::Range::iterator it;
for( it = ves1.begin(); it != ves1.end(); ++ it){
  std::cout << "el # " << *it << std::endl;
  std::vector<moab::Tag> tag_handles;
  rval = mbi.tag_get_tags_on_entity(*it, tag_handles);
  std::vector<moab::Tag>::iterator itv;
  std::map<double, double> ebound_to_data_map;
  for (itv = tag_handles.begin(); itv != tag_handles.end(); ++itv){
    std::string tag_name;
    rval = mbi.tag_get_name(*itv, tag_name);
    // Tokenize the line
    std::vector<std::string> tokens;
    tokenize(tag_name, tokens, delimiter);
    // Find Tally tags whose names end w/ energy bounds
    if(tokens[0] == "TALLY" && tokens.size() > 2){
      // get lower energy bound
      std::vector<std::string> ebounds;
      tokenize(tokens[2], ebounds, delimiter2);
      double lebound; // lower energy bound
      if(ebounds.size() > 2){  //if > 2, number in sci-not so also split by - sign
       lebound = sciToDub(ebounds[0]+"-"+ebounds[1]);
      }
      else{
        lebound = std::atof(ebounds[0].c_str());
      }
      if (it == ves1.begin()){
        collected_ebounds.push_back(lebound);
      }
      // get tally data (scalar flux tag)
      double single_group_flux;
      rval = mbi.tag_get_data(*itv, &(*it), 1, &single_group_flux);
      MB_CHK_SET_ERR(rval, "Could not get single group flux");
      // map lower e_bound to tally data
      ebound_to_data_map[lebound] = single_group_flux;
      std::cout << "single grp flux, lebound " << single_group_flux << ", " << lebound << std::endl;
      // delete scalar tag
      //rval = mbi.tag_delete(*itv);
      //MB_CHK_SET_ERR(rval, "Could not delete scalar tag");
    }
  }
  if (it == ves1.begin()){
    std::sort(collected_ebounds.begin(), collected_ebounds.end());
    std::vector<double>::iterator its;
    int i = 1;
    for(its = collected_ebounds.begin(); its != collected_ebounds.end(); ++its){
      group_to_ebound[i] = *its;
      std::cout << "grp to bound " << i << ", " << group_to_ebound[i] << std::endl;
      ++i;
    }
//    std::cout << "group 1, 175  "<< group_to_ebound[0] << ", " << group_to_ebound[175] << std::endl;
  
  }
  // create vector tag
  int num_e_groups = 175;//collected_ebounds.size();
  std::vector<double> groupwise_flux(num_e_groups);
  std::vector<double> get_groupwise_flux(num_e_groups);
  for (int j = 0 ; j < num_e_groups; ++j){
    double ebound = group_to_ebound[j+1];
    groupwise_flux[j] = ebound_to_data_map[ebound];
  }
  rval = mbi.tag_set_data(flux_tag, &(*it), 1, &groupwise_flux[0]);//MB_CHK_ERR(rval);
  MB_CHK_SET_ERR(rval, "Could not set vector flux tag");
  rval = mbi.tag_get_data(flux_tag, &(*it), 1, &get_groupwise_flux[0]);//MB_CHK_ERR(rval);
  MB_CHK_SET_ERR(rval, "Could not set vector flux tag");
  int grp = 25;
  std::cout << "vector flux, size" << get_groupwise_flux[grp] << ", " << get_groupwise_flux.size() << std::endl;
  int length;
  rval = mbi.tag_get_length(flux_tag, length);
  std::cout << "length" << length << std::endl;


}
rval = mbi.write_mesh("collapsed.h5m");

}
