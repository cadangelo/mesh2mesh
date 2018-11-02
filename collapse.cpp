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

const char* delimiter = "_";
const char* delimiter2 = "-";

moab::Range::iterator it;
for( it = ves1.begin(); it != ves1.end(); ++ it){
  std::cout << "el # " << *it << std::endl;
  std::vector<moab::Tag> tag_handles;
  rval = mbi.tag_get_tags_on_entity(*it, tag_handles);
  std::vector<double> collected_ebounds;
  std::vector<moab::Tag>::iterator itv;
  for (itv = tag_handles.begin(); itv != tag_handles.end(); ++itv){
    std::string tag_name;
    rval = mbi.tag_get_name(*itv, tag_name);
    // Tokenize the line
    std::vector<std::string> tokens;
    tokenize(tag_name, tokens, delimiter);
    if(tokens[0] == "TALLY" && tokens.size() > 2){
      std::vector<std::string> ebounds;
      tokenize(tokens[2], ebounds, delimiter2);
      double ebound;
      if(ebounds.size() > 2){
       ebound = sciToDub(ebounds[0]+"-"+ebounds[1]);
      }
      else{
        ebound = std::atof(ebounds[0].c_str());
      }
      collected_ebounds.push_back(ebound);
    }
  }
  std::sort(collected_ebounds.begin(), collected_ebounds.end());
  std::vector<double>::iterator its;
  for(its = collected_ebounds.begin(); its != collected_ebounds.end(); ++its){
    std::cout << "sorted " << *its << std::endl;
  }
}

}
