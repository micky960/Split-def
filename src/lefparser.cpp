#include<iostream>
#include<fstream>
#include<cstdlib>
#include<string>
#include<sstream>
#include<vector>
#include<algorithm>
#include<stdlib.h>
#include<cassert>
#include<cmath>
#include<numeric>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/make_shared.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/lexical_cast.hpp>
#include"lefparser.h"

void parse_macro(std::ifstream&, std::unordered_map <std::string, B>& b_list);
void parse_cell_bound(std::ifstream& lef, std::string line, std::unordered_map <std::string, B>& b_list);

//void parse_cur_macro(std::ifstream&, std::string, std::unordered_map <std::string, B>& b_list);
//void get_pin(std::ifstream&, std::unordered_map <std::string, B>& b_list);

void lef_main(std::unordered_map <std::string, B>& b_list){

	std::ifstream lef("../files/NangateOpenCellLibrary.lef");
	std::string line, word;

    std::cout<<"START LEF"<<std::endl;

	do{
		getline(lef, line);
	}while(line.find("END FreePDK45_38x28_10R_NP_162NW_34O")==-1);

	//std::cout<<line<<std::endl;
	//getline(lef, line);
	parse_macro(lef, b_list);
    std::cout<<"END LEF"<<std::endl;
}

//int main(){
//
//	std::ifstream lef("files/NangateOpenCellLibrary.lef");
//	std::string line, word;
//
//	do{
//		getline(lef, line);
//	}while(line.find("END PROPERTYDEFINITIONS")==-1);
//
//	//std::cout<<line<<std::endl;
//	getline(lef, line);
//	parse_macro(lef);
//}

void parse_macro(std::ifstream& lef, std::unordered_map <std::string, B>& b_list){
	std::string word, line;
	getline(lef, line);
	while(line.find("END LIBRARY")==-1){
		getline(lef, line);
		if(line.find("MACRO")!=-1){
			//parse_cur_macro(lef, line);	
			parse_cell_bound(lef, line, b_list);	
		}

	}
	//std::cout<<line<<std::endl;
}

void parse_cell_bound(std::ifstream& lef, std::string line, std::unordered_map <std::string, B>& b_list){
	std::string name;
	std::vector<std::string> tokens;
	boost::split(tokens, line, boost::is_any_of(" "));
	name = tokens[1];
	B b;
	std::cout<<"TYPE:"<<name<<std::endl;

	while(line.find("SIZE")==-1)
		getline(lef, line);
	auto pos = line.find("SIZE")+5;
	line = line.substr(pos, line.size()-pos);	
	boost::split(tokens, line, boost::is_any_of(" , \t"));
	b.first = atof(tokens[0].c_str());
	b.second = atof(tokens[2].c_str());
	b_list[name] = b;
	//std::cout<<line<<std::endl;
	//std::cout<<"SIZE:"<<tokens[0]<<" ,BY:"<<tokens[2]<<std::endl;

}

/*void 	parse_cur_macro(std::ifstream& lef,std::string line){
	std::string name;
	std::vector<std::string> tokens;
	boost::split(tokens, line, boost::is_any_of(" "));
	name = tokens[1];
	std::cout<<tokens[1]<<std::endl;
	getline(lef, line);
	while(line.find("PIN VDD")==-1){
		if(line.find("PIN")!=-1){
			get_pin(lef);
		}
		getline(lef, line);
	}	
}

void get_pin(std::ifstream& lef){
	std::string line;
	std::vector<std::string> tokens;
	getline(lef, line);
	while(line.find("POLYGON")==-1){
		getline(lef, line);
	}
	auto pos = line.find("POLYGON")+8;
	line = line.substr(pos, line.size()-pos);	
	boost::split(tokens, line, boost::is_any_of(" , \t"));
	RECT r;
	//std::cout<<line<<std::endl;
	std::pair<int, int> temp;
	for(int i=0; i<tokens.size(); i=i+2){
		temp.first = atoi(tokens[i].c_str());
		temp.second = atoi(tokens[i+1].c_str());
		r.loc.push_back(temp);
	}
	//std::cout<<tokens[0]<<std::endl;

}*/
