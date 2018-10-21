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
#include <list>
#include <omp.h>
#include <boost/algorithm/string.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/make_shared.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/adapted/boost_polygon.hpp>
#include "lefparser.h"

static int cut_nets = 0;

//Define a class that has the data you want to associate to every vertex and edge
struct node{
    node() : is_pin(false), is_via(false), is_via_sep(false), is_open(true) {}
	std::string name;
	int layer;
	int x,y;
	bool is_pin;
	bool is_via;
	bool is_via_sep;
	bool is_open;
	std::string via_line;
};

struct edge{
    edge(): via(false), valid(false) {}
	bool valid;
	bool via;
	int layer;
	double x1,y1, x2, y2;
	std::string metal;
};

struct pin{
    pin(): is_io(false) {}
	std::string name;
	int layer;
	int x,y;
	int x1, y1, x2, y2, x3, y3;
	bool is_io;
};

struct conn{
    double x,y,r;
    std::string p_net;
    std::string dir;
    bool operator<(const conn& a) const { return x < a.x; }
};

typedef struct node NODE;
typedef struct edge EDGE;
typedef struct pin PIN;
typedef struct conn CONN;

//maintain a list of primary IO
std::unordered_map <std::string, PIN> pinlist;

struct cell{
	std::string name_key, type;
	int x1, y1, x2, y2, x3, y3;
	double x,y;
};

std::unordered_map <std::string, B> b_list;

std::ofstream out;
omp_lock_t writelock;

typedef struct cell CELL;
std::unordered_map <std::string, CELL> complist;

typedef std::pair<std::string, std::vector<std::string> > mypair;
std::vector<mypair> netlist, net_lastline;
int split;
std::vector<std::string> keep_layer, track_lines;
typedef std::pair <std::string, std::string> pin_pair;
typedef std::pair <std::string, int> metal_pair;
std::vector<CONN> center;

//Define the graph using above classes
typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, NODE, EDGE > Graph;
//Some typedefs for simplicity
typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t;
typedef boost::graph_traits<Graph>::edge_descriptor edge_t;

// Iterators
typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter;
typedef boost::graph_traits<Graph>::edge_iterator edge_iter;
typedef boost::graph_traits<Graph>::adjacency_iterator adjacency_it;
typedef boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;
//typedef boost::property_map<Graph, boost::edge_index_t>::type EdgeIndexMap;

typedef boost::shared_ptr<std::vector<unsigned long>> vertex_component_map;

struct EdgeInComponent
{ 
    vertex_component_map mapping_;
    unsigned long which_;
    Graph const& master_;

    EdgeInComponent(vertex_component_map m, unsigned long which, Graph const& master) 
        : mapping_(m), which_(which), master_(master) {}

    template <typename Edge> bool operator()(Edge const&e) const {
        return mapping_->at(source(e,master_))==which_
            || mapping_->at(target(e,master_))==which_;
    } 
};

struct VertexInComponent
{ 
    vertex_component_map mapping_;
    unsigned long which_;

    VertexInComponent(vertex_component_map m, unsigned long which)
        : mapping_(m), which_(which) {}

    template <typename Vertex> bool operator()(Vertex const&v) const {
        return mapping_->at(v)==which_;
    } 
};

struct AnyVertex { 
    template <typename Vertex> bool operator()(Vertex const&) const { return true; }
};

typedef boost::filtered_graph<Graph, EdgeInComponent, VertexInComponent> ComponentGraph;

std::vector<ComponentGraph> connected_components_subgraphs(Graph const&g)
{
    vertex_component_map mapping = boost::make_shared<std::vector<unsigned long>>(num_vertices(g));
    size_t num = boost::connected_components(g, mapping->data());

    std::vector<ComponentGraph> component_graphs;

    for (size_t i = 0; i < num; i++)
        component_graphs.push_back(ComponentGraph(g, EdgeInComponent(mapping, i, g), VertexInComponent(mapping, i)));

    return component_graphs;
}

void parse_init(std::ifstream& def);//parse initial lines
void parse_track(std::ifstream& def);//parse track lines
void parse_till_comp(std::ifstream& def);//parse rest lines
void parse_via(std::ifstream& def);//parse via lines
void get_complist(std::ifstream& def);//get components location
void get_pins(std::ifstream& def, int split);//get all the pins 
void get_spnets(std::ifstream& def);//get all the special nets
void get_nets(std::ifstream& def);//get all the nets
void process_net(mypair cur_net, mypair lastline);//process individual net
void parse_till_end(std::ifstream& def);//parse rest lines
void layer_to_keep(int layer);//which layer to keep
void create_graph(std::vector<metal_pair> metal_lines, std::vector<pin_pair> pinlist, mypair lastline);//create graph for a net
void compute_center();
std::vector<pin_pair> pin_inst(std::vector<std::string> pins);//get the pins of an individual net
void print_pin(mypair cur_net, mypair lastline);//print the current pin

int main(int argc, char* argv[]){
	
	std::ifstream def(argv[1]);

    if(def.fail()){
        std::cerr << "Def file not found"<<std::endl;
        exit(1);
    }
	out.open(argv[2]);
    if(out.fail()){
        std::cerr << "Cannot open out file"<<std::endl;
        exit(1);
    }
	std::string line, buf;
	split = atoi(argv[3]);
	
	std::vector<std::string> all_lines;

	assert(argc!=3);

    //parse the lef file to get cell bounds
	lef_main(b_list);	
	
    //which layers to keep
	layer_to_keep(split);
	//parse till component as it is
    std::cout << "parsig till components"<<std::endl;
    parse_till_comp(def);

    //get component list    
    std::cout << "parsig components"<<std::endl;
	get_complist(def);
	
	//get all the IO pins
    std::cout << "parsig pins"<<std::endl;
    get_pins(def, split);
    //omit the special nets; do not process
    std::cout << "parsig spnets"<<std::endl;
	get_spnets(def);
    //parse all the nets
    std::cout << "parsig nets"<<std::endl;
	get_nets(def);

    std::cout << "start splitting"<<std::endl;
    omp_init_lock(&writelock);
    //#pragma omp target
    //#pragma omp parallel for schedule(dynamic)
	for(unsigned i=1; i<netlist.size(); i++){
	    //std::cout<<"#NET "<<i<<std::endl;
	    std::cout<<"NET "<<netlist[i].first<<std::endl;
		process_net(netlist[i], net_lastline[i]);
    }

	//std::cout<<"NET "<<netlist[68].first<<std::endl;
	//process_net(netlist[68], net_lastline[68]);
	
    out<<"END NETS"<<std::endl;
    //parse till end as it is
    parse_till_end(def);

    //std::cout << "# CUT NETS:\t" << cut_nets<<std::endl;
    //compute_center();
    
    omp_destroy_lock(&writelock);

}

void parse_init(std::ifstream& def){
	
	std::string line;
	getline(def, line);

	while(line.find("TRACKS")==-1){
		out<<line<<std::endl;
		getline(def, line);

	}
    //only keep those tracks that are below split layer	
	for(unsigned j=0; j<keep_layer.size(); j++){
		if(line.find(keep_layer[j])!=-1){	
			out<<line<<std::endl;
			track_lines.push_back(line);
			break;
		}
	}
}

void parse_track(std::ifstream& def){
	
	std::string line;
	do{
		getline(def, line);
		for(unsigned j=0; j<keep_layer.size(); j++){
			if(line.find(keep_layer[j])!=-1){
				out<<line<<std::endl;
				track_lines.push_back(line);
				break;
			}
		}
	}while(line.find("TRACKS")!=-1);
}

void parse_till_comp(std::ifstream& def){
	
	std::string line;
	getline(def, line);
	while(line.find("COMPONENTS")==-1){
		out<<line<<std::endl;
		getline(def, line);
	}

    //std::cout<<line<<std::endl;
	out << line <<std::endl;
}

void get_complist(std::ifstream& def){

	std::vector<std::string> tokens;
	CELL temp;
	std::string xy, line;
	while(line.find("END COMPONENTS")==-1){//parse all complist

        //std::cout<<line<<std::endl;
		out<<line<<std::endl;
		temp = {};	
		getline(def, line);
		if(line.find(";")!=-1||line.find("END COMPONENTS")!=-1)
			continue;
	
		boost::split(tokens, line, boost::is_any_of("-"));
		std::stringstream ss;	
		ss<<tokens[1];
		ss>>temp.name_key;//get name of comp
		ss<<tokens[1];
		ss>>temp.type;//type of comp
		boost::split(tokens, line, boost::is_any_of("("));
		boost::split(tokens, tokens[1], boost::is_any_of(" "));
		
		temp.x = atof((tokens[1]).c_str());
		temp.y = atof((tokens[2]).c_str());
		
		
		boost::split(tokens, line, boost::is_any_of(" "));
		std::string dir=tokens[tokens.size()-1];
		//std::cout<<"DIR:"<<tokens[tokens.size()-1]<<std::endl;
        //left bottom corner
		temp.x1 = temp.x;
		temp.y1 = temp.y;
        //right bottom corner
		temp.x2 = temp.x + b_list[temp.type].first*2000;
		temp.y2 = temp.y;
        //left top corner
		temp.x3 = temp.x;
		temp.y3 = temp.y + b_list[temp.type].second*2000;
	
		//std::cout<<"NAME:"<<temp.name_key<<std::endl;
		//std::cout<<"TYPE:"<<temp.type;
		//std::cout<<"X1,Y1: "<<temp.x1<<", "<<temp.y1<<std::endl;
		//std::cout<<"X2,Y2: "<<temp.x2<<", "<<temp.y3<<std::endl;
		
		complist[temp.name_key] = temp;	
	}
	
    //std::cout<<line<<std::endl;
	out<<line<<std::endl;

}

void get_pins(std::ifstream& def, int split){
	
	std::string line, name, dir, x, y, x1, y1, x2, y2;
	std::vector<std::string> tokens;
	PIN temp;

	getline(def, line);
	//std::cout<<line<<std::endl;
	out<<line<<std::endl;
	getline(def, line);
	//std::cout<<line<<std::endl;
    out<<line<<std::endl;
	getline(def, line);
    std::string line1(line);
	//std::cout<<line<<std::endl;
    
	while(line.find("END PINS")==-1){
        out<<line<<std::endl;
		boost::split(tokens, line, boost::is_any_of(" "));
		name = tokens[1];
        if(name == "VDD" || name == "VSS"){
            getline(def, line);
            out<<line<<std::endl;
            getline(def, line);
            continue;
        }
		//std::cout<<name<<std::endl;
		getline(def, line);
        out<<line<<std::endl;
        std::string line2(line);
		std::stringstream ss(line);
        std::string layer;
        ss>>layer;
        ss>>layer;
        ss>>layer;//pin layer
		//std::cout<<"PIN layer:"<<layer.substr(5,layer.size()-5)<<std::endl;
        temp.layer = atoi((layer.substr(5,layer.size()-5)).c_str());
		boost::split(tokens, line, boost::is_any_of("("));
		ss.str(tokens[1]);
		ss>>x1;
		ss>>y1;
		ss.str(tokens[2]);
		ss>>x2;
		ss>>y2;
		getline(def, line);
        out<<line<<std::endl;
        std::string line3(line);
		//std::cout<<atoi(x1.c_str())<<" "<<atoi(y1.c_str())<<std::endl;
		//std::cout<<x1<<" "<<y2<<std::endl;
		boost::split(tokens, line, boost::is_any_of("("));
		ss.str(tokens[1]);
		//std::cout<<tokens[1]<<std::endl;
        //get position of pin
		ss>>x;
		temp.x = atoi(x.c_str());
		ss>>y;
		temp.y = atoi(y.c_str());
		boost::split(tokens, line, boost::is_any_of(")"));
		ss.str(tokens[1]);
		ss>>dir;
		//std::cout<<dir<<std::endl;
		if(dir == "N"){	temp.y = temp.y + (atoi(y1.c_str()) + atoi(y2.c_str()))/2;}
		else if(dir == "S"){ temp.y = temp.y - (atoi(y1.c_str()) + atoi(y2.c_str()))/2;}
		else if(dir == "E"){ temp.x = temp.x + atoi(x2.c_str());}
		else if(dir == "W"){ temp.x = temp.x + atoi(x1.c_str());}
		
		temp.name = name;
		temp.is_io = true;
		std::string key = "PIN "+name;
		//std::cout<<name<<std::endl;
        if(temp.layer <= split || temp.name.find("clk")!=-1 || temp.name.find("CLK")!=-1){
		    pinlist[key] = temp;	
	        //std::cout<<line1<<std::endl<<line2<<std::endl<<line3<<std::endl;
	        //out<<line1<<std::endl<<line2<<std::endl<<line3<<std::endl;
        }
		getline(def, line);
        line1 = line;
		//std::cout<<line<<std::endl;
	}
	//std::cout<<line<<std::endl;
	out<<line<<std::endl;
}

void get_spnets(std::ifstream& def){
	
	std::string line;

	do{
		getline(def, line);
	    //std::cout<<line<<std::endl;
	    //out<<line<<std::endl;
		
	}while(line.find("END SPECIALNETS")==-1);

	//std::cout<<line<<std::endl;
}

void get_nets(std::ifstream& def){

	std::vector<std::string> tokens;
	std::string line, net;

	getline(def, line);
	//std::cout<<line<<std::endl;
	out<<line<<std::endl;
	getline(def, line);
	//std::cout<<line<<std::endl;
	out<<line<<std::endl;
	
	while(line.find("END NETS")==-1){
		mypair temp, temp2;
		//std::cout<<line<<std::endl;
		temp.first = line;
		temp2.first = line;
		while(line.find(";")==-1){
            if(line.find("+")==-1 || line.find("ROUTED")!=-1){
                //std::cout<<line<<std::endl;
			    temp.second.push_back(line);
            }
            else{
                temp2.second.push_back(line);   
            }
			getline(def, line);
		}
		netlist.push_back(temp);//contains the metal lines
        net_lastline.push_back(temp2);//contains the last lines FREQ, PROP etc.
		getline(def, line);
		
	}

}

void process_net(mypair cur_net, mypair lastline){

	std::string line;
	std::vector<std::string> temp, pin_lines, omit_metal_layer;
	metal_pair mp;
	std::vector<metal_pair> metal_lines;	

    if(cur_net.first.find("clk")!=-1){
        print_pin(cur_net, lastline);
        return;
    }
        
	unsigned i;
	for(i=1; i<cur_net.second.size(); i++){//omit first line; it is net name
		if(cur_net.second[i].find("ROUTED")!=-1)
			break;
		else	pin_lines.push_back(cur_net.second[i]);
	}
	for(; i<cur_net.second.size(); i++){
		bool got_split = true;
		for(unsigned j=0; j<keep_layer.size(); j++)
			if(cur_net.second[i].find(keep_layer[j])!=-1){
				//std::cout<<cur_net.second[i]<<std::endl;
				mp.first = cur_net.second[i];
				mp.second = j+1;
				metal_lines.push_back(mp);	
				got_split = false;
				break;
			}
		if(got_split){//if gets splitted, store which metal layers got deleted
			omit_metal_layer.push_back(cur_net.second[i]);
			//std::cout<<"INVALID"<<std::endl;
		}
        if(cur_net.second[i].find("via"+boost::lexical_cast<std::string>(split))!=std::string::npos){
            omit_metal_layer.push_back(cur_net.second[i]);
        }
	}
	std::vector<pin_pair> pinpinlist;
	//std::cout<<"INVALID SIZE: "<<omit_metal_layer.size()<<std::endl;
	//std::cout<<"PIN: "<<cur_net.first<<std::endl;
	
	if(omit_metal_layer.size()==0){//if does not get split then print current net as it is
		print_pin(cur_net, lastline);
        //omp_set_lock(&writelock);
        //cut_nets++;
        //omp_unset_lock(&writelock);
    }
			
	else	{	
		//std::cout<<"OMIT SIZE: "<<omit_metal_layer.size()<<std::endl;
		//std::cout<<"KEPT SIZE: "<<metal_lines.size()<<std::endl;
		pinpinlist = pin_inst(pin_lines);
        std::cout<<cur_net.first<<std::endl;
		create_graph(metal_lines, pinpinlist, lastline);
		//std::cout<<"CALLED"<<std::endl;
	}
}

void layer_to_keep(int layer){

	for(unsigned i=1; i<=layer; i++){
		keep_layer.push_back("metal"+boost::lexical_cast<std::string>(i)+" ");
		//std::cout<<"metal"+boost::lexical_cast<std::string>(i)<<std::endl;;
	}

}

std::vector<pin_pair> pin_inst(std::vector<std::string> pins){
	std::vector<std::string> tokens;
	std::string word;
	std::vector<pin_pair> pinpinlist;
	pin_pair temp;
	for(unsigned i=0; i<pins.size(); i++){
		boost::split(tokens, pins[i], boost::is_any_of("("));
		for(unsigned j=1; j<tokens.size(); j++){
			//std::cout<<tokens[j]<<std::endl;
			std::stringstream ss(tokens[j]);
			ss>>word;
			//std::cout<<word<<std::endl;
			temp.first = word;
			ss>>word;
			temp.second = word;
			//std::cout<<word<<std::endl;
			pinpinlist.push_back(temp);		
		}
	}
	return pinpinlist;
}

void create_graph(std::vector<metal_pair> metal_lines, std::vector<pin_pair> pinpinlist, mypair lastline){
	std::vector<std::string> tokens;
	int x1,x2,y1,y2;
	Graph g;
	std::string word;
	NODE via;
    //for(int i=0; i<metal_lines.size(); i++)
    //    std::cout<<metal_lines[i].first<<std::endl;
    std::unordered_map <std::string, NODE> nodelist, vialist;
	//std::cout<<"KEPT SIZE: "<<metal_lines.size()<<std::endl;
	for(unsigned i=0; i<metal_lines.size(); i++){//pre-process to delete via info
		if(metal_lines[i].first.find("via"+boost::lexical_cast<std::string>(split))!=-1){
			int pos = metal_lines[i].first.find("via"+boost::lexical_cast<std::string>(split));
			metal_lines[i].first = metal_lines[i].first.substr(0, pos);
			//std::cout<<metal_lines[i].first<<std::endl;
		}
	}
	for(unsigned i=0; i<metal_lines.size(); i++){//process to get location
		//std::cout<<metal_lines[i].first<<std::endl;
		boost::split(tokens, metal_lines[i].first, boost::is_any_of("("));
		//std::cout<<tokens.size()<<std::endl;
		if(tokens.size()==3){//metal info
			boost::split(tokens, metal_lines[i].first, boost::is_any_of("("));
			std::stringstream ss(tokens[1]);
			ss>>word;
			x1 = atoi(word.c_str());
			ss>>word;
			y1 = atoi(word.c_str());
			boost::split(tokens, tokens[2], boost::is_any_of(" "));
			if(tokens[1]=="*") x2 = x1;
			else	x2 = atoi(tokens[1].c_str());
			if(tokens[2]=="*") y2 = y1;
			else	y2 = atoi(tokens[2].c_str());
			//std::cout<<"METAL LAYERS: "<<metal_lines[i].second<<std::endl;
			// Create two vertices in that graph
			vertex_t u = boost::add_vertex(g);
			vertex_t v = boost::add_vertex(g);
            //join by metal edge
			edge_t e; bool b;
			boost::tie(e,b) = boost::add_edge(u,v,g);
			g[u].x = x1; g[u].y = y1; g[u].layer = metal_lines[i].second;
			g[v].x = x2; g[v].y = y2; g[v].layer = metal_lines[i].second;
			g[e].x1 = x1; g[e].y1 = y1; g[e].x2 = x2; g[e].y2 = y2; g[e].layer = metal_lines[i].second; g[e].valid = true; g[e].metal = metal_lines[i].first;
			if(metal_lines[i].first.find("via")!=-1){//if via present at the last
				boost::split(tokens, metal_lines[i].first, boost::is_any_of(" "));
				//std::cout<<"VIA LAYER: "<<tokens[tokens.size()-1].substr(3,1)<<std::endl;
				vertex_t w = boost::add_vertex(g);
				g[w].x = x2; g[w].y = y2; g[w].layer = atoi((tokens[tokens.size()-1].substr(3,1)).c_str()); g[w].is_via = true; g[w].is_via_sep = false;
				via.x = x2; via.y = y2; via.layer = atoi((tokens[tokens.size()-1].substr(3,1)).c_str()); via.is_via = true; via.is_via_sep = false;
	                        vialist[boost::lexical_cast<std::string>(x2)+" "+boost::lexical_cast<std::string>(y2)] = via;
			}
			//std::cout<<"X1: "<<x1<<",Y1: "<<y1<<"\tX2: "<<x2<<",Y2: "<<y2<<std::endl;
		}
		else if(tokens.size()==2){//if only via info
			boost::split(tokens, metal_lines[i].first, boost::is_any_of("("));
			std::stringstream ss(tokens[1]);
			ss>>word;
			x1 = atoi(word.c_str());
			ss>>word;
			y1 = atoi(word.c_str());
			boost::split(tokens, metal_lines[i].first, boost::is_any_of(" "));
            //add vertex for via
			vertex_t w = boost::add_vertex(g);
			g[w].x = x1; g[w].y = y1; g[w].layer = atoi((tokens[tokens.size()-1].substr(3,1)).c_str()); g[w].is_via = true; g[w].is_via_sep = true; g[w].via_line =  metal_lines[i].first;
			via.x = x1; via.y = y1; via.layer = atoi((tokens[tokens.size()-1].substr(3,1)).c_str());
			via.is_via = true; via.via_line = metal_lines[i].first; via.is_via_sep  = true;
//metal_lines[i].second;
			vialist[boost::lexical_cast<std::string>(x1)+" "+boost::lexical_cast<std::string>(y1)] = via;
			//std::cout<<"VIA LAYER: "<<tokens[tokens.size()-1].substr(3,1)<<std::endl;
		}

	}//end for loop i

	unsigned v_size = num_vertices(g);

    typedef boost::polygon::polygon_data<int> polygon;
    typedef boost::polygon::polygon_traits<polygon>::point_type point;

    std::unordered_map<std::string, vertex_t> umap;
	for(unsigned j=0; j<pinpinlist.size(); j++){
	    std::cout<<"NAME: "<<pinpinlist[j].first<< " " << pinpinlist[j].second << std::endl;
		for(unsigned i=0; i<v_size; i++){
	        //associate IO pins with vertex
			if(pinpinlist[j].first.find("PIN")!=-1){
                //PIN temp;
                //if(pinlist.find(pinpinlist[j].first+" "+pinpinlist[j].second) != pinlist.end())
    			PIN	temp = pinlist[pinpinlist[j].first+" "+pinpinlist[j].second];
                //std::cout <<"PIN: "<<temp.name<<"X,Y:"<<temp.x<<","<<temp.y<<std::endl;
				vertex_t u = vertex(i, g);
				if(temp.x == g[u].x && temp.y == g[u].y && temp.layer == g[u].layer && temp.layer <= split){
                    //std::cout << "MATCH PIN"<<std::endl;
		            vertex_t p = boost::add_vertex(g);
					g[p].name = "PIN "+temp.name;
					g[p].is_pin = true;
					edge_t e; bool b;
					boost::tie(e,b) = boost::add_edge(u,p,g);
                    g[u].is_open = false;
                    g[p].is_open = false;
					g[e].valid = false;
					break;
				}

			}
			
            //associate COMP pins with vertex
            else{
                CELL c = complist[pinpinlist[j].first];
	            //std::cout<<"NAME: "<<c.name_key<<std::endl;
	  	        vertex_t u = vertex(i, g);
              	point p =  boost::polygon::construct<point>(g[u].x, g[u].y);
        	  	//std::cout<<"PIN: "<<g[u].name<<", POINT X,Y: "<<g[u].x<<","<<g[u].y<<std::endl;
	            //std::cout<<"CELL: "<<c.name_key<<" X,Y: "<<c.x1<<","<<c.y1<<" "<<c.x2<<","<<c.y2<<" "<<c.x3<<","<<c.y3<<" "<<c.x2<<","<<c.y3<<std::endl;


                //create polygon for the cell
      		    point pts[4] = {
        		boost::polygon::construct<point>(c.x1, c.y1),
        		boost::polygon::construct<point>(c.x2, c.y2),
        		boost::polygon::construct<point>(c.x2, c.y3),
        		boost::polygon::construct<point>(c.x3, c.y3),
        		};

        		polygon poly;
        		boost::polygon::set_points(poly, pts, pts+4);
                //(boost::polygon::contains(poly, p, true))? std::cout<<"INSIDE"<<std::endl:std::cout<<"OUTSIDE"<<std::endl;
                
                //if point connects to the cell
        		if(boost::polygon::contains(poly, p, true) && g[u].layer == 1){
                    //std::cout<<"INSIDE"<<std::endl; 
                    if(umap.find(pinpinlist[j].first+" "+pinpinlist[j].second)==umap.end()){//if new pin save it
                        //std::cout<<"PIN CREATED: "<<pinpinlist[j].first<<std::endl;
                        vertex_t t = boost::add_vertex(g);
                        g[t].name = pinpinlist[j].first+" "+pinpinlist[j].second;
                        g[t].is_pin = true;
                        edge_t e; bool b;
                        boost::tie(e,b) = boost::add_edge(u,t,g);
                        g[u].is_open = false;
                        g[t].is_open = false;
                        g[e].valid = false;
                        umap[pinpinlist[j].first+" "+pinpinlist[j].second] = t;
                    }
                    else{//existing pin; then just connect
                        //std::cout<<"COMMON PIN"<<std::endl;
                        vertex_t t = umap[pinpinlist[j].first+" "+pinpinlist[j].second];
                        edge_t e; bool b;
                        boost::tie(e,b) = boost::add_edge(u,t,g);
                        g[u].is_open = false;
                        g[t].is_open = false;
                        g[e].valid = false;
                    }
                }
            }
		}//end for i
	}//end for j



	 //connect the graph
     for(unsigned i=0; i<v_size; i++){
		for(unsigned j=0; j<v_size && j!=i; j++){
			vertex_t u = vertex(i, g);
			vertex_t v = vertex(j, g);
			if(g[u].x==g[v].x && g[u].y==g[v].y){
				if(g[u].layer==g[v].layer){//same metal layer connection 
                    edge_t e; bool b;
                    boost::tie(e,b) = boost::add_edge(u,v,g);
                    g[u].is_open = false;
                    g[v].is_open = false;
                    g[e].x1 = g[u].x; g[e].y1 = g[u].y; g[e].x2 = g[v].x; g[e].y2 = g[v].y; g[e].layer = metal_lines[i].second; g[e].valid = false; g[e].via =false;
				}
				else if((g[u].layer == g[v].layer-1 && (g[u].is_via || g[v].is_via))){//metal layers connected by via
                    edge_t e; bool b;
                    boost::tie(e,b) = boost::add_edge(u,v,g);
                    g[u].is_open = false;
                    g[v].is_open = false;
                    g[e].x1 = g[u].x; g[e].y1 = g[u].y; g[e].x2 = g[v].x; g[e].y2 = g[v].y; g[e].layer = metal_lines[i].second; g[e].valid = false;  g[e].via =false;
				}

				else if((g[u].layer == g[v].layer+1 && (g[u].is_via || g[v].is_via))){//metal layers connected by via
                    edge_t e; bool b;
                    boost::tie(e,b) = boost::add_edge(u,v,g);
                    g[u].is_open = false;
                    g[v].is_open = false;
                    g[e].x1 = g[u].x; g[e].y1 = g[u].y; g[e].x2 = g[v].x; g[e].y2 = g[v].y; g[e].layer = metal_lines[i].second; g[e].valid = false;  g[e].via =false;
				}
			}
				
	
		}
	}

     
    //print current net

    //lock the file before writting
    omp_set_lock(&writelock);
	unsigned count = 0;
	std::vector<int> component_int (boost::num_vertices (g));
	size_t num_components = boost::connected_components (g, &component_int[0]);
    //cut_nets++;
	//std::cout<<"NET NAME: "<<lastline.first<<std::endl;
    cut_nets += boost::connected_components (g, &component_int[0]);
	//std::cout<<"#COMPONENTS: "<<num_components<< "#VERTICES: "<<boost::num_vertices (g)<<std::endl;
	for (auto const& component : connected_components_subgraphs(g))
    {
    bool contain_io=false;
    //print the pin first
	for(auto i = 0; i < boost::num_vertices(g); i++){
	  if (component_int[i] == count)
	    if(g[i].is_pin == true && g[i].name.find("PIN")!=-1){
            contain_io = true;
        }
        

    }
    if(contain_io){//if net contains IO pin, then keep the same net name
        //std::cout << lastline.first<< std::endl;
        out << lastline.first<< std::endl;

    }

    else{//otherwise, new net name
        //std::cout << lastline.first+"_"+std::to_string(count) << std::endl;
        out << lastline.first+"_split_"+std::to_string(count) << std::endl;
    }
    
    //print the PIN name
	for(auto i = 0; i < boost::num_vertices(g); i++){
	  if (component_int[i] == count)
	    if(g[i].is_pin == true){
	    	//std::cout <<"( "<< g[i].name <<" )"<< std::endl;
	    	out <<"( "<< g[i].name <<" )"<< std::endl;
        }
	}

    //print the metal lines
    bool init = true;
    for (auto e :  make_iterator_range(edges(component))){
        if(g[e].valid){
            //std::cout<<g[e].metal<<std::endl;
            if(init){//first line should be ROUTED
                   std::string temp = g[e].metal;
                   if(temp.find("ROUTED")!=-1){
                       out<<g[e].metal<<std::endl;
                       //out << "HABIJABI" << std::endl;
                       //std::cout<<g[e].metal<<std::endl;
                   }
                   else{
                        auto pos = temp.find("NEW");
                        temp = "\t+ ROUTED "+temp.substr(pos+3, temp.size()-pos-3);
                        //std::cout<<temp<<std::endl;
                        out<<temp<<std::endl;
                        //std::cout<<temp<<std::endl;
                    
                   }
                init = false;
            }
            else{    //std::cout<<g[e].metal<<std::endl;
                     out<<g[e].metal<<std::endl;
                     //std::cout<<g[e].metal<<std::endl;
            }
        }
    }


    //print the separate via lines    
	for(auto i = 0; i < boost::num_vertices(g); i++){
	  if (component_int[i] == count)
	    if(g[i].is_via_sep){
            if(init){
                std::string temp = g[i].via_line;
                if(temp.find("ROUTED")!=-1){
                    //std::cout<<temp<<std::endl;
                    out<<temp<<std::endl;
                }
                else{
                        auto pos = temp.find("NEW");
                        temp = "\t+ ROUTED "+temp.substr(pos+3, temp.size()-pos-3);
            		    //std::cout<<temp<<std::endl;
            		    out<<temp<<std::endl;
                }
                init = false;
            }
            else{
                //std::cout << g[i].via_line << std::endl;
                out << g[i].via_line << std::endl;
            }
        }
	}

    //print the last lines; FREQ, PROP etc.
    for(auto i=0; i<lastline.second.size(); i++){
           //std::cout<<lastline.second[i]<<std::endl;
           out<<lastline.second[i]<<std::endl;
    }
    //std::cout<<";"<<std::endl;
    out<<";"<<std::endl;
	//std::cout <<std::endl;
	out <<std::endl;

	count++;

    }//end for component
    //unlock the output file
	
    omp_unset_lock(&writelock);
	count = 0;
    
    //for (auto const& component : connected_components_subgraphs(g))
    //{
    //    //compute center
    //    std::vector<double> xx;
    //    std::vector<double> yy;

    //    for(auto i = 0; i < boost::num_vertices(g); i++){
    //      if (component_int[i] == count && g[i].layer == split){
    //        xx.push_back(g[i].x); 
    //        yy.push_back(g[i].y); 
    //      }
    //    }
    //    std::sort(xx.begin(), xx.end());
    //    std::sort(yy.begin(), yy.end());

    //    //std::cout<<"XX size: "<<xx.size()<<std::endl;
    //    if(xx.size()==0)
    //        continue;

    //    CONN cur_node;
    //    cur_node.x = (xx[0] + xx[xx.size()-1])/2;
    //    cur_node.y = (yy[0] + yy[yy.size()-1])/2;
    //    cur_node.p_net = lastline.first;
    //    if(cur_node.p_net.find("n_24")!=-1){
    //    std::cout << "BOUNDING BOX x1,y1:"<<xx[0]<<","<<yy[0]<<"\t, x2,y2:"<<xx[xx.size()-1]<<","<<yy[yy.size()-1]<<std::endl;
    //    std::cout<<"#COMP:"<<num_components<<std::endl;
    //    }
    //    xx.clear();
    //    yy.clear();
    //    
    //    for(auto i = 0; i < boost::num_vertices(g); i++){
    //        xx.push_back(g[i].x); 
    //        yy.push_back(g[i].y); 
    //    }
    //    
    //    if(xx.size()==0)
    //        continue;
    //    
    //    cur_node.r = sqrt(pow((xx[0] - xx[xx.size()-1]), 2) + pow((yy[0] + yy[yy.size()-1]),2))/2;
    //    omp_set_lock(&writelock);
    //    center.push_back(cur_node);
    //    omp_unset_lock(&writelock);
    //    count++;

    //}
    
    //for(int i=0; i<v_size; i++){
    //    if(g[i].is_open){
    //        vertex_t u = vertex(i, g);
    //        std::pair<adjacency_it, adjacency_it> neighbor = boost::adjacent_vertices(u, g);
    //        //std::cout<<"OPEN END"<<std::endl;
    //        if(g[u].x < g[*neighbor.first].x);

    //        else if(g[u].x > g[*neighbor.first].x);
    //        else if(g[u].y < g[*neighbor.first].y);
    //        else if(g[u].y < g[*neighbor.first].y);
    //    }
    //}

}

void print_pin(mypair cur_net, mypair lastline){

    omp_set_lock(&writelock);
	for(unsigned i=0; i<cur_net.second.size(); i++){
		//std::cout<<cur_net.second[i]<<std::endl;
		out<<cur_net.second[i]<<std::endl;
    }

    for(auto i=0; i<lastline.second.size(); i++){
           //std::cout<<lastline.second[i]<<std::endl;
           out<<lastline.second[i]<<std::endl;
    }
    //std::cout<<";"<<std::endl;
    out<<";\n"<<std::endl;
    omp_unset_lock(&writelock);
	
}

void parse_till_end(std::ifstream& def){//parse rest lines

    std::string line;
    while(getline(def, line)){
        std::cout<<line<<std::endl;
        out<<line<<std::endl;
    }
}

void compute_center(){

    //std::cout << "Inside compute center " <<std::endl;
    std::sort(center.begin(), center.end());

    //#pragma omp parallel for schedule(static)
    for(int i=0; i<center.size(); i++){
        double x = center[i].x, y = center[i].y, r = center[i].r;
        std::string p_net = center[i].p_net;
        //check the left side
        //if(p_net.find("n_24")!=-1){
        for(int j=i-1; j>=0; j--){
            double x2 = center[j].x, y2 = center[j].y;
            std::string p_net2 = center[j].p_net;
            double dist = sqrt(pow(x-x2,2) + pow(y-y2,2));
            if(dist > r)
                break;
            else{
                omp_set_lock(&writelock);
                std::cout << abs(x-x2) + abs(y-y2) <<"\t"<< (p_net == p_net2)<<std::endl;
                omp_unset_lock(&writelock);
            }

        }
        //check the right side
        for(int j=i+1; j<center.size(); j++){
            double x2 = center[j].x, y2 = center[j].y;
            std::string p_net2 = center[j].p_net;
            double dist = sqrt(pow(x-x2,2) + pow(y-y2,2));
            if(dist > r)
                break;
            else{
                omp_set_lock(&writelock);
                std::cout << abs(x-x2) + abs(y-y2) <<"\t"<< (p_net == p_net2)<<std::endl;
                omp_unset_lock(&writelock);
            }

        }
        //}
    }
}
