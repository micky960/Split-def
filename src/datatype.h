struct rect{
	std::vector<std::pair<int, int> > loc;
};

typedef struct rect RECT;
typedef std::pair<double, double> B;
struct cell{
	std::string name_key;
	B b;
	double x,y;
};

//Define a class that has the data you want to associate to every vertex and edge
struct node{
	std::string name;
	int layer;
	int x,y;
	bool is_pin;
};

struct edge{
	bool valid;
	int layer;
	double x1,y1, x2, y2;
	std::string metal;
};

typedef struct node NODE;
typedef struct edge EDGE;

typedef struct cell CELL;


std::unordered_map <std::string, NODE> nodelist, vialist, pinlist;

std::unordered_map <std::string, CELL> complist;

typedef std::pair<std::string, std::vector<std::string> > mypair;
std::vector<mypair> netlist;
int split;
std::vector<std::string> keep_layer, track_lines;
typedef std::pair <std::string, std::string> pin_pair;
typedef std::pair <std::string, int> metal_pair;


