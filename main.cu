
#include"Graph_D.cuh"
using namespace std;

#define UseLatitude
const double beta = 0.4;


void Partition_1(Graph_D_H::Graph* Agraph, int goalHeight)
{
	//Agraph->displayAdjList();
	Graph_D_H::time_Mine time;
	//Agraph->startGPU();
	cout << "partition start " << endl;
	time.updateStart();
	int actualSize = Agraph->partition(goalHeight);
	time.updateEnd();
	Agraph->partitionTime = time.get_microsecond_duration();
	cout << "\t parition using time: " << time.get_microsecond_duration() << endl;
	//Agraph->displayPartition();
	Agraph->writePartition();
	Agraph->makeAdjcentNode_2(0);
	
}         
void Partition_2(Graph_D_H::Graph* Agraph, int goalHeight)
{
	//Agraph->displayAdjList();
	Graph_D_H::time_Mine time;
	//Agraph->startGPU();
	cout << "partition start " << endl;
	time.updateStart();
	int actualSize = Agraph->partition_minumum_cut(goalHeight, beta);
	time.updateEnd();
	Agraph->partitionTime = time.get_microsecond_duration();
	cout << "\t parition using time: " << time.get_microsecond_duration() << endl;

	//Agraph->displayPartition();
	Agraph->writePartition();
	Agraph->makeAdjcentNode_2(0);

}

int main(int argc, char** argv)
{

	string NodeFile, EdgeFile;
	NodeFile = "./Graph/USA/";
	EdgeFile = "./Graph/USA/";
	if (argc != 4)
		return;
	string gname(argv[1]);
	char* endptr;

	int secondArgLong = (int)strtol(argv[2], &endptr, 10);//tree height
	int thirdArgLong = (int)strtol(argv[3], &endptr, 10); //partition method == 1 means latitude == 2 means minimum cut

	int threadPoolSize = 1; //useless
	int TreeHeight = secondArgLong;
	int method = thirdArgLong;
	string graphName = gname;
	NodeFile += graphName;
	EdgeFile += graphName;
	NodeFile += "Node.txt";
	EdgeFile += "Edge.txt";
	cout << NodeFile << endl;
	cout << EdgeFile << endl;

	Graph_D_H::Graph* Agraph = new Graph_D_H::Graph(graphName, NodeFile, EdgeFile, threadPoolSize, TreeHeight);
	Agraph->checklink();
	if(method == 1)
	Partition_1(Agraph, TreeHeight);
	else {
		Partition_2(Agraph, TreeHeight);
	}
	//

}