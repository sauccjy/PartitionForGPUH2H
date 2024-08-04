#pragma once
#include"Graph_D.cuh"

Graph_D_H::Graph::Graph(const string graphName ,const string Nodefile, const string EdgeFile, int threadPoolSize, int changeHeight)
{
	this->graphName = graphName;
	this->NodeFile = Nodefile;
	this->EdgeFile = EdgeFile;
	this->changeHeight = changeHeight;
	time_Mine time;
	cout << "graph node start loading" << endl;
	time.updateStart();
	fstream Node(Nodefile, ios::in | ios::out);
	int tempint = 0;
	Node >> NodeNumber;
	NE_P.assign(NodeNumber,myPair<double>());
	for (int i = 0; i < NodeNumber; i++)
	{
		char s;
		int ID;
		double N, E;
		Node >> s >> ID >> N >> E;
		NE_P[ID - 1].setPairs(ID - 1,N,E);
	}
	Node.close();
	time.updateEnd();
	cout << "\t graph node end and using time:" << time.get_microsecond_duration() << endl;

	cout << "graph edge start loading" << endl;
	time.updateStart();
	fstream Edge(EdgeFile, ios::in | ios::out);
	Edge >>NodeNumber>> EdgeNumber;
	int NodeID1 = 0, NodeID2 = 0, EdgeWeight = 0;
	char s;
	//int tempNodeID = 0;
	int actualEdgeNumber = EdgeNumber;
	adjList.assign( NodeNumber , thrust::host_vector<pairs>() );

	for (int i = 0; i < EdgeNumber; i++)
	{
		Edge >>s >> NodeID1 >> NodeID2 >> EdgeWeight;
		if (NodeID1 == NodeID2)
		{
			actualEdgeNumber--;
			continue;
		}
		adjList[NodeID1 - 1].push_back(pairs(NodeID2 - 1, EdgeWeight));
	}
	EdgeNumber = actualEdgeNumber;
	Edge.close();
	time.updateEnd();
	cout << "\t graph edge finished and using time:" << time.get_microsecond_duration() << endl;
}

Graph_D_H::Graph::~Graph()
{


}

void Graph_D_H::Graph::checklink()
{
	vector<bool> visit(NodeNumber,false);
	priority_queue<int> que;
	int size = 0;
	que.push(0);
	visit[0] = true;
	while (!que.empty()) {
		int ID = que.top();
		que.pop();
		size++;
		for (int i = 0; i < adjList[ID].size(); i++)
		{
			if (visit[adjList[ID][i].first])
				continue;
			que.push(adjList[ID][i].first);
			visit[adjList[ID][i].first] = true;
		}

	}
	if (size == NodeNumber) {
		cout << "connect!"<<endl;
	}
	else {
		cout << "disconnect!"<<endl;
	}
}

