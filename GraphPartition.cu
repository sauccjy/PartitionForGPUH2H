#pragma once
#include"Graph_D.cuh"

namespace Graph_D_H
{
	int Graph_D_H::Graph::partition(const int goalHeight)
	{
		partition_Tree.push_back(pairs(0, NodeNumber));
		ID_hash.assign(NodeNumber, 0);
		TreeHeight = 1;
		thrust::sort(NE_P.begin(), NE_P.end(), myPair_longitude_less_than());
		bool longitudeLast = true;

		int maxSize = NodeNumber;

		long long int lastPartitionTree = 0;
		long long int nowPartitionTree = 1;
		while (TreeHeight < goalHeight)
		{
			long long int tempMax = -1;
			for (long long int i = lastPartitionTree; i < nowPartitionTree; i++)
			{
				int left = partition_Tree[i].first;
				int right = partition_Tree[i].second;
				int mid = (right - left) / 2 + left;
				if (left == right) { //empty
					partition_Tree.push_back(pairs(left, left));
					partition_Tree.push_back(pairs(left, left));
					continue;
				}
				if (right - left == 1) {
					partition_Tree.push_back(pairs(left, right));
					partition_Tree.push_back(pairs(right, right));
					continue;
				}
				if (right - left == 2) {
					partition_Tree.push_back(pairs(left, right - 1));
					partition_Tree.push_back(pairs(right - 1, right));
					continue;
				}
				//left -> mid; mid+1 ->  right
				if (longitudeLast)
				{
					thrust::sort(NE_P.begin() + left, NE_P.begin() + mid, myPair_latitude_less_than());
					thrust::sort(NE_P.begin() + mid, NE_P.begin() + right, myPair_latitude_less_than());
				}
				else
				{
					thrust::sort(NE_P.begin() + left, NE_P.begin() + mid, myPair_longitude_less_than());
					thrust::sort(NE_P.begin() + mid, NE_P.begin() + right, myPair_longitude_less_than());
				}

				partition_Tree.push_back(pairs(left, mid));
				partition_Tree.push_back(pairs(mid, right));

				tempMax = ((right - mid) > tempMax) ? (right - mid) : tempMax;
			}
			TreeHeight++;
			longitudeLast = !longitudeLast;
			lastPartitionTree = nowPartitionTree;
			nowPartitionTree = partition_Tree.size();
			maxSize = tempMax;
		}

		thrust::host_vector<int> Degree(NodeNumber);
		for (int i = 0; i < NodeNumber; i++)
		{
			Degree[i] = adjList[i].size();
		}

		int lowestIndexStart = std::pow(2, TreeHeight - 1) - 1;
		int lowestIndexEnd = std::pow(2, TreeHeight) - 2;
		for (int i = lowestIndexEnd; i >= lowestIndexStart; i--)
		{
			thrust::sort_by_key(thrust::host, Degree.begin() + partition_Tree[i].first, Degree.begin() + partition_Tree[i].second,
				NE_P.begin() + partition_Tree[i].first);
		}
		for (int i = 0; i < NodeNumber; i++)
		{
			ID_hash[NE_P[i].NodeID] = i;
		}
		lowestMaxSize = maxSize;
		return maxSize;
	}

	int Graph::partition_Latitude_first(int goalHeight)
	{
		partition_Tree.push_back(pairs(0, NodeNumber));
		ID_hash.assign(NodeNumber, 0);
		TreeHeight = 1;
		thrust::sort(NE_P.begin(), NE_P.end(), myPair_latitude_less_than());
		bool longitudeLast = true;

		//int maxSize = NE_P.size() / 2 + (int)(NE_P.size() % 2 != 0);
		int maxSize = NodeNumber;

		int lastPartitionTree = 0;
		int nowPartitionTree = 1;
		while (TreeHeight < goalHeight)
		{
			int tempMax = -1;
			for (int i = lastPartitionTree; i < nowPartitionTree; i++)
			{
				int left = partition_Tree[i].first;
				int right = partition_Tree[i].second;
				int mid = (right - left) / 2 + left;
				if (left == right) { //empty
					partition_Tree.push_back(pairs(left, left));
					partition_Tree.push_back(pairs(left, left));
					continue;
				}
				if (right - left == 1) {
					partition_Tree.push_back(pairs(left, right));
					partition_Tree.push_back(pairs(right, right));
					continue;
				}
				if (right - left == 2) {
					partition_Tree.push_back(pairs(left, right - 1));
					partition_Tree.push_back(pairs(right - 1, right));
					continue;
				}
				//left -> mid; mid+1 ->  right
				if (!longitudeLast)
				{
					thrust::sort(NE_P.begin() + left, NE_P.begin() + mid, myPair_latitude_less_than());
					thrust::sort(NE_P.begin() + mid, NE_P.begin() + right, myPair_latitude_less_than());
				}
				else
				{
					thrust::sort(NE_P.begin() + left, NE_P.begin() + mid, myPair_longitude_less_than());
					thrust::sort(NE_P.begin() + mid, NE_P.begin() + right, myPair_longitude_less_than());
				}

				partition_Tree.push_back(pairs(left, mid));
				partition_Tree.push_back(pairs(mid, right));

				tempMax = ((right - mid) > tempMax) ? (right - mid) : tempMax;
			}
			TreeHeight++;
			longitudeLast = !longitudeLast;
			lastPartitionTree = nowPartitionTree;
			nowPartitionTree = partition_Tree.size();
			maxSize = tempMax;
		}

		//thrust::host_vector<int> Degree(NodeNumber);
		//for (int i = 0; i < NodeNumber; i++)
		//{
		//	Degree[i] = CSR_node_OutdegreePoint[NE_P[i].NodeID + 1] - CSR_node_OutdegreePoint[NE_P[i].NodeID];
		//}

		//int lowestIndexStart = std::pow(2, TreeHeight - 1) - 1;
		//int lowestIndexEnd = std::pow(2, TreeHeight) - 2;
		//for (int i = lowestIndexEnd; i >= lowestIndexStart; i--)
		//{
		//	thrust::sort_by_key(thrust::host, Degree.begin() + partition_Tree[i].first, Degree.begin() + partition_Tree[i].second,
		//		NE_P.begin() + partition_Tree[i].first);
		//}
		for (int i = 0; i < NodeNumber; i++)
		{
			ID_hash[NE_P[i].NodeID] = i;
		}
		lowestMaxSize = maxSize;
		return maxSize;
	}



	int  Graph_D_H::Graph::generateBFS(std::map<int, std::map<int, int>>& graph, std::map<int, bool>& visited_s, std::map<int, bool>& visited_t, int s, double beta) {

		if (graph.find(s) == graph.end()) {
			cout << s << " is not in graph" << endl;
			return -1;
		}
		//visited_s[s] = true;
		if (graph.size() < 2)
		{
			cout << s << " is in graph but graph size < 2" << endl;
			return 0;
		}

		//std::cout << "\t In partition graph: (";
		//for (auto& it : graph) {
		//	std::cout << it.first << ",";
		//}
		//cout << ")" << endl;;

		int target_size = (int)(beta * (double)(graph.size()));
		//cout << "target_size: " << target_size << endl;
		std::queue<int> q;
		q.push(s);
		//visited_s[s] = true;
		int markSize = 0;
		//cout << "\t BFS chosen ID: (";
		//cout << s;

		while (!q.empty()) {
			int nodeNow = q.front();
			q.pop();
			if (visited_t[nodeNow]|| visited_s[nodeNow]) {
				continue;
			}
			visited_s[nodeNow] = true;
			//cout << nodeNow << ",";
			markSize++;
			if (markSize > target_size)
				break;
			for (auto& it : graph[nodeNow]) {
				if (!visited_s[it.first]) {
					q.push(it.first);
				}
			}
		}

		//while (!q.empty() && bfs_result.size() < target_size) {
		//	int level_size = q.size();
		//	std::vector<int> current_level;

		//	for (int i = 0; i < level_size; ++i) {
		//		int u = q.front();
		//		q.pop();
		//		bfs_result.push_back(u);
		//		current_level.push_back(u);

		//		for (const auto& v : graph[u]) {
		//			if (!visited_s[v.first]) {
		//				visited_s[v.first] = true;
		//				q.push(v.first);
		//			}
		//		}
		//	}

		//	// Sort current level by degree (number of neighbors)
		//	std::sort(current_level.begin(), current_level.end(), [&graph](int a, int b) {
		//		return graph[a].size() < graph[b].size();
		//		});

		//	// Update BFS result with sorted current level
		//	for (const int& node : current_level) {
		//		if (bfs_result.size() >= target_size) break;
		//		bfs_result.push_back(node);
		//	}
		//}

		//// Mark the nodes in bfs_result as visited
		//for (const int& node : bfs_result) {
		//	visited_s[node] = true;
		//	cout << node<<",";
		//}
		//cout << ")\n";

		return markSize;
	}


	void Graph_D_H::Graph::generateFarthestVertex(std::map<int, std::map<int, int>>& standard_graph, int s, int& f) {
		queue<pair<int, int>> run_que;
		run_que.push(make_pair(s, 1));
		vector<bool> visited_s(NodeNumber, false);
		int tempFarthest = s;
		//visited_s[s] = true;
		int tempMaxHeight = 1;
		//cout << "searching node: " << s << "'th farthest node" << endl;
		while (!run_que.empty()) {
			auto& it = run_que.front();
			int tempID = it.first;
			int tempHeight = it.second;
			run_que.pop();

			if (visited_s[tempID]) continue;

			visited_s[tempID] = true;
			if (tempHeight > tempMaxHeight) {
				tempFarthest = tempID;
				tempMaxHeight = tempHeight;
			}
			else {
				if (standard_graph[tempID].size() < standard_graph[tempFarthest].size() ) {
					tempFarthest = tempID;
					tempMaxHeight = tempHeight;
				}
			}
			//cout << "height: " << tempHeight << " is: " << tempID << endl;
			for (auto& it : standard_graph[tempID]) {
				if (!visited_s[it.first]) {
					run_que.push(make_pair(it.first, tempHeight + 1));
					
				}
					//run_que.push(make_pair(it.first, tempHeight + 1));
			}

		}
		f = tempFarthest;
	}


	bool Graph_D_H::Graph::bfs_for_edmondsKarp(const std::map<int, std::map<int, int>>& residualGraph, int source, int sink, std::map<int, int>& parent) {
		std::map<int, bool> visited;
		std::queue<int> q;
		q.push(source);
		visited[source] = true;
		parent[source] = -1;

		while (!q.empty()) {
			int u = q.front();
			q.pop();

			if (residualGraph.find(u) == residualGraph.end()) continue; // Check if u exists in residualGraph

			for (const auto& v : residualGraph.at(u)) {
				if (!visited[v.first] && v.second > 0) {
					q.push(v.first);
					parent[v.first] = u;
					visited[v.first] = true;
					if (v.first == sink)
						return true;
				}
			}
		}

		return false;
	}

	// Function to implement the Edmonds-Karp algorithm
	int  Graph_D_H::Graph::edmondsKarp(std::map<int, std::map<int, int>>& graph, int source, int sink, std::map<int, std::map<int, int>>& residualGraph) {
		residualGraph = graph; // Initialize residual graph
		std::map<int, int> parent;
		int maxFlow = 0;

		// Augment the flow while there is a path from source to sink
		while (bfs_for_edmondsKarp(residualGraph, source, sink, parent)) {
			int pathFlow = INT_MAX;

			// Find the maximum flow through the path found by BFS
			for (int v = sink; v != source; v = parent[v]) {
				int u = parent[v];
				pathFlow = std::min(pathFlow, residualGraph[u][v]);
			}

			// Update residual capacities of the edges and reverse edges along the path
			for (int v = sink; v != source; v = parent[v]) {
				int u = parent[v];
				residualGraph[u][v] -= pathFlow;
				residualGraph[v][u] += pathFlow;
			}

			// Add path flow to the overall flow
			maxFlow += pathFlow;
		}

		return maxFlow;
	}

	// Function to find the minimum cut using the residual graph
	void  Graph_D_H::Graph::findMinCut(const std::map<int, std::map<int, int>>& residualGraph, int source, std::map<int, bool>& visited) {
		std::queue<int> q;
		q.push(source);
		visited[source] = true;

		while (!q.empty()) {
			int u = q.front();
			q.pop();

			if (residualGraph.find(u) == residualGraph.end()) continue; // Check if u exists in residualGraph

			for (const auto& v : residualGraph.at(u)) {
				if (!visited[v.first] && v.second > 0) {
					q.push(v.first);
					visited[v.first] = true;
				}
			}
		}
	}

	// Main function to find the minimum cut and split the graph
	void  Graph_D_H::Graph::findMinCutAndSplitGraph(std::map<int, std::map<int, int>>& graph, int source, int sink, 
		std::set<int>& source_paritition, std::set<int>& sink_paritition) {
		// Convert undirected graph to directed graph with unit capacity

		std::map<int, std::map<int, int>> residualGraph;
		int maxFlow = edmondsKarp(graph, source, sink, residualGraph);
		//std::cout << "Maximum flow: " << maxFlow << std::endl;

		std::map<int, bool> visited;
		for (const auto& u : residualGraph) {
			visited[u.first] = false;
		}
		findMinCut(residualGraph, source, visited);

		for (const auto& u : graph) {
			if (visited.at(u.first)) {
				source_paritition.insert(u.first);
			}
			else {
				sink_paritition.insert(u.first);
			}
		}
	}






	int Graph_D_H::Graph::partition_minumum_cut(int goalHeight, double beta) {
		PartitionMethod = 2;
		this->Beta = beta;
		partition_Tree.push_back(pairs(0, NodeNumber));
		ID_hash.assign(NodeNumber, 0);
		TreeHeight = 1;
		thrust::host_vector<myPair<double> > NE_P_copy = NE_P;
		int maxSize = NodeNumber;

		int lastPartitionTree = 0;
		int nowPartitionTree = 1;
		while (TreeHeight < goalHeight)
		{
			int tempMax = -1;

			cout << "_________________at layer:" << TreeHeight << "____________" << endl;
			for (int i = lastPartitionTree; i < nowPartitionTree; i++)
			{
				int left = partition_Tree[i].first;
				int right = partition_Tree[i].second;
				//cout << "in partition: <" << left << "." << right << "," << endl;
				if (left == right) { //empty
					partition_Tree.push_back(pairs(left,left));
					partition_Tree.push_back(pairs(left, left));
					continue;
				}
				if (right - left == 1) {
					partition_Tree.push_back(pairs(left, right));
					partition_Tree.push_back(pairs(right, right));
					continue;
				}
				if (right - left == 2) {
					partition_Tree.push_back(pairs(left, right - 1));
					partition_Tree.push_back(pairs(right - 1, right));
					continue;
				}

				//create copy and delete edges which not in this partition
				std::map<int, std::map<int, int>> standard_graph_in_partition;
				map<int, bool> firstPartition, secondPartition;
				int randomNode = -1;
				for (int j = left; j < right; j++) {

					int NodeID = NE_P[j].NodeID;
					randomNode = NodeID;

					standard_graph_in_partition.emplace(NodeID, std::map<int, int>());
					firstPartition.emplace(NodeID, false);
					secondPartition.emplace(NodeID, false);
				}
				
				for (int j = left; j < right; j++) {

					int NodeID = NE_P[j].NodeID;
					randomNode = (adjList[NodeID].size() > adjList[randomNode].size()) ? NodeID : randomNode;
					for (auto& it : adjList[NodeID]) {
						if (standard_graph_in_partition.find(it.first) != standard_graph_in_partition.end())
							standard_graph_in_partition[NodeID].emplace(it.first, 1);
					}
				}
				//displayGraph(standard_graph_in_partition);
				//generate two farthest vertex
				
				int s = -1, t = -1;
				generateFarthestVertex(standard_graph_in_partition, randomNode, s);
				generateFarthestVertex(standard_graph_in_partition, s, t);
				secondPartition.at(t) = true;
				generateBFS(standard_graph_in_partition, firstPartition, secondPartition, s, beta);
				secondPartition.at(t) = false;
				generateBFS(standard_graph_in_partition, secondPartition,firstPartition, t, beta);
				//std::cout << "choose id: " << s << " and " << t << endl;

				//delete such patititon( if exists v in both S and T, put them in mid partition , if not in S and T , push them too)
				standard_graph_in_partition.emplace(NodeNumber + 1, std::map<int, int>());//S
				standard_graph_in_partition.emplace(NodeNumber + 2, std::map<int, int>());//T
				std::set<int> s_partition = {}, t_partition = {}, mid_partiiton = {}, merge_partition = {};
				for (auto& io : firstPartition) {
					int ID = io.first;
					if (io.second && !secondPartition[ID]) {
						s_partition.insert(ID);
					}
					else if (!io.second && secondPartition[ID]) {
						t_partition.insert(ID);
					}
					//else if (!io.second && !secondPartition[ID]) {
					//	mid_partiiton.insert(ID);
					//}
					else {
						mid_partiiton.insert(ID);
						//merge_partition.insert(ID);
					}
				}
				//if (TreeHeight == 1) {
				//	cout << "s size: " << s_partition.size() << " t size: " << t_partition.size() <<" mid size: "<<mid_partiiton.size()
				//		<<" merge size: "<< merge_partition.size()<< endl;
				//}
				//if (s_partition.size() > t_partition.size()) {
				//	t_partition.insert(merge_partition.begin(), merge_partition.end());
				//}
				//else {
				//	s_partition.insert(merge_partition.begin(), merge_partition.end());
				//}

				//std::cout << "vertex in s_partition: ";
				//for (auto& it : s_partition) {
				//	cout << it << ",";
				//}
				//cout << endl;
				//std::cout << "vertex in t_partition: ";
				//for (auto& it : t_partition) {
				//	cout << it << ",";
				//}
				//cout << endl;
				//std::cout << "vertex in mid_partition: ";
				//for (auto& it : mid_partiiton) {
				//	cout << it << ",";
				//}
				//cout << endl;

				queue<int> answer_queue = {};
				map<int, bool> is_visited = {};
				for (auto& io : firstPartition) {
					is_visited.emplace(io.first, false);
				}
				answer_queue.push(s);
				while (!answer_queue.empty()) {
					int nodenow = answer_queue.front();
					answer_queue.pop();
					if (is_visited.find(nodenow)->second) {
						continue;
					}
					is_visited[nodenow] = true;
					for (auto& it : standard_graph_in_partition[nodenow]) {
						if (s_partition.find(it.first) != s_partition.end()) {
							if (!is_visited.find(it.first)->second)
								answer_queue.push(it.first);
							continue;
						}
						else {
							standard_graph_in_partition[it.first].emplace(NodeNumber + 1, 1);
							standard_graph_in_partition[it.first].erase(nodenow);
							standard_graph_in_partition[NodeNumber + 1].emplace(it.first, 1);
						}
					}
				}

				for (auto& its : s_partition) {
					standard_graph_in_partition[its].clear();
					standard_graph_in_partition.erase(its);
				}

				is_visited.clear();
				for (auto& io : firstPartition) {
					is_visited.emplace(io.first, false);
				}
				answer_queue.push(t);
				while (!answer_queue.empty()) {
					int nodenow = answer_queue.front();
					answer_queue.pop();
					if (is_visited.find(nodenow)->second) {
						continue;
					}
					is_visited[nodenow] = true;
					for (auto& it : standard_graph_in_partition[nodenow]) {
						if (t_partition.find(it.first) != t_partition.end()) {
							if (!is_visited.find(it.first)->second)
								answer_queue.push(it.first);
							continue;
						}
						else {
							standard_graph_in_partition[it.first].emplace(NodeNumber + 2, 1);
							standard_graph_in_partition[it.first].erase(nodenow);
							standard_graph_in_partition[NodeNumber + 2].emplace(it.first, 1);
						}
					}
				}

				for (auto& its : t_partition) {
					standard_graph_in_partition[its].clear();
					standard_graph_in_partition.erase(its);
				}
				//displayGraph(standard_graph_in_partition);
				//find minimum cut
				std::set<int> source_paritition = {};
				std::set<int> sink_paritition = {};
				findMinCutAndSplitGraph(standard_graph_in_partition, NodeNumber + 1, NodeNumber + 2, source_paritition, sink_paritition);

				int sourceSizeGap = s_partition.size() + source_paritition.size() - (t_partition.size() + sink_paritition.size());
				
				std::set<int> source_paritition_2 = {};
				std::set<int> sink_paritition_2 = {};
				findMinCutAndSplitGraph(standard_graph_in_partition, NodeNumber + 2, NodeNumber + 1, source_paritition_2, sink_paritition_2);
				int sourceSizeGap2 = s_partition.size() + source_paritition_2.size() - (t_partition.size() + sink_paritition_2.size());

				if (sourceSizeGap > sourceSizeGap2) {
					source_paritition.clear();
					source_paritition.insert(source_paritition_2.begin(), source_paritition_2.end());
					source_paritition_2.clear();
					sink_paritition.clear();
					sink_paritition.insert(sink_paritition_2.begin(), sink_paritition_2.end());
					sink_paritition_2.clear();
				}

				if (source_paritition.find(NodeNumber + 1) != source_paritition.end()) {
					source_paritition.erase(NodeNumber + 1);
					source_paritition.insert(s_partition.begin(), s_partition.end());
					sink_paritition.erase(NodeNumber + 2);
					sink_paritition.insert(t_partition.begin(), t_partition.end());
				}
				else {
					source_paritition.erase(NodeNumber + 2);
					source_paritition.insert(t_partition.begin(), t_partition.end());
					sink_paritition.erase(NodeNumber + 1);
					sink_paritition.insert(s_partition.begin(), s_partition.end());
				}

				int mid = source_paritition.size() + left;

				vector<myPair<double> > NE_P_first = {};
				//thrust::host_vector<myPair<double> > NE_P_second = {};
				for (auto& it : source_paritition) {
					NE_P_first.push_back(myPair<double>(NE_P_copy[it].first, NE_P_copy[it].second, it));
				}
				for (auto& it : sink_paritition) {
					NE_P_first.push_back(myPair<double>(NE_P_copy[it].first, NE_P_copy[it].second, it));
				}
				if (NE_P_first.size() == right - left) {
					thrust::copy(NE_P_first.begin(), NE_P_first.end(), NE_P.begin() + left);
				}

				partition_Tree.push_back(pairs(left, mid));
				partition_Tree.push_back(pairs(mid, right));

				tempMax = ((right - mid) > tempMax) ? (right - mid) : tempMax;

			}
			TreeHeight++;



			lastPartitionTree = nowPartitionTree;
			nowPartitionTree = partition_Tree.size();
			maxSize = tempMax;
		}


		//thrust::host_vector<int> Degree(NodeNumber);
		//for (int i = 0; i < NodeNumber; i++)
		//{
		//	Degree[i] = CSR_node_OutdegreePoint[NE_P[i].NodeID + 1] - CSR_node_OutdegreePoint[NE_P[i].NodeID];
		//}

		//int lowestIndexStart = std::pow(2, TreeHeight - 1) - 1;
		//int lowestIndexEnd = std::pow(2, TreeHeight) - 2;
		//for (int i = lowestIndexEnd; i >= lowestIndexStart; i--)
		//{
		//	thrust::sort_by_key(thrust::host, Degree.begin() + partition_Tree[i].first, Degree.begin() + partition_Tree[i].second,
		//		NE_P.begin() + partition_Tree[i].first);
		//}
		for (int i = 0; i < NodeNumber; i++)
		{
			ID_hash[NE_P[i].NodeID] = i;
		}
		lowestMaxSize = maxSize;
		return maxSize;
	}



	void Graph::makeAdjcentNode_2(int goalHeight)
	{
		Graph_D_H::time_Mine time;

		cout << "start construct Partition Rank tree " << endl;
		time.updateStart();

		int tempHeight = TreeHeight;
		thrust::host_vector<bool> isAdjcent(NodeNumber, true);
		//thrust::host_vector<pairs> range(NodeNumber, pairs(INT_MAX, INT_MIN));
		thrust::host_vector<int> placeOccupy(NodeNumber, 0);

		for (int i = 0; i < NodeNumber; i++) {
			placeOccupy[i] = adjList[i].size();
		}

		NANHash.assign(partition_Tree.size(), myPair<int>());
		nonAdjcentNode.clear();
		vector<double> afs(TreeHeight + 1, 0);
		while (tempHeight > goalHeight)
		{
			int lowestIndexStart = std::pow(2, tempHeight - 1) - 1;
			int lowestIndexEnd = std::pow(2, tempHeight) - 2;
			thrust::host_vector<bool> isAdjcent_cp = isAdjcent;
			for (int i = lowestIndexEnd; i >= lowestIndexStart; i--)
			{
				NANHash[i].first = nonAdjcentNode.size();//index
				NANHash[i].second = 0; //maxSize
				NANHash[i].NodeID = 0;//nonAdjcent Size
				int left = partition_Tree[i].first;
				int right = partition_Tree[i].second;
				//vector<int> outPartitionSize(right - left, 0);
				map<int, int> subgraph;
				for (int j = left; j < right; j++)//construct sub graph
				{
					if (!isAdjcent[NE_P[j].NodeID])
						continue;
					int nodeID = NE_P[j].NodeID;
					subgraph.emplace(nodeID, 0);
				}
				NANHash[i].second = subgraph.size();//record max size
				for (auto& it : subgraph) {//define outsize
					it.second = subgraph.size();
					for (auto& adj : adjList[it.first]) {
						if (!isAdjcent[adj.first]) {
							continue;
						}
						if (subgraph.find(adj.first) == subgraph.end()) {
							it.second++;
							//isAdjcent_cp[it.first] = false;
						}
					}
				}

				for (auto& it : subgraph) {
					if (it.second == NANHash[i].second) {// mark non-adjNode
						isAdjcent_cp[it.first] = false;
						NANHash[i].NodeID++;
						nonAdjcentNode.push_back(TDrank(adjList[it.first].size(), it.first, 0));
					}
					placeOccupy[it.first] = max(placeOccupy[it.first], it.second);

				}
			}
			isAdjcent = isAdjcent_cp;
			afs[tempHeight] = nonAdjcentNode.size();
			tempHeight--;
		}

		time.updateEnd();
		//makePartitionRankTreeTime = time.get_microsecond_duration();
		cout << "\t construct Partition Rank tree end,  using time: " << time.get_microsecond_duration() << endl;
		string latitudePartition = "LatitudePartitionRate.csv";
		if (PartitionMethod == 2) latitudePartition = "MinimumPartitionRate.csv";
		std::fstream heightdata(latitudePartition, ios::in | ios::out | ios::app);

		heightdata << graphName << ",";
		cout << "NodeNumber per layer: \n";
		for (int i = 1; i < afs.size(); i++) {
			heightdata << afs[i] / NodeNumber << ",";
			cout << "\t layer : " << i << " size: " << afs[i] << endl;
		}
		heightdata << "\n";
		heightdata.close();


		Reflash_for_PRT.assign(nonAdjcentNode.size(), -1);
		for (auto i = 0; i < nonAdjcentNode.size(); i++) {
			int ID = nonAdjcentNode[i].second;
			Reflash_for_PRT[ID] = i;
		}
		//Reflash_for_PRT_D = Reflash_for_PRT;

		//allowcate CHTree
		cout << "start Allocate NUB " << endl;
		time.updateStart();
		CHTree.clear();
		CHTreeHash.assign(NodeNumber, myPair<int>());

		for (int i = 0; i < NodeNumber; i++)
		{
			int nodeID = NE_P[i].NodeID;
			CHTreeHash[i].first = 0;//actual size
			CHTreeHash[i].second = placeOccupy[nodeID];//maxSize
			//cout << placeOccupy[i] << endl;
			CHTreeHash[i].NodeID = CHTree.size();//CHTree index
			CHTree.insert(CHTree.end(), CHTreeHash[i].second, myPair<int>(INT_MAX, 0, -1)); //insert maxSize label(outpoint,weight,hub)

			for (auto& it : adjList[nodeID])
			{
				CHTree[CHTreeHash[i].NodeID + CHTreeHash[i].first].setPairs(-1, it.first, it.second);
				CHTreeHash[i].first++;
				if (CHTreeHash[i].first > CHTreeHash[i].second) {
					cout << "NUB error~! ::: NodeID: " << nodeID << " " << CHTreeHash[i].first - CHTreeHash[i].second << endl;
					int nonAdjcentNodeIndex = Reflash_for_PRT[nodeID];
					int inHeight = 1;
					while (1) {
						if (afs[inHeight] < nonAdjcentNodeIndex) {
							break;
						}
						inHeight++;
						//afs[inHeight] > 
					}
					cout << "nonAdjcentNodeIndex : " << nonAdjcentNodeIndex << endl;
					cout << "\t actual adjcent size: " << adjList[nodeID].size() << " nub size£º" << CHTreeHash[i].second << " at height:" << inHeight - 1 << endl;
					//throw("LUB allocate Error!!");
				}
			}
		}

		chSize = CHTree.size();

		time.updateEnd();
		//AllocateLUBTime = time.get_microsecond_duration();
		cout << "\t Malloc LUB end,  using time: " << time.get_microsecond_duration() << endl;
		//cout << "LUB Size: " << ((double)(chSize) * sizeof(Graph_D_H::myPair<int>)) / (1024 * 1024) << "MB" << endl;
		//LUBSize = ((double)(chSize) * sizeof(Graph_D_H::myPair<int>)) / (1024 * 1024);
		cout << "LUB Size: " << ((double)(chSize) * sizeof(Graph_D_H::myPair<int>)) / (1024 * 1024) << "MB" << endl;

	}


	void Graph_D_H::Graph::displayPartition()
	{
		cout << "______________________PARTITION________________________" << endl;
		int tempHeight = 1;
		while (tempHeight <= TreeHeight)
		{
			cout << "height:" << tempHeight << ";  ";
			int lowestIndexStart = std::pow(2, tempHeight - 1) - 1;
			int lowestIndexEnd = std::pow(2, tempHeight) - 2;
			tempHeight++;
			for (int i = lowestIndexStart; i <= lowestIndexEnd; i++)
			{
				int left = partition_Tree[i].first;
				int right = partition_Tree[i].second;
				if (left == right) {
					cout << "()" ;
					continue;
				}
					
				cout << "(";
				for (int j = left; j < right ; j++)
				{
					cout << NE_P[j].NodeID << ",";
				}
				cout  << ")";
			}
			cout << endl;
		}
		cout << "______________________ID_HASH________________________" << endl;
		for (int i = 0; i < NodeNumber; i++)
		{
			cout << "ID:" << i << " hash:" << ID_hash[i] << "; ";
		}
		cout << endl;
	}

	void  Graph_D_H::Graph::writePartition() {

		std::string fileName;
		if (PartitionMethod == 1) {
			fileName = "./Partition/" + graphName + "/LatitudeCut.txt";
		}
		else if (PartitionMethod == 2) {
			fileName = "./Partition/" + graphName + "/MinimumCut.txt";
		}
		cout << fileName << endl;
		int tempHeight = 1;
		std::fstream heightdata;
		heightdata.open(fileName,  ios::out );
		heightdata << TreeHeight << endl;
		while (tempHeight <= TreeHeight)
		{
			int lowestIndexStart = std::pow(2, tempHeight - 1) - 1;
			int lowestIndexEnd = std::pow(2, tempHeight) - 2;
			
			for (int i = lowestIndexStart; i <= lowestIndexEnd; i++)
			{
				int left = partition_Tree[i].first;
				int right = partition_Tree[i].second;
				heightdata << tempHeight << " " << left << " " << right << "\n";
			}
			
			tempHeight++;
		}

		heightdata.close();


		if (PartitionMethod == 1) {
			fileName = "./Partition/" + graphName +  "/LatitudeID.txt";
		}
		else if (PartitionMethod == 2) {
			fileName = "./Partition/" + graphName + "/MinimumID.txt";
		}
		cout << fileName << endl;
		std::fstream heightdata_ID;
		heightdata_ID.open(fileName, ios::out);
		heightdata_ID << NodeNumber << endl;
		for (auto i = 0; i < NE_P.size(); i++) {
			heightdata_ID << i << " " << std::fixed << std::setprecision(1) <<  (double)NE_P[i].first << " " << (double)NE_P[i].second << " " << NE_P[i].NodeID << endl;
		}
		heightdata_ID.close();

	}

}