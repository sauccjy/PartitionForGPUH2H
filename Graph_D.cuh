#pragma once
//#pragma auto_inline(off)
#include<iostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include<thrust/device_vector.h>
#include<thrust/host_vector.h>
#include<thrust/count.h>
#include<fstream>
#include<thrust/sort.h>
#include<thrust/shuffle.h>
#include<thrust/random.h>
#include<thrust/execution_policy.h>
#include<stack>
#include<queue>
#include<algorithm>
#include<chrono>
#include<random>
#include<vector>
#include<map>
#include<set>
#include<unordered_map>
#include<string>
#include <thread>
#include <iomanip>


using namespace std;
using namespace std::chrono;
namespace Graph_D_H
{

	struct time_Mine {
		std::chrono::high_resolution_clock::time_point start;
		std::chrono::high_resolution_clock::time_point end;
		void updateStart()
		{
			start = std::chrono::high_resolution_clock::now();
		}
		void updateEnd()
		{
			end = std::chrono::high_resolution_clock::now();
		}
		long long get_microsecond_duration()
		{
			return std::chrono::duration_cast<microseconds>(end - start).count();
		}

	};

	struct pairs
	{
		int first = -1;
		int second = -1;

		__host__ __device__
		pairs() = default;
		__host__ __device__
		pairs(int first, int second)
			: first(first), second(second)
		{
		}
		__host__ __device__
			pairs(const pairs& other)
			: first(other.first), second(other.second)
		{
		}
		__host__ __device__
		void pairsCopy(const pairs& others)
		{
			first = others.first;
			second = others.second;
		}
		__host__ __device__
		void pairsReset(int first_, int second_)
		{
			first = first_;
			second = second_;
		}
		__host__ __device__
		bool operator<(const pairs& other) const {
			if (first < other.first)
				return true;
			else if (first == other.first)
				return second < other.second;
			else
				return false;
		}
		__host__ __device__
			bool operator>(const pairs& other) const {
			if (first > other.first)
				return true;
			else if (first == other.first)
				return second > other.second;
			else
				return false;
		}
	};

	struct TDrank {
		int first = -1;//Node degree
		int second = -1;//Node ID
		int third = 0;  //continus degree

		__host__ __device__
		TDrank() = default;
		__host__ __device__
		TDrank(int first, int second, int third) : first(first), second(second), third(third){}
		__host__ __device__
		TDrank(const TDrank& others): first(others.first), second(others.second), third(others.third) {}
		__host__ __device__
			void setDetal(int f,int s,int t)
		{
			first = f;
			second = s;
			third = t;
		}
		__host__ __device__
			bool operator<(const TDrank& other) const
		{
			if (third + first != other.third + other.first) {
				return third + first < other.third + other.first;
			}
			else {
				if (first != other.first) {
					return first < other.first;
				}
				else {
					return second < other.second;
				}
			}
		}
		__host__ __device__
			bool operator>(const TDrank& other) const
		{
			if (third + first != other.third + other.first) {
				return third + first > other.third + other.first;
			}
			else {
				if (first != other.first) {
					return first > other.first;
				}
				else {
					return second > other.second;
				}
			}
		}
		__host__ __device__
			void operator=(const TDrank& other)
		{
			this->first = other.first;
			second = other.second;
			third = other.third;

		}
		__host__ __device__
			bool operator==(const TDrank& other)const
		{
			return second == other.second;
		}
	};

	static vector<pairs> CHInfo;
	struct HostRank {
		int x;
		bool operator<(const HostRank& other) const
		{
			if (CHInfo[x].first + CHInfo[x].second != CHInfo[other.x].first + CHInfo[other.x].second)
			{
				return CHInfo[x].first + CHInfo[x].second < CHInfo[other.x].first + CHInfo[other.x].second;
			}
			else
			{
				if (CHInfo[x].first != CHInfo[other.x].first)
				{
					return CHInfo[x].first < CHInfo[other.x].first;
				}
				else
				{
					return x < other.x;
				}
			}
		}
		HostRank(int x) : x(x) {
		}
	};

	template<typename T>
	struct myPair {
		int NodeID = -1;
		T first = 0;
		T second = 0;
		//__host__ __device__
		myPair() = default;
		__host__ __device__
		myPair(const T first, const T second) : first(first), second(second) {}
		__host__ __device__
		myPair(const T first, const T second, int ID) : first(first), second(second) ,NodeID(ID){}
		__host__ __device__
			myPair(const myPair<T>& other ) : first(other.first), second(other.second), NodeID(other.NodeID) {}
		__host__ __device__
		void setPairs( int ID,const T& _first, const T& _second) {
			NodeID = ID;
			first = _first;
			second = _second;
		}

		__host__ __device__
			bool operator<(const myPair<T>& others)const {
			return first < others.first;
		}
		__host__ __device__
			void operator=(const myPair<T>& other)
		{
			this->first = other.first;
			second = other.second;
			NodeID = other.NodeID;

		}
	};

	struct myPair_latitude_less_than
	{
			bool operator()(const myPair<double> &other, const myPair<double> &other1) const
		{
			return other.second < other1.second;
		}
	};

	struct myPair_longitude_less_than
	{
			bool operator()(const myPair<double> &other, const myPair<double> &other1) const
		{
			return other.first < other1.first;
		}
	};

	class Graph {

	public:

		//main info
		int NodeNumber = 0;
		int EdgeNumber = 0;
		string NodeFile = "";
		string EdgeFile = "";
		string graphName = "";

		//adjList and  NE(longitude and latitude)
		thrust::host_vector<thrust::host_vector<pairs> > adjList = { };
		thrust::host_vector<myPair<double> > NE_P = {};//first : N ; second :E ; nodeID : ID;

		thrust::host_vector<int> ID_hash = {};//hash NE_P
		thrust::host_vector<int> ranks = {};//ranks[ID] means ID`s rank,from ID to rank,rank 0 means lowest rank
		int head = -1; //headID,
		//partitionTree
		thrust::host_vector<pairs> partition_Tree = {}; // partition (left,right) index in NE_P , binary tree structure
		int TreeHeight = 1;
		int lowestMaxSize = 0;
		int CHTreeHeight = 0;
		//order and partition
		//Partition Rank Tree
		thrust::host_vector<TDrank> nonAdjcentNode = {};// ID = nonAdjcentNode[rank].second, from rank to ID
		thrust::host_vector<myPair<int>> NANHash = {}; //first : nonAdjcentNode index; second maxSize; nodeID : nonadjcent Size
		thrust::host_vector<int> Reflash_for_PRT = {};
		//adjList and  NE(longitude and latitude)
		//thrust::host_vector<thrust::host_vector<pairs> > adjList = { };
		//thrust::host_vector<myPair<double> > NE_P = {};//first : N ; second :E ; nodeID : ID;

		//thrust::host_vector<int> ID_hash = {};//hash NE_P
		//partitionTree
		//thrust::host_vector<pairs> partition_Tree = {}; // partition (left,right) index in NE_P , binary tree structure
		//int TreeHeight = 1;
		//int lowestMaxSize = 0;

		thrust::host_vector<myPair<int>> CHTreeHash = {};
		thrust::host_vector<bool> visited = {};
		thrust::host_vector<myPair<int>> CHTree = {}; //first: outPoint; second: weight; NodeID: hub;
		long long int chSize = 0;//in byte
		//construct data and time info

		int changeHeight = 0;
		long long int partitionTime = 0;
		int PartitionMethod = 1;
		double Beta = 0.4;

		//fundamental
		Graph(const string graphName, const string Nodefile, const string EdgeFile, int threadPoolSize, int changeHeight);
		~Graph();
		void checklink();
		
		//latitude cut functions
		int partition(int goalHeight);
		int partition_Latitude_first(int goalHeight);
		//minimum cut functions
		int partition_minumum_cut(int goalHeight, double beta);
		int generateBFS(std::map<int, std::map<int, int>>& graph, std::map<int, bool>& visited_s, std::map<int, bool>& visited_t, int s, double beta);
		void generateFarthestVertex(std::map<int, std::map<int, int>>& standard_graph, int s, int& f);
		bool bfs_for_edmondsKarp(const std::map<int, std::map<int, int>>& residualGraph, int source, int sink, std::map<int, int>& parent);
		int edmondsKarp(std::map<int, std::map<int, int>>& graph, int source, int sink, std::map<int, std::map<int, int>>& residualGraph);
		void findMinCut(const std::map<int, std::map<int, int>>& residualGraph, int source, std::map<int, bool>& visited);
		void findMinCutAndSplitGraph(std::map<int, std::map<int, int>>& graph, int source, int sink, std::set<int>& source_paritition, std::set<int>& sink_paritition);
		void makeAdjcentNode_2(int goalHeight);
		void displayPartition();
		void writePartition();

	};

}