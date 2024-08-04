# PartitionForGPUH2H
The Graph Partitioning For GPU-H2H

## Compile

Please run 
```nvcc -O3 -diag-suppress 20012 -arch=sm_60 -o Partition ./*.cu```
to compile our code.


## Data preparation

Graph files for testing like NY, FLA, E can be found at http://www.diag.uniroma1.it/~challenge9/download.shtml.
Directory ./Graph/USA covers a test network 13. 

The format of Edge file Edge.txt is:
```
13 44               //Node Number, Edge Number;
a 1 2 4             //Edge (1, 2) with weight = 4;
......
```

The format of Node file Node.txt is:
```
13                      //Node Number;
v 1 2.0 10.1
v 2 4.0 10.2            //vertex 2 and it's (Longitude, Latitude) = (4.0, 10.2)
......
```


## Run

```
./Partition GraphName PartitionHeight PartitionMethod
```

'GraphName' represents the graph;

'Partition Height' represents total partition Tree Height you used;

'PartitionMethod' determines the partition method, if == 1, then use latitude method, and == 2 uses minimum cut with beta = 0.4;

For example, ```./Partition 13 4 1```
will run latitude-longitude partition method on graph 13 and the partition tree height = 4.
Those results will be recorded at './Partition/13/xxxCut.txt' and './Partition/13/xxxID.txt'. 

Further more, we will calculate NUB size.
The inner vertex rate in each partition layer will be recorded at  'xxxPartitionRate.csv'.
