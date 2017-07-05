# Fast Calculation of Max-Min Fair Rates(MMF) for Multi-commodity Flows in Fat-tree Networks

This repository containts the implementation of various methods devised and published in IEEE Cluster'2015 conference. For details on the backgrounds and imeplementation, please refer to our publication [Cluster15]: http://ieeexplore.ieee.org/document/7307603/


##Abstract

Max-min fairness is often used in the performance modeling of interconnection networks. Existing methods to compute max-min fair rates for multi-commodity flows have high complexity and are computationally infeasible for large networks. In this work, we show that by considering topological features, this problem can be solved efficiently for the fat-tree topology that is widely used in data centers and high performance computing clusters. Using two new algorithms that we developed, we demonstrate it is possible to find the max-min fair rate allocation for multi-commodity flows in fat-tree networks that support tens of thousands of nodes. We evaluate the run-time performance of the proposed algorithms and demonstrate an application.

![alt text][logo]

[logo]: https://github.com/amollah/fast-mmf-fattree/blob/master/443122_xgft.svg "Example Fat-tree Topology"

## Usage

after download/cloning, type 'make' to build the rate calculation rate. The executable to run the tool is "fast_mmf.out".

the command line syntax is as follows:

 _./fast_mmf.out [xgft/pgft] h m1 m2 ... mh w1 w2 ...wh p1 p2 ... ph [GEN/LP/PF/DMODK/OPT] trafficfile_
 
 **h** is the height of the fat-tree. The tool currently supports 2-level and 3-level generalized fat-trees.
 **m1 m2 ... mh** lists the number of  nodes of  children/successor each switch level (excluding leaf level)
 **w1 w2 ...wh** lists the number of parent/predecessor  nodes of each fat-tree level (excluding root level)
**p1 p2 ... ph** lists the number of parallel links used to connect switches of adjacent levels
**GEN** uses the Generic MMF rate calculation algorithm 
**LP** uses the Linear Programming based MMF rate calculation algorithm for fat trees MMF\_PGFT\_LP.
**PF**  uses the Progressive Filling based MMF rate calculation algorithm for fat trees MMF\_PGFT\_PF.
**DMODK**  solves the Progressive Filling based MMF Simple Routing Problem using destination-mod-k routing
**OPT**  uses the optimized MMF rate calculation algorithm based on sub-fat-tree saturation MMF\_PGFT\_OPT.

## Citation
If you use any materials made available here in your work, please refer to our work using the following citation(s)

### Plain Text 
M. A. Mollah, X. Yuan, S. Pakin and M. Lang, "Fast Calculation of Max-Min Fair Rates for Multi-commodity Flows in Fat-Tree Networks," 2015 IEEE International Conference on Cluster Computing, Chicago, IL, 2015, pp. 351-360.
doi: 10.1109/CLUSTER.2015.56

### BibTeX
@INPROCEEDINGS{7307603, \
author={M. A. Mollah and X. Yuan and S. Pakin and M. Lang}, \
booktitle={2015 IEEE International Conference on Cluster Computing}, \
title={Fast Calculation of Max-Min Fair Rates for Multi-commodity Flows in Fat-Tree Networks}, \
year={2015}, \
pages={351-360}, \
keywords={computational complexity;multiprocessor interconnection networks;telecommunication network topology;data centers;fat-tree networks;fat-tree topology;high complexity;high performance computing clusters;interconnection networks;max-min fair rate allocation;max-min fairness;multicommodity flows;performance modeling;run-time performance;topological features;Downlink;Network topology;Resource management;Routing;Topology;Uplink;Vegetation;fat-tree topology;high-performance computing;max-min fairness}, \
doi={10.1109/CLUSTER.2015.56}, \
ISSN={1552-5244}, \
month={Sept},} \
