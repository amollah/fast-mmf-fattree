# Fast Calculation of Max-Min Fair Rates(MMF) for Multi-commodity Flows in Fat-tree Networks

This repository containts the implementation of various methods devised and published in IEEE Cluster'2015 conference. For details on the backgrounds and imeplementation, please refer to our publication [Cluster15]: http://ieeexplore.ieee.org/document/7307603/


##Abstract

Max-min fairness is often used in the performance modeling of interconnection networks. Existing methods to compute max-min fair rates for multi-commodity flows have high complexity and are computationally infeasible for large networks. In this work, we show that by considering topological features, this problem can be solved efficiently for the fat-tree topology that is widely used in data centers and high performance computing clusters. Using two new algorithms that we developed, we demonstrate it is possible to find the max-min fair rate allocation for multi-commodity flows in fat-tree networks that support tens of thousands of nodes. We evaluate the run-time performance of the proposed algorithms and demonstrate an application.


##Usage

after download/cloning, type 'make' to build the rate calculation rate. The executable to run the tool is "fast_mmf.out"

##Citation
If you use any materials made available here in your work, please refer to our work using the following citation(s)

###Plain Text
M. A. Mollah, X. Yuan, S. Pakin and M. Lang, "Fast Calculation of Max-Min Fair Rates for Multi-commodity Flows in Fat-Tree Networks," 2015 IEEE International Conference on Cluster Computing, Chicago, IL, 2015, pp. 351-360.
doi: 10.1109/CLUSTER.2015.56

###BibTeX
@INPROCEEDINGS{7307603, 
author={M. A. Mollah and X. Yuan and S. Pakin and M. Lang}, 
booktitle={2015 IEEE International Conference on Cluster Computing}, 
title={Fast Calculation of Max-Min Fair Rates for Multi-commodity Flows in Fat-Tree Networks}, 
year={2015}, 
pages={351-360}, 
keywords={computational complexity;multiprocessor interconnection networks;telecommunication network topology;data centers;fat-tree networks;fat-tree topology;high complexity;high performance computing clusters;interconnection networks;max-min fair rate allocation;max-min fairness;multicommodity flows;performance modeling;run-time performance;topological features;Downlink;Network topology;Resource management;Routing;Topology;Uplink;Vegetation;fat-tree topology;high-performance computing;max-min fairness}, 
doi={10.1109/CLUSTER.2015.56}, 
ISSN={1552-5244}, 
month={Sept},}
