This folder contains the implementation of three similar notebooks (and their respective implemented .py pipeline) for persistent homology done based on 'road importance' calculated using edge betweeness centrality (see calculate_road_importance folder).

The idea is the same for the three notebooks, but the filtration choice changes.

1. Log approach: Take the log of the importance and refine the first stps as to only get the most important roads. use binning for the rest.
2. Deciles approach: Take 10% (or any equal percentage , for that matter) of streets inside each step and perform homology
3. Smoothen the edge betweeness by using a Gaussian filter. The filtration choice is the same as (1.) 

This pipelines follow on the work of Corcoran & Jones, by also building Persistent Densities 