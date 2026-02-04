Strengthening Multi-hop Channels via Strategic Mesh Connections
==

Artifact DOI
-
https://doi.org/10.1007/978-3-032-07024-1_2

Structure
-
`/contract` includes the demonstrative smart contract in *Solidity*.

`/eval` includes the following program files for simulations.
* `lps.h` builds the LPS graph according to $(p,q)$ as specified in the paper.
* `withGraphs.h` builds the random graphs as specified.
* `simu.cpp` evaluates our concerned metrics given $(p,q,\gamma)$.
* `fees.cpp` simulates user revenues with $(p,q,\gamma)$.
* `dataGen.cpp` generates batched data by triggering `simu.cpp` for the following $(p,q)$ pairs, each for $6$ times:
<p align="center">(5, 23) (7, 23) (11, 23) (17, 23) (19, 23) (5, 47) (11, 47) (13, 47) (19, 47) (23, 47).</p>

How to build
-
```bash
cd eval
g++ -o simu simu.cpp
./simu
```

