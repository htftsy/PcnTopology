Strengthening decentralized payment network by duty channels
==

Structure
-
`/contract` includes the demonstrative smart contract in *Solidity*.

`/eval` includes the following program files for simulations.
* `lps.h` builds the LPS graph according to $(p,q)$ as specified in the paper.
* `withGraphs.h` builds the random graphs as specified in the paper.
* `simu.cpp` evaluates our concerned metrics when inputted $(p,q,\gamma)$ as specified.
* `fees.cpp` simulates user revenues when inputted $(p,q,\gamma)$ as specified.
* `dataGen.cpp` generates batched datas by triggering `simu.cpp` for the following $(p,q)$ pairs, each for $6$ times:
<p style="text-align:center;">(5, 23) (7, 23) (11, 23) (17, 23) (19, 23) (5, 47) (11, 47) (13, 47) (19, 47) (23, 47).</p>

How to build
-
```bash
cd eval
g++ -o simu simu.cpp
./simu
```

