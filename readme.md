Sage implementations for the optimization of the attack parameters and security estimates for the classical and quantum hybrid attack provided in [2] by Thomas Wunderer. Parts of the implementation are based on a previous implementation for the security estimates of [1] by Florian Goepfert, Rachel Player, and Thomas Wunderer.


Licence:
Public domain.

Note: The original code was provided by Thomas Wunderer. Rachel Player has edited the code to make the following changes:
1. provide documentation via comments
2. add additional cost models for BKZ
3. rename variables to align with the notation in [2]
4. estimate security against the hybrid attack of typical homomorphic encryption parameters, documented in Wunderer-analysis-FHE.pdf

Independently and concurrently to the work described in point 4, Son and Cheon [3] also applied a Wunderer-style analysis of the hybrid attack in the FHE parameter setting.

[1] Johannes A. Buchmann, Florian Göpfert, Rachel Player, and Thomas Wunderer. On the hardness of LWE with binary error: Revisiting the hybrid lattice-reduction and meet-in-the-middle attack. In Progress in Cryptology - AFRICACRYPT 2016 - 8th International Conference on Cryptology in Africa, Fes, Morocco, April 13-15, 2016, Proceedings, pages 24 - 43, 2016.

[2] Thomas Wunderer. On the Security of Lattice-Based Cryptography Against Lattice Reduction and Hybrid Attacks. PhD thesis, Darmstadt University of Technology, Germany, 2018.

[3] Yongha Son and Jung Hee Cheon. Revisiting the Hybrid attack on sparse and ternary secret LWE. WAHC 2019, to appear.