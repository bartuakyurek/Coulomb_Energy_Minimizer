# Coulomb_Energy_Minimizer
---

⚠️ This branch is solving Coulomb shapes with Fast Multipole Method (source code is taken at [here](https://github.com/kbrauss/FMM3D/tree/master)). It is supposed to conduct the optimization in O(nlgn) time; however, right now it seems to be slower than expected (possibly not O(nlgn)). 

---

<img width="581" alt="Screen Shot 2023-07-11 at 21 54 33" src="https://github.com/bartuakyurek/Coulomb_Energy_Minimizer/assets/77360680/ecc39cc7-9851-41e6-8d3c-5b652bba5b09">

Left: original decimated mesh with 769 vertices Right: Coulomb shape from iteration 300, using Fast Multipole Method (FFM).
