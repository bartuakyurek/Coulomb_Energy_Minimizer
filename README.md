# Coulomb Energy Minimizer

This repo is an unofficial implementation of the paper *Coulomb shapes: using electrostatic forces for deformation-invariant shape representation*. It currently has two branches: 

    1-) `main`: This branch is a direct (naive) implementation of eqn.1, 4, and 6. It has $O(n^2)$ time complexity.
    
    2-) `fmm`: This branch is the actual implementation of the original paper. 
        It uses Fast Multipole Method as also used by the paper, with $O(nlgn)$ time complexity.

<img width="623" alt="16_decimated-iter-600-comparison" src="https://github.com/bartuakyurek/Coulomb_Energy_Minimizer/assets/77360680/b41a3a0b-95bf-4ab5-96f9-721f98f4bc61">

Left: original decimated mesh (16.obj) Right: Coulomb shape of 600 iterations 


---

⚠️ This naive implementation of Coulomb shapes has $O(n^2)$ time complexity. See `fmm` branch for an $O(nlgn)$ time implementation.

---
