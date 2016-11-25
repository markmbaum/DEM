This repository contains the MATLAB code for an undergraduate [thesis project](https://www.overleaf.com/read/kzsvpynnmtgx) examining the deformation of tectonic scale shear zones. The code is a numerical model that assembles arrays of elastic cylinders and subjects them to simple shear, as shown below (the gif may take a moment to load).

![Movie_friction_1_ncols_16_sh1](Other Things/Movie_friction_1_ncols_16_sh1.gif)

There are also two sub-projects embedded in the modeling code:

The first sub-project is an algorithm for assembling rectangularly confined collections of randomly sized circles. The size of the circles is random but constrained between user defined limits. The algorithm uses simple geometric relationships and root finding methods (mostly multi-dimensional Newton's method) to create a rectangular border of circles then fill it in with more circles, all of them perfectly adjacent.

![demo-gif](Other Things/demo-gif.gif)

The second sub-project is set of polynomials for mapping collections of indices, from an ordered list or array, onto single indices for efficient storage of data relating the collections. A quick write-up of this concept can be found here: [https://www.overleaf.com/read/dpknzcyjxxtt](https://www.overleaf.com/read/dpknzcyjxxtt) (wait and refresh the page if the document doesn't render immediately).
