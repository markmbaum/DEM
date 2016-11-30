This repository contains MATLAB code oringally used for an [undergraduate thesis project](https://www.overleaf.com/read/kzsvpynnmtgx) (give the document time to compile, it'll come through) examining the deformation of tectonic scale shear zones. Significant changes have been made to the code currently stored here, and the most recent commit may not be a stable version.

The code is a numerical model that assembles arrays of elastic cylinders and subjects them to simple shear. The image below shows a trial of the model where the elastic elements are uniformly sized and the friction threshold governing slip between elements is very high, causing almost all of the displacement and rotation to occur on the edges of the array.

![Movie_friction_1_ncols_16_sh1](Other Things/Movie_friction_1_ncols_16_sh1.gif)

There are also two sub-projects embedded in the modeling code.

The first sub-project is an algorithm for assembling rectangularly confined collections of randomly sized circles, which can be used for the simulations. The size of the circles is random but constrained between user defined limits. The algorithm uses simple geometric relationships and root finding methods (mostly multi-dimensional Newton's method) to create a rectangular border of circles then fill it in with more circles, all of them perfectly adjacent.

![demo-gif](Other Things/demo-gif.gif)

The second sub-project is set of polynomials for mapping collections of indices, from an ordered list or array, onto single indices for efficient storage of data relating the collections. The polynomials essentially hash all possible combinations of integer indices in an array to single, unique values. A quick write-up of this concept can be viewed here: [https://www.overleaf.com/read/dpknzcyjxxtt](https://www.overleaf.com/read/dpknzcyjxxtt) (wait and refresh the page if the document doesn't render immediately).
