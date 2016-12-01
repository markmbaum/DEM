This repository contains MATLAB code for an undergraduate thesis project (available [here](https://www.overleaf.com/read/kzsvpynnmtgx)) examining the deformation of tectonic scale shear zones. The code is a numerical model that assembles arrays of elastic cylinders and subjects them to simple shear. The animated images below show trials of the model, viewed from above, with uniformly sized cylinders. Compared to the first, the second image shows a trial with slightly larger elements and a higher friction threshold controlling slip between elements, resulting in more deformation away from the midline of the array.

![gif_friction_p2_ncols_22_sh1](Other Things/gif_friction_p2_ncols_22_sh1.gif)

![gif_friction_p4_ncols_16_sh1](Other Things/gif_friction_p4_ncols_16_sh1.gif)

There are also two sub-projects.

The first sub-project is an algorithm for assembling rectangularly confined collections of randomly sized circles. The size of the circles is random but constrained between user defined limits. The algorithm uses simple geometric relationships and root finding methods (mostly multi-dimensional Newton's method) to create a rectangular border of circles then fill it in with more circles, all of them perfectly adjacent. The resulting array of circles can then be subjected to simple shear like in the animated image above, to more realistically simulate the complexity of fault systems in shear zones. A demonstration of this circle assembling process is shown below.

![demo-gif](Other Things/demo-gif.gif)

The second sub-project is set of polynomials for mapping collections of indices, from an ordered list or array, onto single indices for efficient storage of data relating those collections. A quick write-up and more extensive explanation of this concept can be found here: [https://www.overleaf.com/read/dpknzcyjxxtt](https://www.overleaf.com/read/dpknzcyjxxtt) (wait and refresh the page if the document doesn't render immediately).
