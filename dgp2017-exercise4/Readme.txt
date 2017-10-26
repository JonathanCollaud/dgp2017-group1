			Digital 3D Geometry Processing 

		Exercise 4 Surface Normals, Curvature

	----- 1 Computing Vertex Normals -----

In this first part the main difficulty was to understand and applicate the basics of the mesh naming, 
functions and properties. Once thoses problems overpassed the two first points were quite straight forwardly
computed. While implementing the third one we firstly tried to use the crossproduct to compute the angle.
However this led to lot's of problem which were probably linked to the asin() function. The values were
sometimes out of computable bundary. This problem has been solved simply by using the dot product which
uses the acos() function.


	----- 2.1 Uniform Laplace Operator -----

The first problem encountered, which is a common problem with the points 2.2 and 2.3, is the boundary 
detection. We haven't a efficient way to implement it. However concerning the uniform laplace operator, it was not too hard 
to implement because of the low complexity of the formula.


	----- 2.2 Laplace-Beltrami Curvature -----

In this second part we firstly had to understand what was done in the already computed part. Then the 
major difficulty came from the fact that in this case we do a circulator over the half edges, which
was a good exercice to better understand the general way of getting the wanted informations on the 
mesh.


	----- 2.3 Gaussian Curvature -----

Using your hints and the angle formula computed in the first exercice, we only encountered small issues
while implementig it. For exemple we had a double incrementation problem which finally took us long to
find out.


	----- General comment and sphere artefacts -----

Concerning the artefacts appearing on the sphere while showing the curvature, it is at first astonishing that some
artefacts appears because a sphere has a constant curvature. However looking at the mesh at thoses places we see that
it is formed of pentagon where it is formed of hexagon everywhere else.




Fraction of the total workload done:
Jonathan Collaud : 33%
Antoine Hoffmann : 33%
Valentin Rochat : 33%