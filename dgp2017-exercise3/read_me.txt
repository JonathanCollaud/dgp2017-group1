			Digital 3D Geometry Processing 

		Exercice 3 - Curves

	----- 2 Coding Exercice -----

laplacianSmoothing()

Concerning that first part, the implementation of the smoothing has been done quite easily 
following the given indications. However we also had to uniformly scale the curve to it's 
original length. This has been a more complex point because we counldn't find the perfect 
solution to keep the total length, the distance between each vertex and not do disturb the 
main algorithm while doing it. We finally endend up with a solution which is not perfect 
but we weren't able to find a better one. So the solution we choose is to calculate a center 
of the points and to make them iteratively go a bit further to that center until the length 
has been fully restored. This does keep the total length (up to a certain tolerance) but does
influence the distance between the vertexes and does also influence the relative position of 
each vertex in comparison the the others. This does influence a bit the main algorithm.


osculatingCircle()

For this part, the first complication was to find the center of the circle created by the three
points. We did it as it would be done geometrically, by finding the mid point between two vertexes
and the vector perpendicular to the vector linking them. We did it for two of the three links and
looked at the intersection of the two perpendicular vector passing through the mid points. We also 
included the special cases when one of the perpendicular vector is vertical. Concerning the length 
correction, we used the same algorithm used in the first part. It is also interesting to mention 
that two differents epsilon have been define for each parts because the one for the osculating circle
has to be much smaller.




Fraction of the total workload done:
Jonathan Collaud : 33%
Antoine Hoffmann : 33%
Valentin Rochat : 33%