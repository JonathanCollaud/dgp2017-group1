			Digital 3D Geometry Processing 

		Exercice 2 - Geometry Representations

	----- 2 Coding Exercice -----

In order to implement an equivalent of the marching squares algorithm but on triangles we started by looking at the value of the function at each point. Then we compared the value of each vertex to the other vertex and look if the sign is changing. A particular case has been implement for the case when a vertex is a zero. Then we just had to use a linear approximation to find the position of the zero on each line and like this obtain the point where the funciton is zero. Finally the two points were pushed back in the right array to let the already implemented algorithm draw the wanted lines.



Fraction of the total workload done:
Jonathan Collaud : 33%
Antoine Hoffmann : 33%
Valentin Rochat : 33%