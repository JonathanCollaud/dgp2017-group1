			Digital 3D Geometry Processing 

		Exercise 8 - Printing

	----- Creating the dual mesh -----
In this first part the main problem we encoutered was that we were adding five times each vertex because we werent looking at if it was already added. Which ended up to be not to hard to resolve. For the rest your indications really helped us to do it correctly.


	----- Replacing the edges by cylinders -----
Then the add_primitive function was quite straightforward using your hints. However the create_cylinder_mesh gave us lot's of trouve to find the right way of doing the rotation of the cylindre.
Firstly because it would have been good to be a little clearer on the exact definition of the cylindre (exact position in space and length).
Secondly because our first try was to make it rotate using two perpendicular axis which ended up not working.
Finally we decided to compute the vector perpendicular to the edge and the cylindre and calculating the angle between the same two. Finally we rotated the vector in the plane formed by the two vectors using the perpendicular one.

	----- Verifying 3D printability -----
As we can see in screenshot1.png, the mesh seems to be 3D printable.



Fraction of the total workload done:
Jonathan Collaud : 33%
Antoine Hoffmann : 33%
Valentin Rochat : 33%