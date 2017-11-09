			Digital 3D Geometry Processing 

		Exercise 6 Remeshing

All of the implemented methods work quite fine if we run them one by one. However we can see that some overlap of faces appear in the mesh. Unfortunately we still haven't found where it could come from.


	----- 1 Spliting long edges -----
We implement this first method quite intuitively. The biggest problem we've met is we first wrote "2/3" instead of "2./3." and we've spent some time to find what was the reason of the resulting wild zero. 


	----- 2 Collapsing short edges -----
Once again, we implemented the given algorithm and the result looks pretty good.


	----- 3 Equalizing valences -----
We here had to deal with many multiple conditions and checks to decide if we can flip an edge or not. This part seems tricky at the first look but is not so hard if we implement it step by step.


	----- 4 Tangential smoothing -----
This calculation requires a good ability in understanding basing matrix operations (dot product, inverse, multiplication, ...). A simple switch in the operations order produces a really bad result.


	----- 5 Adaptative remeshing -----
We split this method in three main parts. First we compute the optimal adaptative length. Then we smooth the mesh according to this length and finally we rescale the mesh.



Fraction of the total workload done:
Jonathan Collaud : 33%
Antoine Hoffmann : 33%
Valentin Rochat : 33%