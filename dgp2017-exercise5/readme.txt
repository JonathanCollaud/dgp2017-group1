			Digital 3D Geometry Processing 

		Exercise 5 Smoothing and Feature Enhancement

	----- 1 Explicit smoothing -----

Both functions smooth() and uniform_smooth() were quite intuitive and simple to implement. Moreover, they seem to work properly.


	----- 2 Implicit smoothing -----

We faced a problem that smooth the bunny into a slug. We haven't found what goes wrong in our code until yet.


	----- 3 Feature Enhancement -----

Both of the enhancement are quite aggressive but actually performs an enhancement. We've got the problem that the initial positions are modified during the execution. Hence the difference was inexistent and the enhancement didn't do anything. We then had to found a solution to copy deeply the values of the old positions.


Fraction of the total workload done:
Jonathan Collaud : 33%
Antoine Hoffmann : 33%
Valentin Rochat : 33%