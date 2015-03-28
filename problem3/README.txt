This program was written with PYTHON 3 in mind. It is recommended that this program be run in IDLE, which comes with the Windows installation package for Python 3. To run from command line:
	py problem3.py

The “main” routine looks as follows:

	# Problem 3b    
	compute_leslie_populations(Leslie, x0)
	# Problem 3c
	power_method(A = Leslie, u0 = [1, 1, 1, 1, 1, 1, 1, 1, 1], tol = 1)

	# Problem 3d
	Leslie.get_A()[0][1] /= 2
	compute_leslie_populations(Leslie, x0)
	power_method(A = Leslie, u0 = [1, 1, 1, 1, 1, 1, 1, 1, 1], tol = 1)

Everything will print out to the terminal automatically. There should be no need to modify this code. However, the function compute_leslie_populations takes 2 arguments: Argument 1 is a Leslie matrix. Argument 2 is a vector x0. The function  power_method takes 3 arguments: Argument 1 is a Leslie matrix. Argument 2 is an initial guess u0. Argument 3 is a tolerance. It can be values like 0.1, 0.5, 0.99, 1, or 1.5 or above  in some cases. 
