This program was written with PYTHON 3 in mind. It is recommended that this program be run in IDLE, which comes with the Windows installation package for Python 3. To run from command line:
	py problem2.py

--The following code will randomly generate an x binary stream of length 5, generate matrices A0 and A1, compute vectors y0 and y1, and finally create the y binary stream. Everything will be printed out to the terminal automatically. 

	
	x_stream = c.gen_random_x_stream()
	y_stream = c.gen_y_stream(x_stream)
	c.decode_y_stream(y_stream)

--To generate an x binary stream of length n, simply change the argument for in in the constructor:

	c = ConvolutionalMatrix(n = 5)

--To run the Jacobi and Gauss-Seidel matrices with any matrix, use the following example as a guide:
m = Matrix()
m.set_matrix([[5, -2, 3], [-3, 2, 1], [2, -1, -7]])
m.set_b([-1, 2, 3])
jacobi(matrix = m, guess = [0, 0, 0, 0, 0], tol = 0.0001, decode_binary_stream = 0)
gauss_seidel(matrix = m, guess = [0, 0, 0, 0, 0], tol = 0.0001, decode_binary_stream = 0)

1) Declare a matrix. The constructor can take 0 arguments, 1 argument, or 2 arguments. Argument 1 is matrix A. Argument 2 is vector b.
1a) Alternatively, set A by calling set_matrix on the matrix object.
1b) Alternatively, set b by calling set_b on the matrix object.
2) Call the jacobi and gauss_seidel functions. They each take the following arguments: Argument 1 is the matrix object. Argument 2 is the initial guess. Argument 3 is the tolerance. Argument 4 MUST BE 0 for this to work for a regular matrix. The variable decode_binary_stream is a flag that tells the program to modulo every element in the matrix by 2. It is used for convolutional codes.
