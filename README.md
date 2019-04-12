Finitediff

Finitediff is a library written in Python that will automatically generate the finite difference code of a PDE into the MATLAB language. The user can easily input the PDE system they will be observing in flux form, and the finite difference code is automatically generated, as well as a driver file called driver_fd.m. In the example folder are three examples: Viscous burgers, Nonisentropic Navier-Stokes, and Reactive Navier-Stokes. In each of them, the user was only required to make the driver file in Python, and the rest of the code was automatically generated. Current requirements: Python 3+, Numpy, Scipy.

We will go through a tutorial on how we implemented burgers equation.  Step one is to import the needed packages: symbols from sympy and create_code from time_evolution, which is found in finite_difference_matlab/Core.

>> from sympy import symbols
>> from time_evolution import create_code

Next, we create a symbol for each of our variables.  Since Burger's equation is only one dimensional, we only need to create one variable. For more information about how these symbols work, see Sympy documentation.  We put all the symbols in a list as well as create an empty list for parameters.

>> u = symbols('u')
>> U = [u]
>> parameters = []

We are now ready to input the system in flux form, that is, f_0(U)_t + F_1(U)_x = (B(U)U_x)_x + G(u).  For Burger's equation, (u)_t + (u**2/2)_x - (u)_x = (u)_xx, we input 

>> f0 = [u]
>> f1 = [u**2/2-u]
>> BU = [[1]]
>> G = [0]

(Note that BU is an n x n list, where n is the size of the system.)

Now, we can call the create_code function, passing it the System variables, parameters, and flux form, with the final parameter being the file path. For more information about file paths, see Python's sys module.

>>create_code(U, parameters, f0, f1, G, BU, filePath)

At the point of running this line, the following files have automatically been generated in filePath.
driver_fd.m -- The Matlab driver file
fd_F -- The system defined as a finite difference function
fd_jac -- The jacobian of the system's finite difference code

Because the boundary conditions and initial function to be used in the finite difference study will change based on who uses the code, those are not automatically created by this program.  After generating the above files, and before running them, the user should open driver_fd and add their initial function as well as their boundary conditions. The examples in finite_difference_matlab/Examples can give insight on how to do this.

For questions about this package and its use, requests for added functionality, or request to join the development team, please email,

Jalen Morgan -- jalen.morgan000@gmail.com
David Todd -- rd.todd25@gmail.com




