# Object-Oriented python Pseudo-Spectral Turbulence

## Pseudo-spectral collocation code for two-dimensional turbulence simulations written in object-oriented python.

Once the repository has been cloned or downloaded, check that everything works by running the validation script
```
python valid.py
```
This runs a simulation of the [Taylor-Green Vortex](https://en.wikipedia.org/wiki/Taylor%E2%80%93Green_vortex) for which we have an analytical solution. The output should look similar to this
```
Starting integrating on field.

<<<<<<< HEAD
Iteration 	 100, time 	 0.020, time remaining 	 0.079. TKE: 0.1802  ENS: 23580.54
Iteration 	 200, time 	 0.041, time remaining 	 0.058. TKE: 0.1294  ENS: 16934.10
Iteration 	 300, time 	 0.061, time remaining 	 0.038. TKE: 0.0929  ENS: 12157.35
Iteration 	 400, time 	 0.082, time remaining 	 0.017. TKE: 0.0666  ENS: 8725.782

Execution time for 484 iterations is 1.88 seconds.
The L2-norm of the Error in the Taylor-Green vortex on a 128x128 grid is 1.36e-10.
The Linf-norm of the Error in the Taylor-Green vortex on a 128x128 grid is 2.72e-10.
=======
Iteration 	 100, time 	 0.020, time remaining 	 0.079. TKE: 0.180  ENS: 23580.54
Iteration 	 200, time 	 0.041, time remaining 	 0.058. TKE: 0.129  ENS: 16934.10
Iteration 	 300, time 	 0.061, time remaining 	 0.038. TKE: 0.092  ENS: 12157.35
Iteration 	 400, time 	 0.082, time remaining 	 0.017. TKE: 0.066  ENS: 8725.782

Execution time for 484 iterations is 2.737 seconds.
The L2-norm of the Error in the Taylor-Green vortex on a 128x128 grid is 1.362e-10.
The Linf-norm of the Error in the Taylor-Green vortex on a 128x128 grid is 2.725e-10.
>>>>>>> c6ef2771b957183b7abfd84d82b96eb71875def1
```
Check that you get an error in both norms that are close to machine accuracy.

---

For a description of the theory behind this code, or to run other cases, such as a double shear layer, or decaying isotropic turbulence, look at [this]().

<p align="center">
 <img src="shearlayer.png" width="400"> 
</p>

## Using the code

This repository contains three (as of today) branches, the `master` branch that includes the pseudo-spectral code, the `CDS` branch that contains the high-order (up to 6th) compact difference scheme and another pseudo-spectral method using `pyfftw`. The last branch is not working yet. When you download the repository, you will be by default on the `master` branch. To switch to the central-difference branch use
```
git branch -a
```
to show all the remote branches, and then 
```
git checkout CDS or git checkout fftw
```
to switch to the desired branch. The solver has been implemented such that functions calls perform exact same tasks on each branches.

The simulation is initialized by defining a grid, and specifying the Reynolds number of the flow
```python
flow = Fluid(nx=64, ny=64, Re=1)
flow.init_solver()
flow.init_field("Taylor-Green")
```
Here we have initialized the Taylor-Green vortex. The solver initiation generates all the working arrays and transforms the required initial conditions. Simulations can also be initialized using results from previous simulations (these need to have been saved with `flow.save_vort("PATH/", ID)`)
```python
q = np.genfromtxt("PATH/vort_ID.dat")
flow.time = q[0, 0]
flow.init_field(field=q[1:, :])
```
here we reset the flow timer using the time value saved in the `vort_ID.dat` file. The `finish` time of the simulation must be adjusted accordingly, as well as the `ID` if the field is saved. This allows user-generated field to be used, within the limitations of the method (periodic). The main loop of the solver is called as
```python
# loop to solve
while(flow.time<=finish):
    flow.update()
    if(flow.it % 500 == 0):
        print("Iteration \t %d, time \t %f, time remaining \t %f. TKE: %f, ENS: %f" %(flow.it,
<<<<<<< HEAD
                flow.time, finish-flow.time, flow.tke(), flow.enstrophy()))
=======
              flow.time, finish-flow.time, flow.tke(), flow.enstrophy()))
>>>>>>> c6ef2771b957183b7abfd84d82b96eb71875def1
```
Small simulations can also be run live, that is showing the evolution of the vorticity field
```python
flow.run_live(finish, every=100)
```