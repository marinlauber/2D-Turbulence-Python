# Object-Oriented python Pseudo-Spectral Turbulence

## Pseudo-spectral collocation code for two-dimensional turbulence simulations written in object-oriented python.

Once the repository has been cloned or downloaded, check that everything works by running the validation script
```
python valid.py
```
This runs a simulation of the [Taylor-Green Vortex](https://en.wikipedia.org/wiki/Taylor%E2%80%93Green_vortex) for which we have an analytical solution. The output should look similar to this
```
Starting interating on field.

Iteration  200, time 0.010275, time remaining 0.089725
Iteration  400, time 0.021411, time remaining 0.078589
Iteration  600, time 0.033636, time remaining 0.066364
Iteration  800, time 0.047186, time remaining 0.052814
Iteration 1000, time 0.062385, time remaining 0.037615
Iteration 1200, time 0.079689, time remaining 0.020311
Iteration 1400, time 0.099778, time remaining 0.000222

Execution time for 1403 iterations is 31.427155 seconds.
The L2-norm of the Error in the Taylor-Green vortex on a 128x128 grid is 7.592437e-12.
The Linf-norm of the Error in the Taylor-Green vortex on a 128x128 grid is 1.518541e-11.
```
Check that you get error in both norms that are close to machine accuracy.

---

For a description of the theory behind this code, or to run other cases, such as a double shear layer, or decaying isotropic turbulence, look at [this]().

<p align="center">
    <img src="shearlayer.png" width="400"> 
</p>

## Using the code

The simulation is initialized by defining a grid, and specifying the Reynolds number of the flow
```python
flow = Fluid(nx=64, ny=64, Re=1)
flow.init_field("Taylor-Green")
flow.init_solver()
```
Here we have initialized the Taylor-Green vortex. The solver initiation generates all the working arrays and transforms the required initial conditions. Simulations can also be initialized using results from previous simulations (these need to have been saved with  `flow.save_vort("PATH/", ID)`)
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
    iterr += 1
    if(iterr % 1000 == 0):
        print("Iteration \t %d, time \t %f, time remaining \t %f" %(iterr, flow.time, finish-flow.time))
```
Small simulations can also be run live, that is showing the evolution of the vorticity field
```python
flow.run_live(finish, every=100)
```
