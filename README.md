# Object-Oriented python Pseudo-Spectral Turbulence

## Pseudo-spectral collocation code for two-dimensional turbulence simulations written in object-oriented python.

Once the repository has been cloned or downloaded, check that everythink works by running the validation script
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

For a description of the theory behind this code, or to run other cases, such as a double shear layer, or decaying isotropic trubulence, look at [this]().

![shearlayer](shearlayer.png)
