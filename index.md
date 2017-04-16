## Proposal
### Summary
We are going to implement an interactive, high-resolution fluid simulation.

### Background
The simulation takes discrete time steps, recalculating the velocities and pressures (which are governed by the Navier-Stokes equations) for every pixel each step, which lends itself perfectly to a parallel approach. By rendering each time step as fast as possible, the simulation will appear smoother and more realistic.

### The Challenge
We would ideally like to split up work among processing units to improve speedup. However, the calculations will likely involve memory accesses across units, so we will have to figure out how to distribute work to minimize communication overhead. We will consider using a blocked work distribution, since a pixel will have the largest effect on the pixels physically around it.

In addition, we have to consider how to represent our fluid data. Throughout the simulation, we store various quantities (velocity, pressure, fluid concentration, etc.) at various points in space. We will have to figure out the best way to represent that information in a way that takes advantage of cache locality and shared memory.


### Resources
Our project is inspired by this 
[simulator](https://haxiomic.github.io/GPU-Fluid-Experiments/html5/?q=UltraHigh).

[This paper](https://cg.informatik.uni-freiburg.de/intern/seminar/gridFluids_fluid-EulerParticle.pdf) is useful for understanding the physics behind the fluid simulation. 
[This webpage](http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html) also describes how a fluid simulator would be implemented.

### Goals and Deliverables
What we plan to achieve:
- Figure out how to use OpenGL to display the simulation (using Assignment 2 as a guide)
- Implement a baseline fluid simulator with a single thread of control (use the papers to learn and understand the physics of fluid motion)
- Transform the sequential algorithm to run with CUDA on an NVIDIA GPU (similar to Assignment 2)
- Distribute work among cores by analyzing and identifying workload characteristics of the algorithm.
- Try to further improve performance by considering other techniques like shared memory and cache usage. 

Our final deliverable will be a webpage featuring a widget that will allow users to interact with the fluid using their mouse. Our implemention has to perform well enough for the fluid to flow seamlessly. The baseline CPU implementation will undoubtedly be choppy, so our parallel implementation has to obtain enough speedup and perform well enough for the fluid to flow seamlessly.

We will demo this widget and show speedup graphs during our presentation.

### Platform Choice
We will will run the parallel version of our algorithm of the GHC machines with the NVIDIA GPUs using CUDA. 

### Schedule 
April 13: Start implementing baseline simulator on CPU

April 19: Finish implementing baseline simulator on CPU

April 22 - 25: Translate parallelizable computations of the algorithm into CUDA

April 26 - May 3: Determine optimal distribution of work among cores

May 4 - 9: Further improve performance to achieve acceptable speedup

May 10 - 11: Implement demo widget, Prepare for presentation

