# MATLAB Programming for Virtual Element Methods

If you have any problems, please contact me via the email: terenceyuyue@sjtu.edu.cn

mVEM is a MATLAB software package containing efficient and easy-following codes for various virtual element methods (VEMs) published in the
literature. 

--------------
## MESH

- We present some basic functions to show the polygonal meshes, including marking of the nodes, elements and (boundary) edges.

- For the convenience of computation, some auxiliary mesh data are introduced. 
  The idea stems from the treatment of triangulation in iFEM, which is generalized to polygonal meshes with certain modifications. 

- We also provide a boundary setting function to identify the Neumann and Dirichlet boundaries (See setboundary.m).

- The mesh generation algorithms are also described in detail.

- We present an efficient implementation of the mesh refinement for polygonal meshes. To the best of our knowledge, this is the first publicly available implementation of the polygonal mesh refinement algorithms. We remove the small edges by requiring the one-hanging-node rule.
- nonConvexMesh: generates a typical nonconvex polygonal mesh
- gunelveMesh: generates the Gunelve mesh used in the Cook's membrane problem
- dualMesh: Generate polygonal mesh by establishing the dual mesh of a given Delaunay triangulation.

-------------------------
## Conforming Virtual Element Methods

### Poisson equation (k <= 3)

We describe the details of designing the codes of conforming VEMs for Poisson equation with k up to 3, 
including the computation of elliptic projection and L2 projection, the matrix form of the approximated variational problems as well as the treatment of boundary conditions.
We also show how the errors under L2,H1 and energy norms can be computed using the numerical d.o.f.s.

### Linear elasticity problems (lowest order k = 1)

The VEM of lowest order for linear elasticity problems is introduced for both the displacement-type (or Navier-type) and tensor-type bilinear forms. 

### Plate bending problems (lowest order k = 2)

Three VEMs involved in the literature are described in detail, i.e., the C1, C0 and Morley-type elements.

### Fourth-order singular perturbation problems

 This problem combines the techniques in Poisson equation and plate bending problems.
 
 


---------------
## Nonconforming Virtual Element Methods

### Poisson equation

   - Standard nonconforming VEM in the lowest order case
   - The nonconforming VEM with a continuous treatment on the boundary of the domain

### Linear elasticity

For linear elasticity problems, we mainly consider the (locking-free) virtual elements in the lowest order case. 
   
   - conforming and nonconforming VEMs for the Navier/displacement form  
   - conforming and nonconforming VEMs for the tensor form 
   - nonconforming VEMs for a continuous treatment of the domain boundary
   - locking-free VEMs: Kouhia-Stenberg VEM (Vnc \times Vc),  reduced integration technique based VEMs
   - a standard benchmark: Cook's membrane problem


------------------
## Mixed Virtual Element Methods

 - Mixed vem for Darcy problem in the lowest order H(div)-conforming virtual element
 
 - Mixed vem for Darcy problem in the lowest order lifting H(div)-conforming virtual element
 
 - Divergence-free mixed vem for Stokes problem in the lowest order


------------------
##  Adaptive Virtual Element Methods

  - Mesh refinement of polygonal meshes
  
  - Adaptive VEM for Poisson equation in the lowest order


------------------
##  Variational Inequalities

   We provide the implementation of the the virtual element method for a simplified friction problem, which is a typical elliptic variational inequality of the second kind. The discrete problem is solved by using the Uzawa iteration algorithm, hence the implentation is reduced to that of Poisson equation in each iteration.


------------------
##  3-D Linear Virtual Element Methods

   We present a simple and efficient implementation of the linear virtual element method for the three dimensional Poisson equation. 
The design ideas can be directly generalized to higher-order cases.


-------------------
## Vectorized Implementation

- The vectorized implementatioin for conforming VEMs of Poisson equation has been accomplished in a unified way. 
     
- The current realization may be less straighforward compared to the FEMs.

- Only the calculation of projection matrices is vectorized, and the assembly procedure is still implemented by looping over the elements


-------------------
## Acknowlegement

The author would like to thank Sen Lin and Fang Feng, my schoolmates in Shanghai Jiao Tong University for many helpful discussions.

The author thanks Postdoc Dohyun for pointing out the bugs in PolyMeshRefine.m and provide another version PolyMeshRefine1.m. 
Detailed credits can be found in the M-lint of several m files. 

-----------------
## Citation

If you feel it is helpful for your research, please acknowledge your use by citing:

Y. Yu, mVEM: MATLAB Programming for Virtual Element Methods, https://github.com/Terenceyuyue/mVEM, Shanghai Jiao Tong University, 2019-2022.

    @misc{mvem,
          title = { {mVEM: MATLAB Programming for Virtual Element Methods} },
          url = {https://github.com/Terenceyuyue/mVEM},
          author = {Yu, Y.}
          institution = {Shanghai Jiao Tong University},
          year = {2019-2022},
          }


  PS: 

      程序对区域、边界条件等的处理都较为一般，可容易地移植到其他问题。
      
      编程过程完全遵循有限元的思路。
      
      相比网上有的程序，过程更加清晰，充分使用了 MATLAB 的向量运算。
      
      在计算各种投影的边界项时，部分程序中使用有限元载荷向量的装配技巧避免了繁琐的人为处理。

  TO DO:
  
      - Vectorized implementation
      
      - Time-dependent problems and nonlinear problems
      
      - Various applications: Cahn–Hilliard problem, Stokes–Darcy problem, Navier-Stokes, etc.
      
      - Mesh generation in 3-D 

