# MATLAB Programming for Virtual Element Methods


## Mesh 
 Chapter 1 presents some basic functions to show the polygonal meshes, including marking of the nodes, elements and (boundary) edges.
 
 For the convenience of computation, we introduce some auxiliary mesh data. The idea stems from the treatment of triangulation in iFEM,     which is generalized to polygonal meshes with certain modifications. 
 
 We also provide a boundary setting function to identify the Neumann and Dirichlet boundaries (See setboundary.m).
 

## Poisson VEM (lowest order)
The virtual element method (VEM) of lowest order for Poisson equation is introduced in Chapter 2, including the computation of elliptic projection and L2 projection, the matrix form of the approximate variational problems as well as the treatment of boundary conditions. 

## Linear elasticity VEM (loweset order)
For linear elasticity problems, we consider two kinds of variational formulations, named displacement-type VEM and tensor-type VEM, respectively. One can refer to the relevant chapter given in the mFEM repository beforehand. We provide sparse indices for fast assembling stiffness matrix and load vector as given in the programming of Poisson equation. The only thing is to compute the elliptic projections and L2 projections with the help of the standard assembly algorithm. 


PS: 

     - 网上也有弹性问题 VEM 方法的文章及程序，但与本文的可能不相同 (如 locking free vem，投影限制条件的处理)。
  
     - 可以这么说，本文首次给出了 VEM 方法的 iFEM 方式处理。

        
      
  Undo:
  
      - General programming framework of virtual element methods       
       
      - General boundary conditions for linear elasticity problems
  
      - Locking fee VEM for linear elasticity problems
      
      - Poisson 方程的二阶虚拟元方法

