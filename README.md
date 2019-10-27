# MATLAB Programming for Virtual Element Methods


## Mesh 
 Chapter 1 presents some basic functions to show the polygonal meshes, including marking of the nodes, elements and (boundary) edges.
 
 For the convenience of computation, we introduce some auxiliary mesh data. The idea stems from the treatment of triangulation in iFEM,     which is generalized to polygonal meshes with certain modifications. 
 
 We also provide a boundary setting function to identify the Neumann and Dirichlet boundaries (See setboundary.m).
 
## The VEM of lowest order for linear elasticity problems is introduced in Chapter 3 for both the displacement-type and tensor-type bilinear forms. The sparse assembly indices for vector equation are presented 

## Poisson VEM (lowest order)
The virtual element method (VEM) of lowest order is introduced in Chapter 2, including the computation of elliptic projection and L2 projection, the matrix form of the approximated variational problems as well as the treatment of boundary conditions. 
 
  PS: 
  
      程序对区域、边界条件等的处理都较为一般，可容易地移植到其他问题。
      
      编程过程完全遵循有限元的思路。
      
      相比网上有的程序，过程更加清晰，充分使用了 MATLAB 的向量运算。
      
  Undo:
 
      - 线弹性边值问题的一阶虚拟元方法
      
      - Poisson 方程的二阶虚拟元方法

