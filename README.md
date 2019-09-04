# MATLAB Programming for Virtual Element Methods


## Mesh 
 Chapter 1 presents some basic functions to show the polygonal meshes, including marking of the nodes, elements and (boundary) edges.
 
 For the convenience of computation, we introduce some auxiliary mesh data. The idea stems from the treatment of triangulation in iFEM,     which is generalized to polygonal meshes with certain modifications. 
 
 We also provide a boundary setting function to identify the Neumann and Dirichlet boundaries (See setboundary.m).
 

## Poisson VEM (lowest order)
The virtual element method (VEM) of lowest order for Poisson equation is introduced in Chapter 2, including the computation of elliptic projection and L2 projection, the matrix form of the approximate variational problems as well as the treatment of boundary conditions. 

## Linear elasticity VEM (loweset order)
第三章介绍线弹性问题的虚拟元方法，读者可先阅读 mFEM 的相关章节。对向量型方程，我们给出了类似 Poisson 方程那里的稀疏装配指标，从而整个装配过程完全类似 Poisson 方程进行。线弹性问题有两种类型的变分形式，文中称为位移型 (displacement type) 和张量型 (tensor type)。对这两种形式，分别说明了椭圆投影以及 L2 投影的定义和计算。

PS: 

     - 网上也有弹性问题 VEM 方法的文章及程序，但与本文的可能不相同 (如 locking free vem，投影限制条件的处理)。
  
     - 可以这么说，本文首次给出了 VEM 方法的 iFEM 方式处理。

        
      
  Undo:
  
      - General boundary conditions for linear elasticity problems
  
      - Locking fee VEM for linear elasticity problems
      
      - Poisson 方程的二阶虚拟元方法

