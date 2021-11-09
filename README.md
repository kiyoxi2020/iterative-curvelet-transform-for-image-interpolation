# iterative-curvelet-transform-for-image-interpolation

```
Run demo_ISTA_for_l1_my.m
```

# 问题描述

- L1优化问题基本定义如下：
    
    $$ 
    P_\epsilon :\left\{ \begin{matrix}&\tilde{x}=\arg\min_x \Vert x\Vert_1~~s.t.\Vert Ax-y\Vert_2\leq \epsilon \\ &\tilde{y}=A\tilde{x} \end{matrix} \right.
    $$
    

# 迭代阈值收缩算法求解（Iterative shrinkage threshold algorithm, ISTA）

- 将上述问题转化为一系列无约束优化问题
    
    $$
    P_\lambda:\left\{ \begin{matrix}\tilde{x}_\lambda=\arg\min_x\Vert y-Ax\Vert_2^2+\lambda\Vert x\Vert_1 \\ 
    \tilde{y}_\lambda=A\tilde{x}_\lambda \end{matrix} \right.
    $$
    
    可以通过不断减小$\lambda$使得$P_\lambda$的解逼近$P_\epsilon$的解，$\lambda$的初始值满足
    
    $$⁍$$
    
- 迭代优化式为：
    
    $$
    x=T_\lambda \left(x+A^T(y-Ax)\right)
    $$
    
    其中，$T_\lambda$表示软阈值函数（soft-thresholding function）：
    
    $$
    T_\lambda(x):=\text{sgn}(x)\cdot \max(0, |x|-|\lambda|)
    $$
    
    $\lambda$在初始时较大（只重构少数几个重要的部分），随着迭代的进行，$\lambda$取值逐渐减小，即逐渐开始重构次要的部分
