# 移动机器人导论HW3

> 任云帆
>
> SZ170420117

# 理论部分

叉车的运动模型可以表示为
$$
\dot\xi_I =R(\theta)^{-1} \dot\xi_R = r\dot\phi\left[\begin{matrix}
\cos\theta & -\sin\theta & 0\\
\sin\theta & \cos\theta & 0 \\
0 & 0 & 1

\end{matrix}\right]\left[\begin{matrix}\sin\beta \\ 0 \\ -\frac{\cos\beta}{L_2}\end{matrix}\right]\\
\implies \left[\begin{matrix}\dot x_I \\ \dot y_I\\\dot\theta_I\end{matrix}\right]
=r\dot\phi
\left[\begin{matrix}
\cos\theta\sin\beta\\
\sin\theta\sin\beta\\
-\frac{\cos\beta}{L_2}
\end{matrix}\right]
$$
对于平面移动的机器人，位姿可以用以下向量来表示
$$
p = \left[\begin{matrix}x\\y\\\theta\end{matrix}\right]
$$
对于叉车机器人来说，其未知可以从一个已知位置开始，并将运动进行积分来进行估计。假设时间增量内叉车的输入变化为$\Delta\beta,\Delta s$，那么行走的增量可以表示为
$$
\Delta \theta = -\Delta s \frac{\cos(\beta+\Delta\beta/2)}{L_2} \\
\Delta x = \Delta s\sin(\beta+\Delta\beta/2)\cos(\theta+\Delta\theta/2)\\
\Delta y = \Delta s\sin(\beta+\Delta\beta/2)\sin(\theta+\Delta\theta/2)\\
$$
由此我们可以得到更新过的位置$p'$
$$
p' =f(x,y,\theta,\Delta s, \beta)
=\left[\begin{matrix}x \\ y \\ \theta\end{matrix}\right]+
\left[\begin{matrix}

 \Delta s\sin(\beta+\Delta\beta/2)\cos(\theta-\Delta s \frac{\cos(\beta+\Delta\beta/2)}{2L_2} )\\
\Delta s\sin(\beta+\Delta\beta/2)\sin(\theta-\Delta s \frac{\cos(\beta+\Delta\beta/2)}{2L_2} )\\
-\Delta s \frac{\cos(\beta+\Delta\beta/2)}{L_2} \\
\end{matrix}\right]
$$
接下来建立$p'$的误差模型，已得到里程表位置估计的协方差矩阵$\Sigma_{p'}$。对于运动增量$(\Delta s;\Delta\beta)$，我们假定协方差矩阵$\Sigma_\Delta$
$$
\Sigma_\Delta = covar(\Delta s;\Delta \beta ) = 
\left[\begin{matrix}
k_s|\Delta s|& 0 \\
0&k_b
\end{matrix}\right]
$$
其中$\Delta s$是驱动轮运动的距离，由于转向角度大小与协方差无关，因此协方差矩阵设为一常数，$k_s,k_b$分别为误差常数。我们做出以下假设

* 驱动轮的转动和角度变换是独立的
* 误差的方差正比于各自运动的角度。

假定位置$p$和系统输入$\Delta _{sb} = (\Delta s;\Delta \beta)$不相关。我们用一阶泰勒展开近似函数$f$的微分，我们计算两个雅克比矩阵，$F_p = \nabla_p f$ 和$F_{sb} = \nabla_{\Delta sb}f$
$$
F_p = \nabla_p f = \nabla_p(f^T) = \left[\begin{matrix}
\frac{\partial f}{\partial x} & \frac{\partial f}{\partial y} & \frac{\partial f}{\partial \theta}
\end{matrix}\right]
 = \left[\begin{matrix}
 1 & 0 & -\Delta s\sin(\beta + \Delta\beta / 2)\sin(\theta - \Delta s\frac{\cos(\beta + \Delta\beta/2)}{2L_2})\\
 0 & 1 & \Delta s\sin(\beta + \Delta\beta / 2)\cos(\theta - \Delta s\frac{\cos(\beta + \Delta\beta/2)}{2L_2})\\
 0 & 0 & 1
 \end{matrix}\right]
$$
设$b = \beta + \Delta \beta / 2$
$$
F_{sb = }\left(\begin{array}{cc} 
\cos\left(\theta -\frac{\Delta s\,\cos\left(b \right)}{2\,L_{2}}\right)\,\sin\left(b \right)+\frac{\Delta s\,\sin\left(\theta -\frac{\Delta s\,\cos\left(b \right)}{2\,L_{2}}\right)\,\cos\left(b \right)\,\sin\left(b \right)}{2\,L_{2}} & \Delta s\,\cos\left(\theta -\frac{\Delta s\,\cos\left(b \right)}{2\,L_{2}}\right)\,\cos\left(b \right)-\frac{{\Delta s}^2\,\sin\left(\theta -\frac{\Delta s\,\cos\left(b \right)}{2\,L_{2}}\right)\,{\sin\left(b \right)}^2}{2\,L_{2}}\\ \sin\left(\theta -\frac{\Delta s\,\cos\left(b \right)}{2\,L_{2}}\right)\,\sin\left(b \right)-\frac{\Delta s\,\cos\left(\theta -\frac{\Delta s\,\cos\left(b \right)}{2\,L_{2}}\right)\,\cos\left(b \right)\,\sin\left(b \right)}{2\,L_{2}} & \Delta s\,\sin\left(\theta -\frac{\Delta s\,\cos\left(b \right)}{2\,L_{2}}\right)\,\cos\left(b \right)+\frac{{\Delta s}^2\,\cos\left(\theta -\frac{\Delta s\,\cos\left(b \right)}{2\,L_{2}}\right)\,{\sin\left(b \right)}^2}{2\,L_{2}}\\ -\frac{\cos\left(b \right)}{L_{2}} & \frac{\Delta s\,\sin\left(b \right)}{L_{2}} \end{array}\right)
$$

因此可以由误差传播定律推出误差协方差矩阵的更新函数。
$$
\Sigma_{p'} = \nabla_{p}f\Sigma_p\nabla_pf^T + \nabla_{\Delta sb}f\Sigma_\Delta\nabla_{\Delta s b}f^T
$$

# 实验部分

在完成$\Sigma_{p'}$的递推公式后，我们得到叉车运动的协方差矩阵$\Sigma(t)$。对每一时刻的协方差矩阵进行特征分解，得到$x,y$方向上的特征值和特征向量，将特征向量归一化，作为方向，特征值大小作为长度绘制椭圆。再结合Hw2中的小车运动，基于Python3.7实现了小车运动和误差传递的仿真动画。

