# HW2 设计手册

**题目要求：**

设计基于极坐标的控制器实现叉车定点控制。

叉车运动学
$$
\dot\xi_R = \left[
\begin{matrix}
\dot x_R\\
\dot y_R\\
\dot \theta_R
\end{matrix}
\right]=
r\dot\phi
\left[
\begin{matrix}
\sin\beta\\
0\\
-\frac{\cos\beta}{L_2}
\end{matrix}
\right]
$$
## 基于极坐标的控制器设计

叉车系统的输入为
$$
u = \left[
\begin{matrix}
v_1\\
b
\end{matrix}
\right] \in \R^2
$$
分别表示主动轮的速度和角度,根据运动学模型可以将输入表示为机器人的速度$v_1$和角速度$w_1$
$$
\dot\xi_R = \left[
\begin{matrix}
\dot x_R\\
\dot y_R\\
\dot \theta_R
\end{matrix}
\right]=
R(\theta)\dot\xi_I=

\left[
\begin{matrix}
rv_1\sin b\\
0\\
-\frac{\cos b}{L_2}rv_1
\end{matrix}
\right]\\
\implies
\dot\xi_I = R(\theta)^{-1}\left[
\begin{matrix}
rv_1\sin b\\
0\\
-\frac{\cos b}{L_2}rv_1
\end{matrix}
\right]\\


\left[
\begin{matrix}
\dot x\\\dot y \\ \dot\theta
\end{matrix}
\right]
=rv_1\left[
\begin{matrix}
\cos\theta & -\sin\theta & 0\\
\sin\theta & \cos\theta & 0\\
0 & 0 & 1 
\end{matrix}
\right]\left[
\begin{matrix}
\sin b\\
0\\
-\frac{\cos b}{L_2}
\end{matrix}
\right]
$$

$$
\begin{cases}
\dot x_I =  rv_1\cos\theta  \sin b\\
\dot y_I = rv_1\sin\theta\sin b\\
\dot \theta_I = -rv_1\frac{\cos b}{L_2}
\end{cases}
$$

其中$\theta$表示机器人在世界坐标系下的朝向，$b$表示车轮在车身坐标系的角度，向前的时候为$\frac{\pi}{2}$。

为了方便表示，我们令$\gamma =b + \frac{\pi}{2}$来表示车轮与车身方向的夹角.
$$
\begin{cases}
\dot x_I =  rv_1\cos\theta  \sin (\gamma - \frac{\pi}{2})\\
\dot y_I = rv_1\sin\theta\sin (\gamma - \frac{\pi}{2})\\
\dot \theta_I = -rv_1\frac{\cos \gamma - \frac{\pi}{2})}{L_2}
\end{cases}
$$
考虑到极坐标系的坐标变换函数
$$
\rho = \sqrt{\Delta x^2 + \Delta y^2}\\
\alpha = \arctan2(\Delta y, \Delta x)\\
\beta =  -\theta- \alpha
$$
求导得
$$
\begin{cases}
\dot \rho = \frac{\dot{\Delta x}\Delta x+\dot{\Delta y}\Delta y}{\rho}\\
\dot \beta = -\frac{{\Delta x}}{\rho}(\frac{\dot y}{\Delta x} - \frac{\dot x \Delta y}{\Delta x^2}) \\
\dot \alpha = -\dot\beta-\dot\theta  \\
\end{cases}
,
\begin{cases}
\dot x_I =  rv_1\cos(\theta)  \sin b\\
\dot y_I =  rv_1\sin\theta\sin b\\
\dot \theta_I = -rv_1\frac{\cos b}{L_2}
\end{cases}
$$
代入运动模型得
$$
\dot\rho = rv_1\cos(\theta)  \sin b (-\cos(-\beta)) + (-\sin(\beta))rv_1\sin\theta\sin b\\
\dot\beta = (\cos(-\beta))[(-\cos(-\beta))rv_1\sin\theta\sin b +\sin(-\beta)
 rv_1\cos(\theta)  \sin b
]\\
\dot\alpha = -\dot\beta +rv_1\frac{\cos b}{L_2}
$$
设计线性控制器
$$
\begin{cases}
v_1 = k_\rho \rho\\
b = k_\alpha \alpha + k_\beta \beta
\end{cases}
$$
假设车轮的运动角度满足
$$
b\in[0,\pi]
$$
代入系统闭环误差模型
$$
\begin{cases}
\dot\rho = rk_\rho\rho\cos(\theta)  \sin ( k_\alpha \alpha + k_\beta \beta) (-\cos(-\beta)) + (-\sin(\beta))rk_\rho\rho\sin\theta\sin( k_\alpha \alpha + k_\beta \beta)\\
\dot\beta = (\cos(-\beta))[(-\cos(-\beta))rk_\rho\rho\sin\theta\sin ( k_\alpha \alpha + k_\beta \beta) +\sin(-\beta)
 rk_\rho\rho\cos(\theta)  \sin ( k_\alpha \alpha + k_\beta \beta)
]\\
\dot\alpha = -\dot\beta +rk_\rho\rho\frac{\cos ( k_\alpha \alpha + k_\beta \beta)}{L_2}
\end{cases}
$$
随后进行稳定性分析。设李雅普诺夫函数
$$
\begin{aligned}
V_1 &= \frac{1}{2}\rho^2 +\frac{1}{2}\alpha^2 + \frac{1}{2}\beta^2\\
\dot V_1 &= \rho\dot\rho + \alpha\dot\alpha + \beta\dot\beta\\


\end{aligned}
$$

<img src="/Users/skymac/Desktop/HW2 设计手册.assets/image-20201128135807890.png" alt="image-20201128135807890" style="zoom:40%;" />

其中$\alpha$表示移动机器人当前朝向和目标位置之间的误差，$\theta$表示移动机器人当前朝向在世界坐标系下的角度
$$
\rho = \sqrt{\Delta x^2 + \Delta y^2}\\
\alpha = \arctan2(\Delta y, \Delta x)\\
\beta = -\theta - \alpha
$$
于是得到在新的坐标系中移动机器人的一个系统描述
$$
\begin{cases}
\dot \rho = \frac{\dot{\Delta x}\Delta x+\dot{\Delta y}\Delta y}{\rho}\\
\dot \beta = -\frac{{\Delta x}^2}{\rho^2}(\frac{\dot y}{\Delta x} - \frac{\dot x \Delta y}{\Delta x^2}) \\
\dot \alpha = -\dot\beta-\dot\theta  \\
\end{cases}
,
\begin{cases}
\dot x_I =  r\dot\phi\cos(\theta)  \sin b\\
\dot y_I =  r\dot\phi\sin\theta\sin b\\
\dot \theta_I = -r\dot\phi\frac{\cos b}{L_2}
\end{cases}
$$

$$
k_1 = -\cos(-\beta ),k_2 = -\sin(-\beta)
\\
\begin{cases}
\dot \rho =\dot x k_1 + \dot y k_2\\
\dot \beta = \frac{-1}{\rho}[\dot y k_1 -\dot xk_2]\\
\dot \alpha =  -\dot\beta - \dot\theta
\end{cases}
$$

稳定性分析

$$