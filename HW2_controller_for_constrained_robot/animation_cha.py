import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from math import *


def cartesian2polar(x, y, theta):
    rho = sqrt(x*x + y*y)
    beta = atan2(y,x)
    alpha = theta -beta
    return rho, beta, alpha

def polar2cartesian(rho, beta, alpha, cur_pos, tar_pos):
    x_I = rho * cos( beta)
    y_I = rho * sin( beta)
    theta_I =  (beta + alpha)
    XI = np.array([x_I, y_I,theta_I])
    theta_tar = target_pos[2]
    R = np.array([[cos(theta_tar), sin(theta_tar), 0],
                  [-sin(theta_tar), cos(theta_tar), 0],
                  [0, 0, 1]
                  ])
    x_G = np.linalg.inv(R) @ XI - target_pos
    # x_G = (R) @ (XI)
    # x_G = XI
    return 0-x_G[0], 0-x_G[1], x_G[2]
    # return x_G[0], x_G[1], x_G[2]

def rad2deg(rad):
    return rad/pi*180

def deg2rad(deg):
    return deg/180*pi

def rotation(x,theta_tar):
    if(x.sum() == 0):
        return x
    x_ = np.ones(3)
    if(x.shape[0] != 3):
        x_[:2] = x
    else:
        x_ = x

    R = np.array([[cos(theta_tar), sin(theta_tar), 0],
                  [-sin(theta_tar), cos(theta_tar), 0],
                  [0, 0, 1]
                  ])
    return R @ x_

def moveSteChache(v1, b ,rho,beta,alpha,dt):
    r = 0.1
    l = 0.3
    theta_ = - alpha - beta
    k1 = -cos(-beta)
    k2 = -sin(-beta)
    dot_x = r * v1 * cos(theta_) * sin(b)
    dot_y = r * v1 * sin(theta_) * sin(b)
    dot_t = r * v1 * (-cos(b)/l)

    dot_rho =dot_x * k1 + dot_y * k2
    dot_beta = (-1/rho) * (dot_y*k1 - dot_x*k2)
    dot_alpha = -dot_beta - dot_t

    rho = rho + dt * dot_rho
    beta = beta + dt *dot_beta
    alpha = alpha + dt * dot_alpha
    return rho, beta, alpha

def limiting(theta, lower, upper):
    if(lower>upper):
        t = lower
        lower = upper
        upper = t
    if(theta > upper):
        theta = upper
    elif(theta<lower):
        theta = lower
    else:
        theta = theta

    return theta

def runOnceCha(current_pos,target_pos):
    global xs,ys,thetas,bs
    tile =rotation(target_pos - current_pos, target_pos[2])

    kp = 2.5
    ka = 4
    kb = -2
    dt = 0.1

    x = current_pos[0]
    y = current_pos[1]
    theta = current_pos[2]
    x_tar = target_pos[0]
    y_tar = target_pos[1]
    theta_tar = target_pos[2]
    rho, beta, alpha = cartesian2polar(tile[0], tile[1], current_pos[2] - target_pos[2])
    # delta = target_pos - current_pos
    # rho, beta, alpha = coorTrans(delta[0], delta[1], current_pos[2]-target_pos[2])
    print(rho, rad2deg(beta), rad2deg(alpha),rad2deg(theta))

    while (rho > 0.01):# or abs(beta) > 0.01 or abs(alpha) > 0.01):
        """
        控制器
        """
        # v1 = kp * rho  * cos(alpha)
        # w = ka * alpha +  (alpha - kb * beta) * kp * sin(alpha) * cos(alpha) / alpha
        b = pi/2 + ka * alpha + kb * (beta) # (alpha - kb * beta) * kp * sin(alpha) * cos(alpha) / alpha
        b = limiting(b, deg2rad(90-80), deg2rad(90+80))
        v1 = kp * rho
        """
        控制步长
        """
        rho, beta, alpha = moveSteChache(v1, b, rho, beta, alpha, dt)
        # rho, beta, alpha = moveStep(v1, w, rho, beta, alpha, dt)
        """
        可视化
        """
        x_, y_, theta_ = polar2cartesian(rho, beta, alpha, current_pos, target_pos)
        # print(rho, rad2deg(beta), rad2deg(alpha), rad2deg(theta_), rad2deg(b))
        xs.append(x_)
        ys.append(y_)
        bs.append(b)
        thetas.append(theta_ - target_pos[2])
    xs = np.array(xs)
    ys = np.array(ys)
    bs = np.array(bs)
    thetas = np.array(thetas)
    return xs, ys

def plotRect(pos,recSize,theta):
    pos = np.array(pos)
    recSize = np.array(recSize)
    ps = np.zeros((4,2))
    ps[0,:] = np.array([-recSize[0]/2,recSize[1]/2])
    ps[1,:] = np.array([recSize[0]/2,recSize[1]/2])
    ps[2,:] = np.array([recSize[0]/2,-recSize[1]/2])
    ps[3,:] = np.array([-recSize[0]/2,-recSize[1]/2])

    for i in range(0,4):
        ps[i,:] = rotation(ps[i,:],(theta))[:2]
    lines = np.zeros((4,2,100))
    ps[:, 0] = ps[:, 0] + pos[0]
    ps[:, 1] = ps[:, 1] + pos[1]
    # p0-p1
    lines[0,0,:] = np.linspace(ps[0,0], ps[1,0],100)
    lines[0,1,:] = np.linspace(ps[0,1], ps[1,1],100)
    # p0-p2
    lines[1,0,:] = np.linspace(ps[1,0], ps[2,0],100)
    lines[1,1,:] = np.linspace(ps[1,1], ps[2,1],100)
    # p2-p3
    lines[2,0,:] = np.linspace(ps[2,0], ps[3,0],100)
    lines[2,1,:] = np.linspace(ps[2,1], ps[3,1],100)
    # p1-p3
    lines[3,0,:] = np.linspace(ps[0,0], ps[3,0],100)
    lines[3,1,:] = np.linspace(ps[0,1], ps[3,1],100)
    # plt.axis('equal')
    # for i in range (0, 4):
    #     plt.plot(lines[i,0], lines[i,1], color='r')
    return lines

def plotCha(pos,theta,b):
    b =  deg2rad(90) - b
    theta = deg2rad(90) - theta# +deg2rad(180)
    pos = np.array(pos)
    width = 0.5
    length = 1
    posWheel = np.zeros(2)
    posWheel[0] = pos[0]+sin(theta) * 0.3
    posWheel[1] = pos[1]+cos(theta) * 0.3
    lines1 = plotRect(pos,[width,length],theta)
    lines2 = plotRect(posWheel,[0.1,0.2], theta+b)
    return lines1, lines2
    # for i in range(0, 4):
    #     plt.axis('equal')
    #     plt.plot(lines1[i, 0], lines1[i, 1], color='r')
    #     plt.plot(lines2[i, 0], lines2[i, 1], color='b')
    # plt.axis('equal')
    #     # for i in range (0, 4):
    #     #     plt.plot(lines[i,0], lines[i,1], color='r')



fig = plt.figure(figsize=(10,10),tight_layout=True)
xs = []
ys = []
thetas=[]
bs = []
xdata, ydata = [], []
plt.xlim(-5,5)
plt.ylim(-5,5)
ln1, = plt.plot([], [], 'r')
ln2, = plt.plot([], [], 'r')
ln3, = plt.plot([], [], 'r')
ln4, = plt.plot([], [], 'r')
ln5, = plt.plot([], [], 'b')
ln6, = plt.plot([], [], 'b')
ln7, = plt.plot([], [], 'b')
ln8, = plt.plot([], [], 'b')


def update_points(frame):
    global xs,ys,thetas

    xdata.append(xs[frame])
    ydata.append(ys[frame])
    lines1, lines2 = plotCha([xs[frame],ys[frame]],thetas[frame],-bs[frame])
    ln.set_data  (xdata,ydata)

    ln1.set_data(lines1[0, 0], lines1[0, 1])
    ln2.set_data(lines1[1, 0], lines1[1, 1])
    ln3.set_data(lines1[2, 0], lines1[2, 1])
    ln4.set_data(lines1[3, 0], lines1[3, 1])

    ln5.set_data(lines2[0, 0], lines2[0, 1])
    ln6.set_data(lines2[1, 0], lines2[1, 1])
    ln7.set_data(lines2[2, 0], lines2[2, 1])
    ln8.set_data(lines2[3, 0], lines2[3, 1])
    return ln, ln1, ln2, ln3, ln4, ln5, ln6, ln7, ln8,


if(__name__ == "__main__"):
    current_pos = np.array([4, 0, deg2rad(0)])
    target_pos =  np.array([0, 0, deg2rad(0)])

    runOnceCha(current_pos, target_pos)
    # plt.plot(xs, ys)
    # # plt.scatter(xs[:10], ys[:10])
    # plotCha([0,0],deg2rad(30),deg2rad(0))
    # plt.axis('equal')
    # plt.show()
    ln, = plt.plot(xs[0], ys[0], "g")
    plt.grid(ls="--")
    # 开始制作动画
    ani = animation.FuncAnimation(fig, update_points, xs.shape[0], interval=20, blit=True)
    # animation.PillowWriter('cur0_tar60.gif', writer='imagemagick', fps=10)
    ani.save('cur0_tar620.gif', writer='imagemagick', fps=100)
    plt.show()