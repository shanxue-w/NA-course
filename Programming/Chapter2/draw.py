import numpy as np
import matplotlib.pyplot as plt
import os
import typing

def drawB(filename: str):
    """
    绘制文件中的数据。
    filename ./data/xxx.txt
    output ./figures/xxx.png
    """
    data = np.loadtxt(filename, delimiter=',')
    x = data[:, 0]
    y = data[:, 1]
    plt.plot(x, y)
    x1 = np.linspace(-5, 5, 1000)
    y1 = 1 / (1 + x1 ** 2)
    plt.plot(x1, y1, 'r--', label='$\\frac{1}{1+x^2}$')
    plt.legend()
    plt.savefig('./figures/' + filename.split('/')[-1].replace('.txt', '.png'))
    plt.close()

def drawC(filename: str):
    """
    绘制文件中的数据。
    filename ./data/xxx.txt
    output ./figures/xxx.png
    """
    data = np.loadtxt(filename, delimiter=',')
    x = data[:, 0]
    y = data[:, 1]
    plt.plot(x, y)
    x1 = np.linspace(-1, 1, 1000)
    y1 = 1 / (1 + 25 * x1 ** 2)
    plt.plot(x1, y1, 'r--', label='$\\frac{1}{1+25x^2}$')
    plt.legend()
    plt.savefig('./figures/' + filename.split('/')[-1].replace('.txt', '.png'))
    plt.close()

def draw(filename: str):
    """
    绘制文件中的数据。
    filename ./data/xxx.txt
    output ./figures/xxx.png
    """
    data = np.loadtxt(filename, delimiter=',')
    x = data[:, 0]
    y = data[:, 1]
    plt.plot(x, y)
    plt.savefig('./figures/' + filename.split('/')[-1].replace('.txt', '.png'))
    plt.close()

"""
绘制 B 的图像。
"""
drawB("./data/B_Newton2.txt")
drawB("./data/B_Newton4.txt")
drawB("./data/B_Newton6.txt")
drawB("./data/B_Newton8.txt")

"""
绘制 C 的图像。
"""
drawC("./data/C_Newton5.txt")
drawC("./data/C_Newton10.txt")
drawC("./data/C_Newton15.txt")
drawC("./data/C_Newton20.txt")

"""
绘制 D 的图像。
"""
draw("./data/D_Hermite.txt")
draw("./data/D_Hermite_prime.txt")

"""
绘制 E 的图像。
"""
draw("./data/E_Newton1.txt")
draw("./data/E_Newton2.txt")

"""
绘制 F 的图像。
"""
draw("./data/F_Bezier_10.txt")
draw("./data/F_Bezier_40.txt")
draw("./data/F_Bezier_160.txt")