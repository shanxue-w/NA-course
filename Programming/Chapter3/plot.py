import matplotlib.pyplot as plt
import numpy as np

def drawA(filename:str):
    data = np.loadtxt(filename, delimiter=',')
    x = data[:,0]
    y = data[:,1]
    plt.plot(x, y)
    x1 = np.linspace(-1, 1, 500)
    y1 = 1.0/(1.0+25.0*x1*x1)
    plt.plot(x1, y1, 'r--', label='$f(x)=\\frac{1}{1+25x^2}$')
    plt.legend()
    plt.savefig('./figure/' + filename.split('/')[-1].replace('.txt', '.png'))
    plt.close()

drawA('result/A_6.txt')
drawA('result/A_11.txt')
drawA('result/A_21.txt')
drawA('result/A_41.txt')
drawA('result/A_81.txt')


def drawC(filename:str):
    data = np.loadtxt(filename, delimiter=',')
    x = data[:,0]
    y = data[:,1]
    plt.plot(x, y)
    x1 = np.linspace(-5, 5, 1000)
    y1 = 1.0/(1.0+x1*x1)
    plt.plot(x1, y1, 'r--', label='$f(x)=\\frac{1}{1+x^2}$')
    plt.legend()
    plt.savefig('./figure/' + filename.split('/')[-1].replace('.txt', '.png'))
    plt.close()

drawC('result/C_3.txt')
drawC('result/C_2.txt')

def drawE(filename:str):
    data = np.loadtxt(filename, delimiter=',')
    x = data[:,0]
    y = data[:,1]
    plt.plot(x, y)
    plt.savefig('./figure/' + filename.split('/')[-1].replace('.txt', '.png'))
    plt.close()
drawE('result/E_heart_B10.txt')
drawE('result/E_heart_B40.txt')
drawE('result/E_heart_B160.txt')
drawE('result/E_heart_B500.txt')
drawE('result/E_heart_PP10.txt')
drawE('result/E_heart_PP40.txt')
drawE('result/E_heart_PP160.txt')
drawE('result/E_heart_PP500.txt')
drawE('result/E_2_B100.txt')
drawE('result/E_2_PP100.txt')


# 1.000,2.0000,3.0000
# plot 3d data
def drawE3d(filename:str):
    data = np.loadtxt(filename, delimiter=',')
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z)
    plt.savefig('./figure/' + filename.split('/')[-1].replace('.txt', '.png'))
    plt.close()
drawE3d('result/E_3_B100.txt')
drawE3d('result/E_3_PP100.txt')


def draw(filename:str):
    data = np.loadtxt(filename, delimiter=',')
    x = data[:,0]
    y = data[:,1]
    plt.plot(x, y)
    plt.savefig(filename.split('/')[-1].replace('.txt', '.png'))
    plt.close()
# draw('log.txt')

# # read file log.txt
# with open('log.txt', 'r') as f:
#     lines = f.readlines()

# # split lines into two lists
# x = []
# y = []
# for line in lines:
#     line = line.strip()
#     if line:
#         x.append(float(line.split(',')[0]))
#         y.append(float(line.split(',')[1]))

# # plot the data
# plt.plot(x, y)
# plt.xlabel('x')
# plt.ylabel('y')
# plt.show()