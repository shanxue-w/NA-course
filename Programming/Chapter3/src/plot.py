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
    plt.axis("equal")
    plt.savefig('./figure/' + filename.split('/')[-1].replace('.txt', '.png'))
    plt.close()
drawE('result/E_heart_B10.txt')
drawE('result/E_heart_B40.txt')
drawE('result/E_heart_B160.txt')
drawE('result/E_heart_PP10.txt')
drawE('result/E_heart_PP40.txt')
drawE('result/E_heart_PP160.txt')
drawE('result/E_heart_cum_B10.txt')
drawE('result/E_heart_cum_B40.txt')
drawE('result/E_heart_cum_B160.txt')
drawE('result/E_heart_cum_PP10.txt')
drawE('result/E_heart_cum_PP40.txt')
drawE('result/E_heart_cum_PP160.txt')
drawE('result/E_2_B160.txt')
drawE('result/E_2_PP160.txt')
drawE('result/E_2_cum_B160.txt')
drawE('result/E_2_cum_PP160.txt')



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
drawE3d('result/E_3_B160.txt')
drawE3d('result/E_3_PP160.txt')
drawE3d('result/E_3_Ball160.txt')
drawE3d('result/E_3_Ball_cum160.txt')
drawE3d('result/E_3_PPBall160.txt')
drawE3d('result/E_3_PPBall_cum160.txt')
drawE3d('result/E_3_BallProj160.txt')
drawE3d('result/E_3_BallProj_cum160.txt')
drawE3d('result/E_3_PPBallProj160.txt')
drawE3d('result/E_3_PPBallProj_cum160.txt')
drawE3d('result/EXACT160.txt')


def draw(filename:str):
    data = np.loadtxt(filename, delimiter=',')
    x = data[:,0]
    y = data[:,1]
    plt.plot(x, y)
    plt.axis("equal")
    plt.savefig(filename.split('/')[-1].replace('.txt', '.png'))
    plt.close()


def drawF(filename:str):
    data = np.loadtxt(filename, delimiter=',')
    x = data[:,0]
    length = len(data[0]) - 1
    i = 1
    while i*(i+1)//2 < length:
        i += 1
    fig, axs = plt.subplots(i, i, figsize=(i*3, i*3))
    
    idx = 1
    for j in range(i):
        for k in range(i-j):
            y = data[:,idx]
            idx += 1
            ax = axs[k+j, j]
            ax.plot(x, y)
    # 隐藏多余的子图
    for j in range(i):
        for k in range(j):
            axs[k, j].axis('off')

    # 调整布局
    plt.tight_layout()
    plt.savefig('./figure/' + filename.split('/')[-1].replace('.txt', '.png'))


drawF('result/F_1.txt')
drawF('result/F_2.txt')

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