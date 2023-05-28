# -*_ coding: utf-8 -*_

import math
import numpy as np

"""
Created on Mar 15 2023
@author: Yilong Wu（Uyoin：https://github.com/uyoin）

homework:后方交会法代码实现
"""

def rotation_matrix(INS):
    """计算旋转矩阵。

    根据给定的姿态角（omega, phi, kappa），计算旋转矩阵，用于描述相机坐标系和地面坐标系之间的转换关系。

    参数：
        INS (list): 一个包含三个元素的列表，分别表示INS测出的姿态角（omega, phi, kappa），单位为度。

    返回：
        M (numpy.ndarray): 一个3x3的二维数组，表示旋转矩阵。
    """
    omega = INS[0]
    phi = INS[1]
    kappa = INS[2]

    # 将角度转换为弧度
    omega = math.radians(omega)
    phi = math.radians(phi)
    kappa = math.radians(kappa)

    # 根据旋转矩阵的公式，计算各个元素的值
    m11 = math.cos(phi) * math.cos(kappa)
    m12 = math.sin(omega) * math.sin(phi) * math.cos(kappa) + math.cos(omega) * math.sin(kappa)
    m13 = -math.cos(omega) * math.sin(phi) * math.cos(kappa) + math.sin(omega) * math.sin(kappa)
    m21 = -math.cos(omega) * math.sin(kappa)
    m22 = -math.sin(omega) * math.sin(phi) * math.sin(kappa) + math.cos(omega) * math.cos(kappa)
    m23 = math.cos(omega) * math.sin(phi) * math.sin(kappa) + math.sin(omega) * math.cos(kappa)
    m31 = math.sin(phi)
    m32 = -math.sin(omega) * math.cos(phi)
    m33 = math.cos(omega) * math.cos(phi)

    # 将各个元素组合成一个矩阵，并返回
    M = np.array([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]])
    return M


def error_equations(GNSS, INS, M, P, G, f):
    """构建误差方程组。

      根据给定的GNSS坐标、INS欧拉角、旋转矩阵、像点坐标、地面控制点坐标和相机焦距，计算误差方程组中的误差向量epsilon和系数矩阵B。

      参数：
          GNSS (list): 一个包含三个元素的列表，分别表示GNSS测出的坐标（XL, YL, ZL），单位为米。
          INS (list): 一个包含三个元素的列表，分别表示INS测出的姿态角（omega, phi, kappa），单位为度。
          M (numpy.ndarray): 一个3x3的二维数组，表示旋转矩阵。
          P (list): 一个包含两个元素的列表，分别表示像点坐标（XB, YB），单位为毫米。
          G (list): 一个包含三个元素的列表，分别表示地面控制点坐标（XA, YA, ZA），单位为米。
          f (float): 相机焦距，单位为毫米。

      返回：
          B (numpy.ndarray): 一个数组，表示误差方程组中的系数矩阵B。
          epsilon (numpy.ndarray): 一个数组，表示误差方程组中的误差向量epsilon。
      """
    # 初始化
    epsilon = []
    B = []

    XL = GNSS[0]
    YL = GNSS[1]
    ZL = GNSS[2]

    omega = INS[0]
    phi = INS[1]
    kappa = INS[2]

    # 通过for循环遍历P,G各的坐标，来计算各点的系数矩阵和误差矩阵
    for P,G in zip(P.values(), G.values()):

        XA = G['X']
        YA = G['Y']
        ZA = G['Z']

        XB = P['x']
        YB = P['y']

        deX = G['X'] - XL
        deY = G['Y'] - YL
        deZ = G['Z'] - ZL

        # 计算r,s,q
        r = M[0][0] * (XA - XL) + M[0][1] * (YA - YL) + M[0][2] * (ZA - ZL)
        s = M[1][0] * (XA - XL) + M[1][1] * (YA - YL) + M[1][2] * (ZA - ZL)
        q = M[2][0] * (XA - XL) + M[2][1] * (YA - YL) + M[2][2] * (ZA - ZL)

        # 计算矩阵B的元素
        b11 = f / np.dot(q,q) * (r * (-M[2][2] * deY + M[2][1] * deZ) - q * (-M[0][2] * deY + M[0][1] * deZ))
        b12 = f / np.dot(q,q) * (r * (math.cos(phi) * deX) + (math.sin(omega) * math.sin(phi) * deY) - (
                math.cos(omega) * math.sin(phi) * deZ))
        b13 = -f / q * (M[1][0] * deX + M[1][1] * deY + M[1][2] * deZ)
        b14 = f / np.dot(q,q) * (r * M[2][0] - q * M[0][0])
        b15 = f / np.dot(q,q) * (r * M[2][1] - q * M[0][1])
        b16 = f / np.dot(q,q) * (r * M[2][2] - q * M[0][2])
        b21 = f / np.dot(q,q) * (s * (-M[2][2] * deY + M[2][1] * deZ) - q * (-M[1][2] * deY + M[1][1] * deZ))
        b22 = f / np.dot(q,q) * (s * (math.cos(phi) * deX) + (math.sin(omega) * math.sin(phi) * deY) - (
                math.cos(omega) * math.sin(phi) * deZ))
        b23 = f / q * (M[0][0] * deX + M[0][1] * deY + M[0][2] * deZ)
        b24 = f / np.dot(q,q) * (s * M[2][0] - q * M[1][0])
        b25 = f / np.dot(q,q) * (s * M[2][1] - q * M[1][1])
        b26 = f / np.dot(q,q) * (s * M[2][2] - q * M[1][2])

        # 计算矩阵epsilon的元素
        J = XB - 0 + f * (r / q)
        K = YB - 0 + f * (s / q)


        B.append([b11, b12, b13, b14, b15, b16])
        B.append([b21, b22, b23, b24, b25, b26])
        epsilon.append((np.array([J])))
        epsilon.append((np.array([K])))

    # 系数矩阵B，误差矩阵epsilon转换为numpy数组，并返回
    B = np.array(B)
    epsilon = np.array(epsilon)
    return B, epsilon


def solve_unknowns(B, epsilon):
    """求解未知数delta。

    给定误差方程组中的系数矩阵B和误差向量epsilon，计算未知数delta。

    参数：
        B (numpy.ndarray): 一个数组，表示误差方程组中的系数矩阵B。
        epsilon (numpy.ndarray): 一个数组，表示误差方程组中的误差向量epsilon。

    返回：
        delta (numpy.ndarray): 一个数组，表示误差方程组中的未知数delta。
    """
    delta = np.dot(np.linalg.inv(np.dot(B.T, B)), np.dot(B.T, epsilon))

    return delta


def correction(delta, INS, GNSS):
    """计算改正数和改正后的值。

    给定误差方程组中的未知数delta，以及初始的相机姿态INS和位置GNSS，计算改正数和改正后的值。

    参数：
        delta (numpy.ndarray): 一个数组，表示误差方程组中的未知数delta。
        INS (list): 一个列表，表示初始的相机姿态[omega, phi, kappa]。
        GNSS (list): 一个列表，表示初始的相机位置[XL, YL, ZL]。

    返回：
        INS_cor (list): 一个列表，表示改正后的相机姿态[omega_cor, phi_cor, kappa_cor]。
        GNSS_cor (list): 一个列表，表示改正后的相机位置[XL_cor, YL_cor, ZL_cor]。
    """

    omega_cor = delta[0][0]
    phi_cor = delta[1][0]
    kappa_cor = delta[2][0]
    XL_cor = delta[3][0]
    YL_cor = delta[4][0]
    ZL_cor = delta[5][0]

    omega = INS[0]
    phi = INS[1]
    kappa = INS[2]

    XL = GNSS[0]
    YL = GNSS[1]
    ZL = GNSS[2]

    GNSS_cor = [XL + XL_cor, YL + YL_cor, ZL + ZL_cor]
    INS_cor = [omega + omega_cor * (180 / math.pi), phi + phi_cor * (180 / math.pi), kappa + kappa_cor * (180 / math.pi)]
    return INS_cor, GNSS_cor


def space_resection(GNSS, INS, P, G, f, iteration, threshold):
    """迭代求解空间后方交会法。

       给定初始的相机位置GNSS和姿态INS，以及物相的坐标P，地面控制点坐标G，焦距f，迭代次数iteration和收敛阈值threshold，使用空间后方交会法迭代求解相机位置和姿态。

       参数：
           GNSS (list): 一个列表，表示初始的相机位置[XL, YL, ZL]。
           INS (list): 一个列表，表示初始的相机姿态[omega, phi, kappa]。
           P (numpy.ndarray): 一个数组，表示物相的坐标矩阵。
           G (numpy.ndarray): 一个数组，表示地面控制点的矩阵。
           f (float): 一个浮点数，表示相机的焦距。
           iteration (int): 一个整数，表示迭代次数。
           threshold (float): 一个浮点数，表示收敛阈值。

       返回：
           GNSS (list): 一个列表，表示迭代后的相机位置[XL, YL, ZL]。
           INS (list): 一个列表，表示迭代后的相机姿态[omega, phi, kappa]。
       """
    # 计算旋转矩阵M
    M = rotation_matrix(INS)

    i = 0
    # 进入迭代循环，直到满足收敛条件或达到最大迭代次数
    while True:
        # 构建系数矩阵B，常数向量L
        B, epsilon = error_equations(GNSS, INS, M, P, G, f)
        # 求解未知数向量X
        delta = solve_unknowns(B, epsilon)
        # 计算改正数和改正后的值
        INS_cor, GNSS_cor = correction(delta, INS, GNSS)

        # 更新原有的INS和GNSS值
        GNSS = GNSS_cor
        INS = INS_cor

        print("当前迭代次数：", i)
        print("相机位置：", GNSS)
        print("相机姿态：", INS)
        # 增加迭代次数
        i += 1

        # 判断是否满足收敛条件或达到最大迭代次数
        if np.linalg.norm(delta) < threshold or i > iteration:
            print("达到条件，结束迭代")
            # 返回相机位置和姿态
            return GNSS, INS



def mian(GNSS, INS, P, G, f, iteration, threshold):
    space_resection(GNSS, INS, P, G, f, iteration, threshold)


if __name__ == '__main__':

    GNSS = [1009.923, 1038.056, 649.6]
    INS = [0.0, 0.0, 102.83]
    P = {
        'A': {'x': 86.421,'y': -83.977},
        'B': {'x': -100.916,'y': 92.582},
        'C': {'x': -98.322,'y': -89.161},
        'D': {'x': 78.812, 'y': 98.123}
    }
    G = {
        'A': {'X': 1268.102,'Y': 1455.027,'Z': 22.606},
        'B': {'X': 732.181,'Y': 545.344,'Z': 22.299},
        'C': {'X': 1454.553,'Y': 731.666,'Z': 22.649},
        'D': {'X': 545.245, 'Y': 1268.232,'Z': 22.336}
    }

    f = 152.916
    iteration = 20
    threshold = 0.001

    mian(GNSS, INS, P, G, f, iteration, threshold)
