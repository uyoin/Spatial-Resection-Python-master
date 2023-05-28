# Space-Resection-Algorithms
### [Video-Bilibili](https://www.bilibili.com/video/BV1L84y1M7AG/#reply167187697408) | [Online-PPT](https://kdocs.cn/l/cjniSJD1c8EG)  
<br/>

本仓库包含了用于实现空间后方交会法（Space Resection）的Python代码，这是一种重要的摄影测量学计算方法。代码实现了一种用于处理GNSS/INS系统的数据，并计算相机的位置和姿态的方法。

## 项目描述

本项目实现了以下几个主要功能：

- 旋转矩阵计算：根据GNSS/INS系统给出的姿态角计算旋转矩阵，描述相机坐标系和地面坐标系之间的关系。
- 误差方程构建：根据相机和地面控制点的坐标构建误差方程，形成误差矩阵和系数矩阵。
- 未知量求解：通过求解误差方程得到未知量delta。
- 改正值计算：根据未知量delta计算出的改正值，对相机的位置和姿态进行改正。
- 空间后方交会法：基于以上的步骤，通过迭代的方法实现空间后方交会法。

代码提供了完整的计算流程，从旋转矩阵的计算，到误差方程的构建，再到未知量的求解，最后计算改正值，实现空间后方交会法的全过程。

## 使用方法

首先，需要在python环境下安装numpy库，可以使用以下命令安装：

```bash
pip install numpy
```
运行代码需要在Python 3.6以上版本的环境下运行。

运行代码示例：
```python
python space_resection.py
```

如果你有任何问题，建议或者反馈，欢迎通过 "Issues" 功能与作者联系。
如果你觉得这个项目对你有帮助，欢迎 "Star" 支持作者！