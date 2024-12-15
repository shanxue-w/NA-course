## 样条项目简介

建议多核编译，要不然有点慢。
使用`make -j4`可以一键编译运行项目，并生成项目报告。`make`是由以下几个组成
```bash
make all && make run && make plot && make doxygen
```

部分文件说明：
- `src/`
    - `include/`: 项目包含所有头文件
      - `PPoly.hpp`: 分段多项式类
      - `PPoly.tpp`: 分段多项式类模板实现
      - `PPInterpolate.hpp`: PP-form样条插值类
      - `PPInterpolate.tpp`: PP-form样条插值类模板实现
      - `BSpline.hpp`: B样条类
      - `BSpline.tpp`: B样条类模板实现
      - `BSInterpolate.hpp`: B样条插值类
      - `BSInterpolate.tpp`: B样条插值类模板实现
      - `Curve.hpp`: 曲线拟合类
      - `Curve.tpp`: 曲线拟合类模板实现
      - `BallFunction`: 球面函数类，主要包含球极投影和球面坐标两种。
      - `MinLeastSquare`: 线性最小二乘法，分析收敛阶会用到。
    - `src/`: 部分类的具体实现，对模板类的实现进行了具体化。
      - `PPInterpolate.cc`: 实例化了1,2,3阶样条插值类
      - `BSInterpolate.cc`: 实例化了1,2,3阶B样条插值类
      - `Curve.cc`: 实例化了2,3阶曲线拟合类
    - `A.cc`: 问题A的代码实现
    - `C.cc`: 问题C的代码实现
    - `D.cc`: 问题D的代码实现
    - `E.cc`: 问题E的代码实现
    - `F.cc`: 问题F的代码实现
    - `test.cc`: 测试高阶样条的代码
    - `HighOrder.cc`: 测试高阶样条的代码
    - `OrderAnalysis.cc`: 分析收敛阶的代码
    - `SplineTest.cc`: 更多样条的测试样例
    - `testjson.cc`: 测试`jsoncpp`的代码
    - `control.json`: `testjson.cc`的输入文件
    - `plot.py`: 画图的python脚本
- `figure/`: 存放所有的图片
- `doc/`
  - `design_html`: `doxygen`生成的html文档
    - `index.html`: 主页，包含所有数学公式
  - `design_pdf`: `doxygen`生成的pdf文档
    - `design.pdf`: 文档，包含所有数学公式
  - `report.tex`: 项目报告tex文件
  - `report.pdf`: 项目报告pdf文件
  - `Makefile`: 生成报告的Makefile
- `bin/`: 存放所有二进制文件
- `Makefile`: 项目的Makefile
- `README.md`: 项目说明文档
- `Doxyfile`: `doxygen`配置文件
- `mainpage.dox`: `doxygen`主页


## 使用

关于具体的使用。总体的框架是这样的，在`doxygen`生成的文档中也可以看出，这里具体解释一下。
```cc
std::vector<double> x = {0, 1, 2, 3, 4, 5};
std::vector<double> y = {0, 1, 2, 3, 4, 5};
std::vector<double> boundary_condition = {0,0};

PPInterpolate<3, double> pp(x,y,0,boundary_condtion);
BInterpolate<3, double>  bs(x,y,0,boundary_condtion);
```

其中`PPInterpolate`和`BInterpolate`里的第一个参数为阶数，即采用几阶样条来进行插值，默认使用的数据类型是`double`，所以此处`double`可以不写，当然也可以自己改成`long double`,`mpf_class`等数据类型。

构造函数的第三个参数为`method`，默认为0，有如下说明：
- $N=1$时，就算输入也不考虑。
- $N=2$时，0代表周期样条，1代表起点导数给定。
- $N=3$时，0代表周期样条，1代表完全样条，2代表自然样条
- $N \geq 4$时，0代表周期样条，1代表所有起点处的导数值给定。

这些需要记住，后续`Curve`拟合中的`method`参数和上述的参数描述完全一致。

`boundary_condition`不是很严格，长度也不需要严格匹配，只需要前$N-1$个和要求吻合就行。