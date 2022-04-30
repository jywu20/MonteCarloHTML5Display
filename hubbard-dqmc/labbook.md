# 2022.4.28

从旧代码里面找来了一个Hubbard DQMC的程序`hubbard-old-code.jl`，看看能不能修改得generic一些。

先测试一下能不能用：取$U=1$，输出为
```
T -1.4777923120276673 0.004471842143327594
Mott 0.20837948635182507 0.000888249213972672
```
第一列为平均值，第二位为标准差。

在`hubbard-old-code-renamed.jl`中修改一些名称，使得后续的模块化更容易进行。
修改了一部分，取$U=2$，有

```
Progress: 100%|████████████████████████████████████████████████████████████████████| Time: 0:14:43
T -1.431692758497201 0.0051304951942906215
Mott 0.1799663182632188 0.0011940302855397808
```

初步计划：
- 将主循环移动到`main.jl`
- 构建格点的移动到`lattice.jl`
- 理想情况下，可以先写一个generic的版本用来计算$\ee^{- \Delta \tau T}$，然后写一个专属于正方晶格的版本做checkboard decomposition。

# 2022.4.30

目前还没有做的事情：
- 将自旋上和自旋下的各种函数合并（是否有这个必要？性能overhead？）

提交
```
U = 4.0
β = 4.0
```
的任务，slurm id为1262582。输出为`slurm-1262582.out`，应该是没啥问题。

将checkboard decomposition的无用代码（它们应该在计算$\exp(-\Delta \tau \bold{T})$乘以别的矩阵时出现，现在我先把$\exp(-\Delta \tau \bold{T})$算好存起来，那用不用checkboard decomposition没有差别）注释掉，重新做计算，扩大规模。

目前的规模：
```julia
L = 10 
lattice = SquareLattice2D(L)
n_sites = lattice.n_sites

t = 1.0
U = 4.0
β = 4.0
n_imtimes = 100
n_wrap = 10

n_heat = 1000
n_sweep = 4000

n_bin = 10
```

更新服务器上的代码，跑一下试一试。slurm id为1263693。

似乎`println`更新不及时；加入了一些`flush(stdout)`语句，再看看。
此时
```
t = 1.0
U = 2.0
β = 4.0
```
slurm id为1265098。

1265098超时了。看来最好还是从$L=4$开始做。

可以做的事情：
- 计算$L=10$时扫一次要多久

运行
```julia
@time sweep!(model, 1)
```
得到
```
 33.939589 seconds (3.29 M allocations: 12.192 GiB, 2.80% gc time)
```
确实很长……

按照
```julia
L = 10 
lattice = SquareLattice2D(L)
n_sites = lattice.n_sites

t = 1.0
U = 4.0
β = 4.0
n_imtimes = 100
n_wrap = 10

n_heat = 1000
n_sweep = 4000

n_bin = 10
```
的规模，总计扫了$1000 + 10 \times 4000 = 41000$次，需要的时间是
```julia
33.939589 * 41000 / 3600
```
也就是386个小时！预热的时间也要9小时。

我觉得我还是试一试不同的$L$下扫一次需要的时间为好。

| $L$ | 时间(s) |
| :------ | :------ |
| 4 | 0.735203 |
| 5 | 1.568580 |
| 6 | 2.041254 |
| 7 | 7.743210 |
| 8 | 12.953749 |
| 9 | 20.645465 |
| 10 | 33.939589 |

这意味着$L=4$的情况也需要8小时。

在加入磁化率计算之后，使用如下小规模计算：
```julia
#region Model parameters

t = 1.0
U = 2.0
β = 4.0

#endregion

#region Problem size

L = 4 
lattice = SquareLattice2D(L)
n_sites = lattice.n_sites
site_list = lattice.site_list

n_imtimes = 100
n_wrap = 10

#endregion

#region Sweeping control

n_heat = 100
n_sweep = 400
n_bin = 10

#endregion
```

似乎磁化那里有问题。设置$U=0$，只扫描一个bin，得到
```
Square lattice 2D Hubbard simple DQMC


t   =    1.0
U   =    0.0
β   =    4.0

n_sites   =   16
n_imtimes =   100
n_wrap    =   10

n_heat    =   100
n_sweep   =   400
n_bin     =   1

Heating up finished.


====================================================
E_kin            double_occ    magnetization
====================================================
-1.4993291872    0.2500000000    0.0507602959 

====================================================

E_kin      =   -1.4993291872038983 ± NaN
double_occ =   0.25000000000000255 ± NaN
mag        =   0.050760295887215943 ± NaN
```
然而`mag`那里应该是0.01才对。