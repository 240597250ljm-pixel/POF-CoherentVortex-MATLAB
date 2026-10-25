# POF平面涡旋与水平尺度间传递 MATLAB程序包

**v1.0-rc4：与修订稿Methods/SI一致的公开方法参考实现；尚未认证为论文全部图件的数值复现版。**

本次与10月24日版正文及SI对齐了平面Liutex、P80涡旋识别、有限窗口Voronoi、水平二分量Gaussian粗粒化、近场/背景区域统计、强Liutex R90与三阶纵向结构函数S3、帧级移动块Bootstrap等核心方法。`+pof/`为可复用方法，`pof.analyzeFrame`可将经QC的单帧物理速度场串接为完整单帧参考计算，MATLAB中执行`run_pof('methods')`查看。

依赖MATLAB R2022b+和Image Processing Toolbox。

```matlab
run_pof('status')
run_pof('test')
run_pof('example')
```

上述合成测试不等于重新计算2935帧论文原始实验，也不保证逐点复现主文和补充材料全部图件。实验数据不在包内。

**速度预平滑已与终稿说明对齐：** 正文II E和SI 1.1明确说明：各帧先扣除速度空间均值，再对U、V分别进行σ=0.8矢量网格的Gaussian预平滑（离散对称核±4σ），随后以FWHM ℓ执行单独的Gaussian粗粒化（±3σ核、3σ有效域排除）。`horizontalFlux`和`analyzeFrame`现在默认执行0.8网格预平滑，指定`preSmoothSigmaCells=0`可做不预平滑的敏感性比较。Liutex R场本身的0.8网格平滑是另一独立步骤。SI S6/S7的旧81775目录与主文81773目录之间的溯源问题仍独立保留。没有声称逐点复现冻结图件。

0604原始VEC U/V已是m/s；0605为pixel位移，标定一致已由实验人员确认。程序单位接口不会重复换算。

`reference/`为原始参考程序，不属于安全统一入口；历史47份源码单独压缩归档，避免误运行旧数据口径。仓库已经公开，但尚未选择代码再使用许可证；见`LICENSE_PENDING.md`。

**本轮公式核查改动：** 对齐Fig04/05生产源码，涡旋中心取连通域几何质心；Liutex R场平滑采用原始生产脚本的4σ离散截断；通量Gaussian核及有效域保留3σ。随机Voronoi零模型默认不在每个实验帧重置同一个随机种子。参见 `docs/METHOD_AUDIT.md`。
