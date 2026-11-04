# RCG Paper Figure/Table Formatting — source snapshot

Source: User's Notion page “RCG Paper Figure/Table Formatting”.
Retrieved: 2026-11-04. The Notion page remains authoritative if this copy becomes outdated.

# RCG manuscript — 作图规范 v1.1（整合工作版）

原工作版更新日期：2026-11-04。整合日期：2026-11-04。本文以 8 页《RCG manuscript — 作图规范 v1.0（工作版）》为主体，合并原 Notion 页面「RCG Paper Figure/Table Formatting」的 Color Code 内容。投稿时仍须按目标期刊要求调整尺寸与字号。主图与补图遵守同一规范。

<b>快速导航</b>

<table-of-contents />

<toggle heading="h2">
分组名称与配色

### 总原则与旧规则冲突的处理

颜色方案由图的比较目的决定，并在同一张图及其图例中保持一致；同时使用文字、固定位置或点形，避免仅靠颜色识别组别。

1. <b>同一基础模型内部比较不同处理条件</b>：Baseline 为蓝色（PyMOL <code>marine</code>），steering 为橙色（<code>orange</code>），RCG 为蓝绿色（<code>aquamarine</code>），RCG + steering 为粉紫色（<code>violet</code>）。
2. <b>同图比较不同基础模型及其处理条件</b>：采用原 Notion 页面「different model」的模型色系及下表的精确 HEX：AF3 棕色系、Boltz-1 蓝色系、Boltz-2 绿色系；非 steering 较浅，steering 较深。跨模型比较优先分面；必须叠加时仍使用下文规定的模型点形。
3. <b>旧 PDF 配色的明确替换</b>：v1.0 的 Okabe–Ito 四色表（Baseline <code>#0072B2</code>、steering <code>#009E73</code>、RCG <code>#E69F00</code>、RCG + steering <code>#CC79A7</code>）及“同一处理条件在不同模型中保持相同颜色”“两组比较保留上述原色”的要求，与 Notion 配色冲突，自本整合版起不再作为当前默认规则。两组比较应保留其所属情境的颜色，不因删减其他组而重新配色。处理条件名称、顺序及下文所有非配色规则继续有效。
4. <b>参考结构与界面色</b>：跨模型配色表中的 Ground truth <code>#CCCCCC</code> 用于一般参考数据；结构图中不同对象仍按「结构图规范」使用各自的灰色。白底、文字、轴线、连接线等界面颜色仍按「默认视觉样式」执行。

### 跨模型图：模型色系与 HEX

| 基础模型 | 图内完整组名 | HEX |
| --- | --- | --- |
| Ground truth | Ground truth | <code>#CCCCCC</code> |
| AF3 | AF3 | <code>#DACAB0</code> |
| AF3 | AF3 + steering | <code>#BCA987</code> |
| AF3 | AF3 + RCG | <code>#E6A45E</code> |
| AF3 | AF3 + RCG + steering | <code>#A96628</code> |
| Boltz-1 | Boltz-1 | <code>#C7E3EE</code> |
| Boltz-1 | Boltz-1 + steering | <code>#8AB2BE</code> |
| Boltz-1 | Boltz-1 + RCG | <code>#56B4E9</code> |
| Boltz-1 | Boltz-1 + RCG + steering | <code>#176A9C</code> |
| Boltz-2 | Boltz-2 | <code>#C6E5C6</code> |
| Boltz-2 | Boltz-2 + steering | <code>#8FB88F</code> |
| Boltz-2 | Boltz-2 + RCG | <code>#66C98B</code> |
| Boltz-2 | Boltz-2 + RCG + steering | <code>#237B4B</code> |

Boltz-2 的两种 RCG 颜色仅作为预留色码，不能据此在当前结果图中画出未开展的结果。Notion 的 <code>Boltz 1x</code>／<code>Boltz 2x</code> 名称在本表按 steering 明暗规则映射；正式图内组名统一采用本规范的 <code>+ steering</code> 写法。

### 同模型图：处理条件颜色

<table header-row="true">
<colgroup>
<col />
<col />
<col />
<col width="121.6015625" />
<col />
</colgroup>
<tr>
<td>
处理条件
</td>
<td>
面板已标明基础模型时的图内简写
</td>
<td>
颜色
</td>
<td>
PyMOL 内置颜色
</td>
<td></td>
</tr>
<tr>
<td>
Baseline
</td>
<td>
<code>Baseline</code>
</td>
<td>
蓝色
</td>
<td>
<code>marine</code>
</td>
<td>
<code>#007282</code>
</td>
</tr>
<tr>
<td>
steering
</td>
<td>
<code>(+ steering)</code>
</td>
<td>
橙色
</td>
<td>
<code>orange</code>
</td>
<td>
<code>#E69F00</code>
</td>
</tr>
<tr>
<td>
RCG
</td>
<td>
<code>(+ RCG)</code>
</td>
<td>
蓝绿色
</td>
<td>
<code>aquamarine</code> 
</td>
<td>
<code>#009E73</code>
</td>
</tr>
<tr>
<td>
RCG + steering
</td>
<td>
<code>(+ RCG + steering)</code>
</td>
<td>
粉紫色
</td>
<td>
<code>violet</code>
</td>
<td>
<code>#CC79A6</code>
</td>
</tr>
</table>

配色预览：

[The palette preview image is omitted from this public text snapshot; exact HEX values are copied above.]

### PyMOL 中使用跨模型精确颜色

可将以下指令放入 PyMOL 启动文件 <code>~/.pymolrc</code>，使自定义颜色随启动加载；编辑命令为 <code>nano ~/.pymolrc</code>。同模型方法色可直接使用上表的 PyMOL 内置名称。

```javascript
set_color af3, [218, 202, 176]
set_color af3_steering, [188, 169, 135]
set_color af3_rcg, [230, 164, 94]
set_color af3_rcg_steering, [169, 102, 40]

set_color boltz1, [199, 227, 238]
set_color boltz1x, [138, 178, 190]
set_color boltz1_rcg, [86, 180, 233]
set_color boltz1x_rcg, [23, 106, 156]

set_color boltz2, [198, 229, 198]
set_color boltz2x, [143, 184, 143]
set_color boltz2_rcg, [102, 201, 139]
set_color boltz2x_rcg, [35, 123, 75]
```

### 各基础模型的完整名称与纳入范围

| 基础模型 | Baseline | steering | RCG | RCG 与 steering 联用 |
| --- | --- | --- | --- | --- |
| AF3 | AF3 | AF3 + steering | AF3 + RCG | AF3 + RCG + steering |
| Boltz-1 | Boltz-1 | Boltz-1 + steering | Boltz-1 + RCG | Boltz-1 + RCG + steering |
| Boltz-2 | Boltz-2 | Boltz-2 + steering | 当前不纳入 | 当前不纳入 |

### 标签使用规则

- 正文首次出现写 <code>AlphaFold 3 (AF3)</code>；图内统一使用 <code>AF3</code>。
- 统一使用 <code>Boltz-1</code>、<code>Boltz-2</code>，以及小写 <code>steering</code>。
- 联用组固定写成 <code>+ RCG + steering</code>。
- 已在面板标题注明基础模型时，组名可以简写为 <code>Baseline</code>、<code>+ steering</code>、<code>+ RCG</code>、<code>+ RCG + steering</code>；跨模型总图使用完整名称。
- 四组总览顺序固定为 Baseline → + steering → + RCG → + RCG + steering。
- 两组比较沿用本节对应情境的颜色。
- 未开展或不适用的条件标为 <code>Not evaluated</code> 或省略，不画成零值。

### 跨模型编码

跨模型比较优先使用分面，顺序固定为 AF3、Boltz-1、Boltz-2。必须叠加在同一坐标系时，点形固定如下：

| 基础模型 | 点形 |
| --- | --- |
| AF3 | 圆形 <code>o</code> |
| Boltz-1 | 方形 <code>s</code> |
| Boltz-2 | 三角形 <code>^</code> |

点形保留给基础模型，不在其他图中改作处理条件或 seed 的编码。
</toggle>

<toggle heading="h2">
字体、字号与图幅

所有字体统一使用 Arial；下表字号均指最终排版尺寸下的字号。单栏图默认宽 90 mm，双栏图默认宽 180 mm。多面板图的宽度指整张图的宽度。

| 文字位置 | 典型内容 | 字号 | 字重及位置 |
| --- | --- | --- | --- |
| 整图顶部总标题 | <code>Effects of RCG on AF3 predictions</code> | 11 pt | 加粗，整图左对齐；用于独立展示版 |
| 子图上方标题 | <code>DockQ、Binder conformation</code> | 9 pt | 加粗，面板左对齐 |
| 子图编号 | <code>a、b、c</code> | 10 pt | 加粗，小写，面板左上方 |
| X 轴标题 | <code>D-Cα chirality correctness (%)</code> | 9 pt | 常规，轴下方居中 |
| Y 轴标题 | <code>Target-aligned binder RMSD (Å)</code> | 9 pt | 常规，轴左侧 |
| X／Y 轴刻度数字 | <code>0、0.5、1.0</code> | 8 pt | 常规 |
| 分类轴上的组名 | <code>Baseline、(+ RCG)</code> | 8 pt | 常规，优先水平排列 |
| 图例文字 | 各组名称 | 8 pt | 常规，黑色文字配色块或点形 |
| 图例标题 | <code>Model</code>、<code>Condition</code> | 8 pt | 加粗；含义明确时省略 |
| 连续变量色条标题／刻度 | 名称及数值 | 9／8 pt | 常规，注明单位 |
| 图内统计标注 | <code>n = 39</code>、均值、Δ、置信区间 | 8 pt | 常规 |
| 结构图标签 | 链名、残基名、体系编号 | 8 pt | 常规 |
| 结构图下方指标文字 | <code>DockQ = 0.650</code> | 8 pt | 常规，各面板位置一致 |
| 比例尺文字 | <code>10 Å</code> | 8 pt | 常规 |
| 整图底部短说明 | <code>Dots represent individual systems.</code> | 7.5 pt | 常规，左对齐 |
| 正文中的完整 figure caption | 图题及完整说明 | 9 pt | 放在正文中，首句可加粗 |

补充要求：

- 论文图件默认保留子图标题，整图标题放在正文 figure caption 中。单独发给合作者或用于讨论时，可以增加 11 pt 的整图顶部标题。
- 图内底部短说明控制在 1–2 行；完整方法、统计和筛选规则写入正文图注。
- 分类组名过长时换行，不通过不断缩小字体解决。
- “坐标轴标题”和“刻度数字”分开设置：前者 9 pt，后者 8 pt。
- 除底部短说明外，常规文字不低于 8 pt。
- 合图时保留各面板相同的实际字号，不分别任意缩放。例如，先按 360 mm 宽、9 pt 文字作图，再缩到 180 mm 宽，文字实际只有 4.5 pt。应直接按最终宽度制作；若必须整体缩放，同步计算字体和线宽。
- 以上是团队工作版标准。投稿时以目标期刊要求为准；例如 Nature 图件规范对最终图件另有较小字号要求，不能把本表直接视为其投稿规范。
</toggle>

<toggle heading="h2">
指标英文名称与定义

箭头表示优劣方向：↑ 越高越好，↓ 越低越好。箭头属于展示提示，不属于指标名称本身。

| 指标 | 统一英文名称 | 单位／方向 | 定义 |
| --- | --- | --- | --- |
| D-Cα 手性正确率 | D-Cα chirality correctness | %，↑ | 正确的目标 D-Cα 位点数，占预先指定目标 D-Cα 位点数的比例 |
| 全部 D-Cα 正确的体系比例 | Systems with all D-Cα sites correct | %，↑ | 所有目标 D-Cα 位点均正确的体系数，占纳入评估体系数的比例 |
| D-Cα 角度裕量通过率 | D-Cα margin pass rate | %，↑ | 达到指定手性角度裕量的目标 D-Cα 位点比例；当前阈值为定向角 q ≥ 15° |
| Target 自身构象偏差 | Target self-aligned backbone RMSD | Å，↓ | Target 单独拟合到参考 target 后，其对应主链原子的 RMSD |
| Binder 自身构象偏差 | Binder self-aligned backbone RMSD | Å，↓ | Binder 单独拟合到参考 binder 后，其对应主链原子的 RMSD |
| Target 对齐后的 binder 偏差 | Target-aligned binder backbone RMSD | Å，↓ | 仅使用 target 进行拟合，将同一变换应用于 binder 后计算其主链 RMSD |
| DAA 链自身构象偏差 | DAA-chain self-aligned backbone RMSD | Å，↓ | 含目标 DAA 的链分别独立拟合到对应参考链后，计算主链 RMSD |
| 对接质量 | DockQ | 0–1，↑ | 综合原生接触恢复、界面 RMSD 和 ligand RMSD 的对接质量评分 |
| 整体局部距离一致性 | Global heavy-atom lDDT | 0–1，↑ | 在指定共同重原子范围内，根据预测与参考原子间距离的一致性评分 |
| 界面局部距离一致性 | Interface heavy-atom lDDT | 0–1，↑ | 使用指定跨链重原子对计算的局部距离一致性评分 |
| 所选肽键约束违反率 | Selected peptide-bond constraint violation rate | %，↓ | 违反允许距离范围的所选肽键约束数，占实际评估约束数的比例 |
| 所选残基几何约束违反率 | Selected residue-geometry constraint violation rate | %，↓ | 违反允许范围的所选残基内部几何约束数，占实际评估约束数的比例 |

### 坐标轴命名

轴标题采用“英文名称 + 单位”，例如 <code>D-Cα chirality correctness (%)</code>、<code>Binder self-aligned backbone RMSD (Å)</code>、<code>Target-aligned binder backbone RMSD (Å)</code>、<code>DockQ</code>。名称较长时允许换行。只有在面板标题和图注明确完整含义时，坐标轴才可以简写为 <code>RMSD (Å)</code>。
</toggle>

<toggle heading="h2">
指标计算与汇总口径

指标计算口径按以下规则记录，不能仅统一名称而省略定义。

- <b>手性判定</b>：采用经过目标手性和原子顺序校正的定向角 q，当前分析以 q > 0° 判为方向正确。q = 0° 不计为正确。角度裕量 q ≥ 15° 单独报告，不能替代一般手性正确率。无相应手性中心的残基不纳入分母。
- <b>位点覆盖</b>：同时保存目标位点数、可评估位点数和正确位点数。缺失位点不计入正确数，也不能标成“已确认手性反转”；图注说明其对分母的处理。
- <b>手性汇总</b>：默认跨体系统计采用体系等权，即先计算每个体系的正确率，再计算体系间平均值。合并所有位点计算的结果另标为 <code>Pooled-site D-Cα chirality correctness (%)</code>。
- <b>RMSD</b>：主链原子统一为 N、Cα、C、O；使用对应共同原子和不允许镜像的刚体拟合。记录 target／binder 的链身份与分配规则。
- <b>Binder 两种 RMSD</b>：self-aligned 用于描述自身构象；target-aligned 同时受到自身构象和相对位置影响。两者不混用名称。
- <b>多界面汇总</b>：沿用现有四组分析时，体系级 DockQ 取纳入界面的均值；target／binder RMSD 取界面级结果的中位数；DAA-chain RMSD 取相关链独立拟合结果的均值。另有分析口径时须单独标注。
- <b>DockQ</b>：记录软件版本、评分模式、链映射和界面范围。全部界面与 DAA 相关界面的结果分开命名。
- <b>当前 lDDT 实现</b>：使用共同重原子、不同残基、参考距离约 15 Å 邻域，以及 0.5／1／2／4 Å 的距离误差阈值；按原子对平均。Interface 版本仅保留跨链原子对。图注注明 <code>distance-only, pair-weighted implementation</code>，不与其他实现直接混合。
- <b>几何约束</b>：记录具体约束列表与上下界。当前肽键距离约束按其参数中的均值 ± 3σ 判断；残基几何约束按各自上下界判断。分母为实际评估的约束数，无可评估约束时记为 NA。
- 当前 lDDT 实现不含化学对称性处理和立体化学处罚；上述几何违反率也只覆盖列出的约束，不能直接称为“整体物理合理率”。这些定义对应现有分析代码中的指标计算与汇总（原 PDF 所列内部脚本路径在公开副本中省略）。
</toggle>

<toggle heading="h2">
统计图展示规则

每张图首先说明一个点代表什么，再说明统计和选模方式。

| 项目 | 统一要求 |
| --- | --- |
| 观测单位 | 明确为体系、体系 × seed、位点或界面 |
| 配对比较 | 使用相同体系；连线只连接同一体系的对应结果 |
| 多 seed | 说明按 seed 分面、体系内汇总或其他处理方式；不能把重复 seed 当成独立体系 |
| 输出选择 | 写明最终导出模型、置信度选模、指定优化状态或其他规则 |
| Oracle 选模 | 使用参考结构指标选模时，明确标注，例如 <code>DockQ oracle selection</code> |
| 同一结构的指标 | 所有指标来自同一个选定结构文件 |
| 缺失结果 | 记录为 NA 并注明原因；配对分析使用对应可比较样本，不填成零 |
| 箱线图 | 中线为中位数，箱体为第 25–75 百分位，须线为 1.5×IQR 范围内的最远观测值 |
| 原始点 | 保留展示；隐藏箱线图的重复离群点符号时，原始点仍需完整呈现 |
| 误差条 | 明确是 SD、SEM 或 95% CI；统一一张图内的含义 |
| 置信区间 | 对体系级比较，通常在体系层面重采样，并保留组间配对关系 |
| 变化量 | 统一为 Δmetric = 新条件 − 对照条件 |

ΔDockQ > 0 表示提高，ΔRMSD \< 0 表示降低。正确率的差值使用百分点，例如从 80% 到 90% 写为 <code>+10 percentage points</code>。

### 数值显示精度

| 数值类型 | 默认显示精度 |
| --- | --- |
| 正确率、通过率、违反率 | 1 位小数 |
| RMSD | 2 位小数 |
| DockQ、lDDT | 3 位小数 |
| 样本量、位点数 | 整数 |

统计计算保留完整精度。效应小于默认显示精度时，可以增加小数位，但同一面板保持一致。
</toggle>

<toggle heading="h2">
坐标轴、线条、数据点与布局

### 默认视觉样式

| 元素 | 默认设置 |
| --- | --- |
| 背景 | 白色 <code>#FFFFFF</code> |
| 文字和主要坐标轴 | 深灰 <code>#222222</code> |
| 坐标轴线宽 | 0.7 pt |
| 数据线／轮廓线 | 1.0 pt |
| 配对连接线 | 0.6 pt，浅灰 <code>#B8B8B8</code> |
| 参考线 | 0.7 pt，灰色虚线 |
| 数据点大小 | 默认约 3.5 pt；同类图保持一致 |
| 数据点透明度 | 默认 0.85 |
| 分布填充透明度 | 默认 0.25 |
| 网格 | 默认关闭；必要时仅保留淡灰色辅助网格 |
| 图例边框 | 默认关闭 |
| 顶部和右侧边框 | 常规统计图默认隐藏 |

### 坐标范围与布局

- 同一指标、用于直接比较的面板尽量使用相同坐标范围。
- 正确率默认 0–100%；DockQ 和本文 lDDT 默认 0–1。为防止边界点被裁切，可以留少量显示边距。
- RMSD 原始值从 0 附近展示并覆盖全部数据；ΔRMSD 图保留正负范围和零参考线。
- 局部放大图明确标注放大范围，同时保留完整范围的信息。
- 同一指标的方法 A／B 对照散点图，x、y 范围一致，并加入 y = x 参考线。
- 数值轴不添加改变读数的抖动；分类轴可进行适量横向抖动。
- 分类组名优先水平排列，长名称换行。
- 同一行面板保持相同绘图区高度；同一列面板保持相同宽度。
- 长标签或多面板导致拥挤时，优先调整布局、留白或图高，保持字号不变。
- 大量数据集中在 0% 或 100% 时，保留原始点和数量信息，避免密度平滑造成误解。
</toggle>

<toggle heading="h2">
结构图规范

结构图沿用本规范中与比较情境对应的处理条件颜色，并固定参考结构与显示方式。

| 对象 | 默认颜色／表示 |
| --- | --- |
| Reference target | 浅灰 <code>#C7C7C7</code> |
| Reference binder | 深灰 <code>#555555</code> |
| Predicted target | 中灰 <code>#969696</code> |
| Predicted binder | 对应比较情境的处理条件颜色 |
| 错误手性位点 | 黑色叉号、标记圈或局部放大，并提供图例 |

- 同一案例的各处理组使用相同视角、缩放、裁切范围和结构表示方式。
- 展示复合物相对位置时，按 target 对齐，binder 随同一变换移动。
- 展示 binder 自身构象时，单独拟合 binder，并注明 <code>Binder self-aligned</code>。
- 同类结构统一采用相同 cartoon、tube 或 stick 表示。
- 结构图下方指标使用 8 pt，各面板保持相同顺序，例如 DockQ 在前、RMSD 在后。
- 专门展示 target 构象差异时，可以调整预测 target 的颜色，并在图例中明确说明。
- Ramachandran 图固定使用 φ／ψ 坐标和 −180° 至 180° 范围，注明 D 残基参考及位点筛选条件。
</toggle>

<toggle heading="h2">
完整图注模板

每张图的完整图注都应包含样本、选择和统计信息。替换方括号内容：

```txt
Figure X. [Main finding or comparison].
[Panel descriptions.] Each point represents [one system / one system–seed pair / one residue]. Analyses include [N] matched systems and [S] seeds per system. Structures were selected using [selection rule]. [Metric definitions and alignment procedure.] Boxes indicate the median and interquartile range; whiskers extend to the most extreme observations within 1.5×IQR. [Definition of error bars or confidence intervals.] [Missing-data handling and relevant evaluation restrictions.]
```

图内底部只保留简短说明，例如：<code>Dots: individual systems. Boxes: median and interquartile range.</code>
</toggle>

<toggle heading="h2">
交付文件

- 交付文件保留可编辑图件、预览图和对应数据。
- 统计图导出 PDF 或 SVG，保留可编辑文字并嵌入字体。
- 讨论预览导出 PNG，300 dpi；提高 DPI 不会改变 pt 字号，也不能修复源图字体过小的问题。
- 结构渲染图按最终使用尺寸导出，至少 300 dpi，建议 600 dpi。
- 统一 RGB 色彩空间、白色背景。
- 文件名包含图号、内容和版本，例如 <code>Fig05_AF3_RCG_DockQ_v01.pdf</code>。
- 同时保留作图脚本、源数据和选模清单。
- 输出后按最终宽度检查：文字是否清楚、是否裁切、点是否完整、图例是否一致。
- 配色和标签集中维护；修改规范后同步更新主图与补图。
</toggle>

