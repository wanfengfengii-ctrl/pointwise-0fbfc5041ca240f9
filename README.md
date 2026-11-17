# 晶格材料实验室原胞加工账本

Python 3.11，仅标准库。实验记录整理人员提供晶格、原子及已确认的周期键；
加工员登记完整送件和更正，下游模拟准备人员读取精确几何成果或离库核账。
本软件处理原胞，不猜键、不取最短像、不进行化学或结构优化。

## 一天的加工流程

1. 接收产品ID和完整送件，为原始容器建立独立快照。
2. 把晶格长度精确换算为nm，保留原单位审计材料。
3. 核验可逆晶格和半开坐标，建立原子身份索引。
4. 核验周期键端点、身份和整数偏移，提交新不可变修订及当前指针。
5. 下游按明确产品和修订固定任务，任务不追随当前指针。
6. 计算体积、原子位置和逐键账，核对身份覆盖与物理值。
7. 在同一事务提交任务、成果头和成果体，下游只读取完整成果。
8. 按固定修订和整数S派生超胞，每个源原子按商类产生半开副本、
   每条源键在每个起点副本上交接整数偏移，整体持久保存后按身份回审计。
9. 保存更正并分别计算新旧版本，关闭重开，导出旧版证据供离库复算。

## 运行

从项目根目录执行：

```sh
python3.11 -B verify.py
python3.11 -B demo.py
```

verify发现并执行全部tests/test_*.py，断言失败返回非零；随后执行
`supercell_cases.workflow()` 的联合精确核账（新鲜发布、关闭重开、
离库复算三组证据）；`python3.11 -B supercell_cases.py --verify`
是同时重跑原回归的等价脚本入口。
demo默认在测试临时目录建库，实际接纳样例、发布、登记更正、重开并导出旧成果；
它用Fraction从JSON证据独立复算，然后输出POINTWISE_SCENARIO单行JSON。
随后以同一临时库发布斜超胞（S=[[2,1,0],[0,1,0],[0,0,1]]），关闭重开，
再只用导出JSON离库复算，输出N9场景。
不需要安装、联网或准备数据。运行清单见.pointwise-runtime.json。

直接核验超胞契约本身：

```sh
python3.11 -B supercell_cases.py --fixtures
python3.11 -B supercell_cases.py --verify
python3.11 -B supercell_cases.py --demo
```

对已有数据库演示指定修订：

```sh
python3.11 -B demo.py --database ./lab.sqlite --product specimen --revision 1
```

此方式不装载样例或更正来源，不发布原胞成果，只在该固定修订上发布斜超胞、
重开读取并离库核账。默认临时demo不修改项目样例。数据库初始化不灌样例、
不重置用户数据；兼容V1库重复打开保持所有修订、当前指针和成果，不兼容模式
抛ValueError。打开旧V1库时把`tasks.calculation`的CHECK从仅primitive-v1
扩为同时接受supercell-v1：仅重建该表结构，任务ID、AUTOINCREMENT序列、
成果行与成果体逐字节不变，重建在单个事务内完成。

## 公开API和输入

```python
import json
from crystal import Archive, Tasks, Results

with Archive("lab.sqlite") as archive:
    with open("examples/day_cell.json", encoding="utf-8") as stream:
        cell = json.load(stream)
    revision = archive.put("specimen", cell)
    result_id = Tasks(archive).run_geometry("specimen", revision)
    result = Results(archive).read(result_id)
    evidence_json = Results(archive).export(result_id)
    assert result["bonds"][0]["length_squared"].numerator == 61

    S = [[2, 1, 0], [0, 1, 0], [0, 0, 1]]
    super_id = Tasks(archive).build_supercell("specimen", revision, S)
    supercell = Results(archive).read(super_id)
    assert supercell["source"]["calculation"] == "supercell-v1"
    supercell_json = Results(archive).export(super_id)
```

Archive(path=":memory:")支持上下文管理和幂等close()。
put(product_id,cell)提交完整替换版本，返回从1起递增的产品内修订号；
相同内容再次提交也登记新修订，摘要不用于自动去重。
read_revision(product_id,revision)返回规范nm完整cell字典；
current_revision(product_id)只查询当前指针。
Tasks(archive).run_geometry(product_id,revision)返回已提交成果整数ID；
每次正常运行可新增成果，必须显式指定修订。
Tasks(archive).build_supercell(product_id,revision,S)返回超胞成果整数ID；
S为3×3整数（非bool）列表/元组，元素绝对值≤3，1≤|det S|≤32，
源修订至多40原子、80键，空集合法。非法参数抛ValueError，
合法但缺失来源抛KeyError，多个参数同时非法时不保证检查顺序。
Results(archive).read(result_id)读取持久成果；
export(result_id)返回完整JSON文本，不是文件路径。

cell恰含unit、lattice、sites、bonds四个必填字段：

- unit为nm、angstrom或pm，长度比例为1、1/10、1/1000，规范单位为nm。
- lattice是按行书写的3×3数组，数学上列为基向量A，行列式非零；
  负行列式合法。只缩放晶格，不缩放分数坐标或周期偏移。
- sites是{id,fractional}列表。id为非空字符串，fractional为[0,1)内三元组。
- bonds是{id,u,v,offset}列表。端点引用本修订原子，offset严格为三个Python整数。
  键(id,u,v,t)连接u和整数平移t处的v。

有理数接受非布尔整数、Fraction或两个非布尔整数构成的列表/元组[n,d]，
d不得为0。浮点、bool、None、缺字段及多余字段拒绝。
偏移不接受Fraction或有理数对，即使其值为整数。
产品、原子、键ID按原字符串比较，不裁剪或大小写折叠；原子和键分别要求唯一。
原子列表与键列表可为空，但任何键端点必须存在。零位移键合法。
同坐标不同原子、重复端点不同键ID、自像键、反向键均独立保留。
数值输出按身份排序。规范cell及所有JSON中的有理数为约分正分母整数对；
周期偏移、版本和身份编号保持整数。

## 结果、口径和证据

read返回独立嵌套容器，数值是Fraction：
source={product_id,revision,digest,calculation}，calculation为primitive-v1；
units={length:nm,length_squared:nm2,volume:nm3}；
volume为|det A|；
positions逐项为{id,cartesian}，位置为A f；
bonds逐项为{id,u,v,offset,displacement,length_squared}，
位移为A(f_v+t-f_u)，长度平方为位移内积。
不使用浮点容差、开方或坐标去重。

手算样例A=[[1,1/2,0],[0,1,0],[0,0,2]]，
u=(0,0,0)，v=(1/2,1/2,0)，t=(0,1,0)：
体积2，v位置(3/4,1/2,0)，位移(5/4,3/2,0)，长度平方61/16。
这些值实际经过任务发布、读取和导出，tests/test_geometry.py逐项核对。

export包格式为crystal-primitive-evidence、version=1。
source包含product_id、revision、digest、完整规范cell与原始送件raw_cell；
result包含与read相同的完整字段，但Fraction编码为[n,d]。
原始送件保留单位和送件顺序，Fraction转换为整数对，输入元组转换为JSON列表。
摘要为规范cell的SHA-256：字段固定、原子和键按ID排序，
json.dumps(ensure_ascii=True,sort_keys=True,separators=(",",":"),allow_nan=False)，
ASCII字节；摘要不包含任务ID、成果ID、时间或当前指针。
源码tests/test_export.py包含只依赖json、hashlib、Fraction的独立复算参照，
不会调用产品的矩阵、几何或摘要函数。导出前验证原始与规范来源一致，
并验证原子、键的逐项覆盖和数值。

## 超胞表示与交接不变量

超胞成果read恰含source、units、matrix、cell、origins、volume、positions、bonds：

- source与V1四字段同形，calculation为supercell-v1；digest与cell、raw_cell
  始终指向固定源原胞摘要，不指向派生结构。
- matrix为整数行列表S本身。cell为Archive.read_revision同形的规范nm字典，
  其中晶格按行等于A S（A为源规范晶格），有理数为约分正分母整数对。
- origins.sites逐项{id,source_id,translation}，origins.bonds逐项
  {id,source_id}；id引用派生cell，source_id引用固定源原子或源键，
  translation为整数三元列表。
- volume=|det(A S)|；positions为AS乘派生半开坐标，bonds逐项
  {id,u,v,offset,displacement,length_squared}，位移为
  AS(y_v+T-y_u)，长度平方为位移内积。排序按身份，口径与V1相同。

派生不变量（不使用浮点或容差）：

1. 半开副本：每个源原子f恰有D=|det S|个派生副本y∈[0,1)^3，配整数来源n
   满足S y=f+n；同一源原子的D个n在Z^3/SZ^3中两两不同，恰好覆盖整个商群。
   不同源原子即使坐标相同也各自产生独立身份，不焊接。
2. 键偏移交接：每条源有向键在它每个起点副本上恰产生一条派生键，
   每条为(源键ID, 起点副本)的单一组合，不合并自像键、反向键、零位移键
   或重复端点键。派生终点取该源键终点原子的唯一副本，T严格为整数，
   且S T=offset+start−end，于是源位移A(f_v+t−f_u)被精确保留。
3. 身份独立：派生原子与键ID是各自非空唯一字符串，命名规则不属于契约；
   下游只能经origins回指源身份。

持久保存固定源修订、S、派生cell、origins和完整账本；read/export只解码
校验，不在读取时重新枚举副本或重建超胞。每次读取都做独立精确核账：
核对摘要、AS=A S、体积、S y=f+n、副本商类覆盖、源键起点覆盖、
端点来源身份、T的整数性与位移/长度平方。任一证据缺失、重复或错配都拒绝。

export格式为crystal-supercell-evidence、version=1；source与V1导出完全同形
（含固定源规范cell与raw_cell），result与read同形，仅把Fraction编码为
约分正分母整数对。matrix、origins和派生cell全部来自持久体。

手算样例A=单位阵，a的f=(0,3/4,0)，b=(a,a,(0,1,0))，
S=[[2,1,0],[0,1,0],[0,0,1]]：两副本的(n,y)为
((1,0,0),(1/8,3/4,0))与((2,0,0),(5/8,3/4,0))；从n=(1,0,0)终点n=(2,0,0)，
T=(−1,1,0)；从n=(2,0,0)终点n=(1,0,0)，T=(0,1,0)；两键位移均为(0,1,0)。
tests/test_supercell.py与demo的N9场景用Fraction从导出JSON独立复算这组值。

负行列式合法：商类键用±adj(S)模|det S|的规范代表判定，不依赖左右手符号；
派生位移仍走V1的A S与半开坐标，故反向晶格与V1几何一致。

## 错误与失败边界

非法产品、结构、数值、单位和参数抛ValueError。
合法但不存在的产品/正整数修订/成果抛KeyError。
计算或数据库写入异常原样传播；本次修订与指针、任务与成果各自原子发布，
失败不留下空产品或半成品。返回值、原始输入的后续修改不写回归档。
旧修订及成果不因更正、读取或重开改变。证据缺失或来源不一致拒绝导出。
超胞派生在发布事务之外完成；派生或写入任一步失败都不新增任务、成果、
来源映射或原胞修订，也不移动当前指针。S非法（含bool、Fraction、浮点、
非整数、越界、奇异或|det S|>32）或源超过40原子/80键抛ValueError；
合法S但产品/修订缺失抛KeyError。读取或导出发现持久体被篡改（含派生cell、
origins、账本任何字段与源不一致）抛ValueError。
API不提供删除旧来源、隐式latest、并发调度或崩溃恢复协议。
archive.connection是底层SQLite连接；外部SQL破坏数据不属于正常写入API，
但读取会拒绝缺失成果体和头体错配。测试用真实SQLite触发器在写入中途拒绝，
验证事务回滚，并不替换计算函数。

## 代码与回归

接纳：intake/input_model/exact_numbers/units/input_digest。
几何：lattice/rational_matrix/sites/cell_geometry。
周期图及核账：bonds/periodic_graph/bond_geometry/bond_audit。
修订：archive/revision_store/revision_codec/schema/database。
成果：tasks/executor/results/export/result_store/result_codec。
超胞：supercell_matrix/supercell/supercell_codec/supercell_audit。
所有模块均接入公开API，输入和结果值对象使用不可变元组与冻结记录。

N1、N3在test_intake；N2在test_geometry；N4、N5在test_bonds；
N6在test_tasks；N7在test_revisions；N8在test_export。
test_codec补充相同业务边界的精度、完整字段和模式版本验证。
test_supercell覆盖手算两副本、五种S（含负行列式与非对角）、
空源、参数边界、重复发布、事务中途失败、关闭重开、更正后旧成果冻结、
持久体篡改拒绝以及旧V1库升级，并用只依赖标准库的独立复算核对导出证据。
完整套件还覆盖固定seed=1729的有限三角晶格、左右手晶格、
真实发布/更正中途失败、大整数有理数、重复初始化及修改返回容器后重开。
测试创建的数据库仅位于测试临时目录，结束后清理。
