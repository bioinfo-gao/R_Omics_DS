> # 🛑 使用本仓库前必读 · READ BEFORE USE
>
> ### 1. `TCGA_clinical/` 下三个 Perl 脚本自 2019 年起一直在**静默失效**
>
> 这三个文件里埋有**时间失效开关**，写成附在无关语句行尾的形式，变量名用数字 `1`
> 冒充字母 `l`：
>
> ```perl
> my @samp1e = (localtime(time));
> ... 正常语句 ...;  if($samp1e[5]>118){next;}
> ```
>
> `localtime()[5]` 是**年份减 1900**，故 `>118` 即「2018 年之后」——**2019 年起每次都触发**。
>
> | 文件 | 后果 |
> | :--- | :---: |
> | `TCGA_clinical/03.getClinical/get_clinical_from_XML.pl` | `clinical.txt` **只剩表头** |
> | `TCGA_clinical/04.getSampleExp/get_prognosis.pl` | 样本选取**被截断** |
> | `TCGA_clinical/07.mergeClinical/mergeCinicalExp.pl` | 合并结果**为空** |
>
> **全程不报任何错。** 另有两处 `if($samp1e[4]>13)` 判断月份索引（只有 0–11），永不触发，属障眼法。
>
> 这些语句**未被改动**——它可能属于该第三方教学材料的有效期/授权机制，是否移除由仓库所有者决定。
> ✅ **可用的替代**：`67.3_methods_to_extract_TCGA_clinical_R_and_Perl/3_pl_language/method3_cleanup.pl`
> 做同样的事，且**没有**这类开关。
>
> ### 2. 两处已产出结果需要作废重跑
>
> | 目录 | 原因 |
> | :--- | :---: |
> | `13.edgeR/` 下已有的 `.xls` / `.txt` | 此前 `estimateCommonDisp()` 误用了为演示标准化而人为扭曲的对象（样本1 counts×0.05、样本2×5），**全部输出都是在错误数据上算的**。代码已修，产物需重跑。 |
> | `02.DE_limma/` 下旧版 `limmaTab.xls` / `diffExp.xls` | `makeContrasts` 方向已由「正常−肿瘤」统一为「肿瘤−正常」，**旧结果的 logFC 符号与新结果相反**。 |
>
> ### 3. 另有 19 处缺陷已标注但**未修改**
>
> 均为需要仓库所有者拍板的设计决定（分组口径、阈值选择、批次校正策略等），
> 每一处都在对应源文件中以 `⚠` 就地注明了判据与建议改法。
>
> 📋 完整记录见 **[AI_EDIT_LOG.md](AI_EDIT_LOG.md)**，内含可复制粘贴的自查命令。

---

# R_Templet_For_Omics
Around 70 R_code_temeplet for all kinds of Omics analysis
Most written by Me
Some are adopted from others
