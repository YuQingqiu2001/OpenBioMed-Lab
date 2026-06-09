# LLM 深度推理 vs NLP 提取 — 具体范例

> 参考: `paper-research` SKILL.md "绝对禁令" 章节

## 同一篇论文的两种分析对比

以 **Nature Medicine** 发表的 AAV 基因治疗 HoFH 论文 (PMID: 42243546) 为例。

### ❌ NLP 提取（禁止）

```json
{
  "claims": [{
    "claim": "Homozygous familial hypercholesterolemia (HoFH) is a rare autosomal disease characterized by severely elevated low-density lipoprotein cholesterol...",
    "evidence": "Homozygous familial hypercholesterolemia (HoFH) is a rare autosomal...",
    "synthesis": "基于摘要描述",
    "strength": 2
  }],
  "mechanism_cascade": "More than 80% of patients with HoFH carry LDLR mutations.",
  "hidden_axis": "待全文阅读后发现隐藏组织轴",
  "concept_innovation": "待全文评估概念创新性"
}
```

**问题**: 只复制了摘要第一句作为 claim；evidence 和 claim 完全一样；mechanism 是孤立的一句背景陈述；隐藏轴和概念创新全是"待…"——根本没有推理。

### ✅ LLM 深度推理（必须的标准）

```json
{
  "claims": [{
    "claim": "NGGT006（AAV8-密码子优化LDLR）在动物模型中可剂量依赖性地降低LDL-C",
    "evidence": "Ldlr-/-小鼠中验证；肝脏LDLR蛋白表达恢复；明确的剂量-反应关系",
    "synthesis": "动物模型是完全Ldlr敲除而非人源化点突变——可能高估了疗效，因为完全缺失LDLR的肝脏对转基因产物无免疫原性。人患者的部分功能突变可能引发对野生型LDLR的免疫应答",
    "strength": 3,
    "uncertain": "AAV递送的LDLR在人体中能否达到足够高的肝脏转导效率？预先存在的AAV8中和抗体会排除多少患者？"
  }],
  "mechanism_cascade": "AAV8-LDLR静脉注射→肝细胞表面AAVR识别衣壳→clathrin介导内吞→核内体逃逸→核内episomal dsDNA→密码子优化LDLR mRNA转录→rER翻译→高尔基体折叠/糖基化→细胞膜LDLR表达→LDL结合→网格蛋白内吞→LDLR再循环→溶酶体降解LDL→游离胆固醇→SREBP2-SCAP的Insig保留→SREBP2不能进入高尔基被切割→HMGCR转录↓→内源性胆固醇合成↓。**双重降脂机制**：增加清除(LDLR-LDL内吞) + 抑制合成(SREBP2反馈)",
  "hidden_axis": "基因治疗的降脂效果来自两个协同途径而非单纯补足LDLR。这解释了为何基因剂量低于预期也可达到临床意义的LDL-C降低。'通路放大效应'在代谢性疾病基因治疗中可能是一个普遍原则——修复一个节点通过内源性反馈放大效果",
  "concept_innovation": "从血友病(FVIII/FIX)扩展到常染色体隐性代谢病的首次AAV应用；Phase 1即纳入替代终点(LDL-C%)支持加速审批"
}
```

**差异总结**:

| 维度 | NLP | LLM 推理 |
|---|---|---|
| Claims | 复制摘要 | 自己的话 + 评估 + 不确定性 |
| 机制 | 一句背景 | 完整级联：18步从病毒到胆固醇反馈 |
| 隐藏轴 | "待全文" | 独立发现"通路放大效应"原理 |
| 跨文献 | 无 | 类比血友病 AAV 经验 |
| 用户收到的是 | "索引"（找到了但没有学习） | "学习"（真的理解和推理了） |
