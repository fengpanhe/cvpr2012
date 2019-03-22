# cvpr2012
论文 A learning-based framework for depth ordering 的复现

## 依赖准备

### 1. vlfeat

下载 [VLFeat](http://www.vlfeat.org/download/vlfeat-0.9.21-bin.tar.gz) 并解压到目录 `lib/vlfeat`


## Synthetic Dataset

1. 生成数据集 10 张

```matlab
generateSyntheticDataset(10);
```

2. 在 Synthetic Dataset 训练

```matlab
[precision1, precision2] = sdSSVMTrainTest(1);
```

3. TJ Test

```matlab
[ssvm_precision, mrf_precision] = TJTest('resources/SSVMmodel/ssvm.model', 'resources/SyntheticDataset', 'sd2', 1);
```
