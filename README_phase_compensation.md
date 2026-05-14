# THz MIMO 雷达运动抖动相位补偿代码说明

本代码用于补偿 THz MIMO 雷达在二维运动扫描过程中，由硬件或机械抖动引起的 RX 通道相位偏差。默认输入数据为复数频域数据：

```matlab
data(rx, tx, freq, y, x)
```

其中 `rx` 是接收通道，`tx` 是发射通道，`freq` 是频率采样点，`y/x` 是二维扫描位置索引。

## 文件说明

- `phaseCompensateMIMO.m`：核心相位估计与补偿函数。
- `demo_phaseCompensateMIMO.m`：合成数据示例，演示相位抖动注入、补偿和可视化对比。

## 核心算法

当前函数默认使用 PGA-like 相位梯度自聚焦方法：

```matlab
opts.method = 'pga';
```

PGA 的核心思想不是直接估计绝对相位，而是先估计孔径方向的相位梯度，再通过积分恢复相位误差。代码首先对频率维做 IFFT，得到距离像数据；随后选择能量最强的若干距离单元作为强散射点支撑：

```matlab
rangeData = ifft(estData, [], 3);
```

然后沿二维扫描孔径的 x/y 方向计算相邻位置的复数相关：

```matlab
prodX = pgaData(:,:,:,:,2:end) .* conj(pgaData(:,:,:,:,1:end-1));
prodY = pgaData(:,:,:,2:end,:) .* conj(pgaData(:,:,:,1:end-1,:));
```

其相位分别表示 x/y 方向的相位梯度。代码再将梯度相对参考 RX 通道去公共项，并通过加权积分恢复 `phaseEst(rx,y,x)`。最后施加反相位：

```matlab
dataCorr(rx,:,: ,y,x) = data(rx,:,: ,y,x) .* exp(-1j * phaseEst(rx,y,x))
```

函数也保留了原来的参考通道相干法：

```matlab
opts.method = 'coherence';
```

该方法以一个 RX 通道作为参考通道，默认 `refRx = 1`。对每一个扫描位置 `(y, x)`，将其他 RX 通道与参考 RX 通道在 TX 和频率维度上做复数相关：

```matlab
c(rx,y,x) = sum(data(rx,:,: ,y,x) .* conj(data(refRx,:,: ,y,x)))
```

相关结果的相位：

```matlab
phaseEst(rx,y,x) = angle(c(rx,y,x))
```

即为该 RX 通道相对参考 RX 的相位偏差估计。随后代码会对相位图进行空间相位展开、加权平滑，并施加补偿：

```matlab
dataCorr(rx,:,: ,y,x) = data(rx,:,: ,y,x) .* exp(-1j * phaseEst(rx,y,x))
```

## 使用方式

最简单调用方式：

```matlab
opts = struct();
opts.refRx = 1;
opts.dimOrder = 'rx_tx_freq_y_x';
opts.freqIdx = [];

[dataCorr, phaseEst, info] = phaseCompensateMIMO(data, opts);
```

常用参数：

- `opts.refRx`：参考 RX 通道编号。
- `opts.method`：相位估计方法，默认 `'pga'`；也可设为 `'coherence'`。
- `opts.freqIdx`：用于估计相位的频点范围，例如 `12:84`。
- `opts.pgaInputDomain`：PGA 输入域，默认 `'frequency'`；如果第三维已经是
  SAR `range_fft` 输出的距离单元，应设置为 `'range'`，避免再次 IFFT。
- `opts.phaseModel`：PGA 相位模型，默认 `'rx_relative'`；SAR 联动实验可使用
  `'aperture_common'` 或 `'mimo_aperture'` 来补偿公共孔径相位误差。
- `opts.nominalGradientX/Y`：可选的名义孔径相位梯度，用于从相邻孔径相位中扣除
  已知几何相位，只估计残余误差相位。
- `opts.removeLinearPhase`：是否去除估计相位中的常数/线性相位项。SAR autofocus
  中线性相位主要对应图像平移，通常不作为散焦误差补偿。
- `opts.pgaRangeBins`：PGA 使用的距离单元编号；为空时自动选择强散射距离单元。
- `opts.pgaNumRangeBins`：自动选择强散射距离单元时使用的数量。
- `opts.smoothWindow`：二维相位图平滑窗口，例如 `[7 7]`。
- `opts.unwrapSpatial`：是否在扫描平面上展开相位。
- `opts.removeMeanPhase`：是否去除每个 RX 的平均相位，只保留位置相关抖动。
- `opts.minMagnitudePercentile`：低幅度位置阈值百分位，弱信号位置会被降权。

## 输出说明

- `dataCorr`：补偿后的复数频域 MIMO 数据，维度顺序与输入一致。
- `phaseEst`：估计到的相位补偿量，维度顺序与输入一致，但 TX 和频率维为单例维度。
- `info.phaseEstStandard`：标准顺序下的相位估计，尺寸为 `rx x y x`。
- `info.coherenceBefore`：补偿前各 RX 与参考 RX 的相干性。
- `info.coherenceAfter`：补偿后各 RX 与参考 RX 的相干性。

## Demo 图像含义

运行：

```matlab
demo_phaseCompensateMIMO
```

会显示 6 个子图：

- `True phase`：合成数据中注入的真实 RX 相位偏差。
- `Estimated phase`：算法估计出的 RX 相位偏差。
- `Residual phase error`：估计相位与真实相位之间的残差。
- `Channel coherence`：补偿前后各 RX 通道相对参考 RX 的相干性。
- `Coherent sum before`：补偿前跨 RX/TX/频率的相干叠加幅度。
- `Coherent sum after`：补偿后跨 RX/TX/频率的相干叠加幅度。

图像类子图的横轴为 `Scan x index`，表示二维扫描平面的 x 方向采样序号；纵轴为 `Scan y index`，表示 y 方向采样序号。柱状图横轴为 RX 通道编号，纵轴为相对参考 RX 的归一化相干性。

脚本还会额外生成一张校正前/后成像效果对比图：

- `Image before compensation (dB)`：相位补偿前的二维成像结果。
- `Image after compensation (dB)`：相位补偿后的二维成像结果。
- `Image gain after compensation (dB)`：补偿后图像相对补偿前的 dB 增益分布。

这里的示例成像使用跨 RX、TX 和选定频点的复数相干叠加：

```matlab
imgBeforeComplex = squeeze(sum(sum(sum(dataJittered(:, :, opts.freqIdx, :, :), 1), 2), 3));
imgAfterComplex = squeeze(sum(sum(sum(dataCorr(:, :, opts.freqIdx, :, :), 1), 2), 3));
```

随后将幅度图归一化到最大值并转成 dB 显示。脚本同时打印 `Image contrast before/after`，该指标定义为图像幅度标准差除以均值，用于粗略评价补偿后图像聚焦和对比度是否提升。

## 适用条件与限制

PGA-like 方法适合估计随扫描位置变化的相位误差。由于它从相位梯度积分得到相位，因此不能唯一恢复每个 RX 的绝对常量相位；若需要补偿固定通道相位，可使用 `opts.method='coherence'` 或单独做通道标定。

该方法仍假设相位误差主要表现为“每个 RX 通道、每个扫描位置”的频率近似无关相位偏差。若抖动等效为明显的距离变化，导致相位误差随频率线性变化，则需要扩展为距离/频率相关相位补偿模型。
