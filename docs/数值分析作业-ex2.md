# 数值分析作业

## 题目

* 题目1   
* 题目2   
* 题目3   

* 习题二    第2题   完成
* 习题二    第3题   完成
* 习题二    第4题   完成
* 数值实验二   第一题   完成

### 题目1

$$
$\left(\begin{array}{ccccc}
2 & 0 & 0 & 0 & 0\\
-\frac{1}{2} & \frac{\sqrt{15}}{2} & 0 & 0 & 0\\
0 & -\frac{2\,\sqrt{15}}{15} & \frac{2\,\sqrt{14}\,\sqrt{15}}{15} & 0 & 0\\
0 & 0 & -\frac{\sqrt{14}\,\sqrt{15}}{28} & \frac{\sqrt{14}\,\sqrt{209}}{28} & 0\\
0 & 0 & 0 & -\frac{2\,\sqrt{14}\,\sqrt{209}}{209} & \frac{2\,\sqrt{195}\,\sqrt{209}}{209}
\end{array}\right)$
$$



```matlab
初始矩阵 A:
     4    -1     0     0     0
    -1     4    -1     0     0
     0    -1     4    -1     0
     0     0    -1     4    -1
     0     0     0    -1     4


计算 L(1, 1): 2.0000
     2     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0

计算 L(2, 1): -0.5000
    2.0000         0         0         0         0
   -0.5000         0         0         0         0
         0         0         0         0         0
         0         0         0         0         0
         0         0         0         0         0

计算 L(3, 1): 0.0000
    2.0000         0         0         0         0
   -0.5000         0         0         0         0
         0         0         0         0         0
         0         0         0         0         0
         0         0         0         0         0

计算 L(4, 1): 0.0000
    2.0000         0         0         0         0
   -0.5000         0         0         0         0
         0         0         0         0         0
         0         0         0         0         0
         0         0         0         0         0

计算 L(5, 1): 0.0000
    2.0000         0         0         0         0
   -0.5000         0         0         0         0
         0         0         0         0         0
         0         0         0         0         0
         0         0         0         0         0


计算 L(2, 2): 1.9365
    2.0000         0         0         0         0
   -0.5000    1.9365         0         0         0
         0         0         0         0         0
         0         0         0         0         0
         0         0         0         0         0

计算 L(3, 2): -0.5164
    2.0000         0         0         0         0
   -0.5000    1.9365         0         0         0
         0   -0.5164         0         0         0
         0         0         0         0         0
         0         0         0         0         0

计算 L(4, 2): 0.0000
    2.0000         0         0         0         0
   -0.5000    1.9365         0         0         0
         0   -0.5164         0         0         0
         0         0         0         0         0
         0         0         0         0         0

计算 L(5, 2): 0.0000
    2.0000         0         0         0         0
   -0.5000    1.9365         0         0         0
         0   -0.5164         0         0         0
         0         0         0         0         0
         0         0         0         0         0


计算 L(3, 3): 1.9322
    2.0000         0         0         0         0
   -0.5000    1.9365         0         0         0
         0   -0.5164    1.9322         0         0
         0         0         0         0         0
         0         0         0         0         0

计算 L(4, 3): -0.5175
    2.0000         0         0         0         0
   -0.5000    1.9365         0         0         0
         0   -0.5164    1.9322         0         0
         0         0   -0.5175         0         0
         0         0         0         0         0

计算 L(5, 3): 0.0000
    2.0000         0         0         0         0
   -0.5000    1.9365         0         0         0
         0   -0.5164    1.9322         0         0
         0         0   -0.5175         0         0
         0         0         0         0         0


计算 L(4, 4): 1.9319
    2.0000         0         0         0         0
   -0.5000    1.9365         0         0         0
         0   -0.5164    1.9322         0         0
         0         0   -0.5175    1.9319         0
         0         0         0         0         0

计算 L(5, 4): -0.5176
    2.0000         0         0         0         0
   -0.5000    1.9365         0         0         0
         0   -0.5164    1.9322         0         0
         0         0   -0.5175    1.9319         0
         0         0         0   -0.5176         0


计算 L(5, 5): 1.9319
    2.0000         0         0         0         0
   -0.5000    1.9365         0         0         0
         0   -0.5164    1.9322         0         0
         0         0   -0.5175    1.9319         0
         0         0         0   -0.5176    1.9319


最终的 L 矩阵:
    2.0000         0         0         0         0
   -0.5000    1.9365         0         0         0
         0   -0.5164    1.9322         0         0
         0         0   -0.5175    1.9319         0
         0         0         0   -0.5176    1.9319


前代法解 Ly = b 后的 y 向量:
    5.0000
   11.6190
   13.4563
   13.9576
    8.9162


回代法解 L'x = y 后的 x 向量 (方程组解):
    4.6154
    8.4615
    9.2308
    8.4615
    4.6154
```



### 题目2

#### 1 高斯消去法

```matlab
function x = gaussian_elimination(H, d)
    % 高斯消去法 (不使用主元)
    n = length(d);
    A = [H, d];
    
    % 消元过程
    for i = 1:n-1
        for j = i+1:n
            m = A(j,i) / A(i,i);
            A(j,:) = A(j,:) - m * A(i,:);
        end
    end
    
    % 回代过程
    x = zeros(n,1);
    x(n) = A(n,n+1) / A(n,n);
    for i = n-1:-1:1
        x(i) = (A(i,n+1) - A(i,i+1:n) * x(i+1:n)) / A(i,i);
    end
end

```



#### 2 列主元高斯消去法

```matlab
function x = gaussian_elimination_pivot(H, d)
    % 列主元高斯消去法
    n = length(d);
    A = [H, d];
    
    for i = 1:n-1
        % 找到列中的最大元素并交换行
        [~, max_idx] = max(abs(A(i:n,i)));
        max_idx = max_idx + i - 1;
        if max_idx ~= i
            A([i max_idx], :) = A([max_idx i], :);
        end
        
        % 正常的高斯消去
        for j = i+1:n
            m = A(j,i) / A(i,i);
            A(j,:) = A(j,:) - m * A(i,:);
        end
    end
    
    % 回代过程
    x = zeros(n,1);
    x(n) = A(n,n+1) / A(n,n);
    for i = n-1:-1:1
        x(i) = (A(i,n+1) - A(i,i+1:n) * x(i+1:n)) / A(i,i);
    end
end

```



#### 3 全主元高斯消去法

```matlab
function x = gaussian_elimination_full_pivot(H, d)
    % 全主元高斯消去法
    n = length(d);
    A = [H, d];
    p = 1:n;
    
    for i = 1:n-1
        % 找到矩阵中绝对值最大的元素作为主元
        [~, max_idx] = max(abs(A(i:n,i:n)), [], 'all', 'linear');
        [max_row, max_col] = ind2sub([n-i+1, n-i+1], max_idx);
        max_row = max_row + i - 1;
        max_col = max_col + i - 1;
        
        % 行交换
        if max_row ~= i
            A([i, max_row], :) = A([max_row, i], :);
        end
        
        % 列交换
        if max_col ~= i
            A(:, [i, max_col]) = A(:, [max_col, i]);
            p([i, max_col]) = p([max_col, i]);
        end
        
        % 正常的高斯消去
        for j = i+1:n
            m = A(j,i) / A(i,i);
            A(j,:) = A(j,:) - m * A(i,:);
        end
    end
    
    % 回代过程
    x = zeros(n,1);
    x(n) = A(n,n+1) / A(n,n);
    for i = n-1:-1:1
        x(i) = (A(i,n+1) - A(i,i+1:n) * x(i+1:n)) / A(i,i);
    end
    
    % 恢复列交换后的顺序
    x(p) = x;
end

```



#### 4 平衡加权高斯消去法

```matlab
function x = balanced_weight_gaussian_elimination(H, d)
    % 平衡加权高斯消去法
    n = length(d);
    A = [H, d];
    
    for i = 1:n-1
        % 计算平衡因子
        w = max(abs(A(i:n,i:n)), [], 2) ./ max(abs(A(i:n,:)), [], 2);
        [~, max_idx] = max(w);
        max_idx = max_idx + i - 1;
        
        % 行交换
        if max_idx ~= i
            A([i, max_idx], :) = A([max_idx, i], :);
        end
        
        % 正常的高斯消去
        for j = i+1:n
            m = A(j,i) / A(i,i);
            A(j,:) = A(j,:) - m * A(i,:);
        end
    end
    
    % 回代过程
    x = zeros(n,1);
    x(n) = A(n,n+1) / A(n,n);
    for i = n-1:-1:1
        x(i) = (A(i,n+1) - A(i,i+1:n) * x(i+1:n)) / A(i,i);
    end
end

```

#### 5 主程序

```matlab
% Hilbert矩阵阶数
n_values = [5, 10, 20, 40];

% 方法列表，包含四种高斯消去法以及 MATLAB 的左除法
methods = {'gaussian_elimination', 'gaussian_elimination_pivot', 'gaussian_elimination_full_pivot', 'balanced_weight_gaussian_elimination', 'matlab_backslash'};

% 存储各方法的误差和时间
results = zeros(length(n_values), length(methods));
times = zeros(length(n_values), length(methods));

% 循环处理不同的矩阵阶数
for i = 1:length(n_values)
    n = n_values(i);
    H = hilb(n);  % 生成Hilbert矩阵
    d = (n + 1 - (1:n))';  % 定义右端项d

    % 使用不同的求解方法
    for j = 1:length(methods)
        if strcmp(methods{j}, 'matlab_backslash')
            % 如果是左除方法
            tic;  % 开始计时
            x = H \ d;  % 直接用左除解方程
            elapsed_time = toc;  % 结束计时
        else
            % 调用自定义的求解方法
            tic;
            x = feval(methods{j}, H, d);
            elapsed_time = toc;
        end
        
        % 记录时间消耗
        times(i, j) = elapsed_time;
        fprintf('n = %d, Method: %s, Time: %.4f seconds\n', n, methods{j}, elapsed_time);
        
        % 将解与 MATLAB 的左除结果进行比较
        x_exact = H \ d;  % 通过左除法得到精确解
        error_norm = norm(x - x_exact);  % 计算误差
        fprintf('Error norm compared to left division: %.4e\n', error_norm);
        
        % 记录误差
        results(i, j) = error_norm;
    end
end

% 显示最终结果
disp('误差对比：');
disp(results);
disp('时间消耗对比：');
disp(times);

```



#### 6 输出结果

> ```
> n = 5, Method: gaussian_elimination, Time: 0.0028 seconds
> Error norm compared to left division: 1.5996e-08
> ```
>
> ```
> n = 5, Method: gaussian_elimination_pivot, Time: 0.0036 seconds
> Error norm compared to left division: 2.2737e-12
> ```
>
> ```
> n = 5, Method: gaussian_elimination_full_pivot, Time: 0.0207 seconds
> Error norm compared to left division: 2.7330e-08
> ```
>
> ```
> n = 5, Method: balanced_weight_gaussian_elimination, Time: 0.0040 seconds
> Error norm compared to left division: 4.4169e-09
> ```
>
> ```
> n = 5, Method: matlab_backslash, Time: 0.0001 seconds
> Error norm compared to left division: 0.0000e+00
> ```
>
> 
>
> ```
> n = 10, Method: gaussian_elimination, Time: 0.0013 seconds 
> Error norm compared to left division: 3.5218e+04
> ```
>
> ```
> n = 10, Method: gaussian_elimination_pivot, Time: 0.0021 seconds
> Error norm compared to left division: 1.8109e-06
> ```
>
> ```
> n = 10, Method: gaussian_elimination_full_pivot, Time: 0.0038 seconds
> Error norm compared to left division: 1.8568e+02
> ```
>
> ```
> n = 10, Method: balanced_weight_gaussian_elimination, Time: 0.0022 seconds
> Error norm compared to left division: 1.2116e+03
> ```
>
> ```
> n = 10, Method: matlab_backslash, Time: 0.0000 seconds
> Error norm compared to left division: 0.0000e+00
> ```
>
> 
>
> ```
> n = 20, Method: gaussian_elimination, Time: 0.0001 seconds
> 警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  1.276108e-19。
> Error norm compared to left division: 1.4478e+12
> ```
>
> ```
> n = 20, Method: gaussian_elimination_pivot, Time: 0.0004 seconds
> 警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  1.276108e-19。
> Error norm compared to left division: 8.5371e+11
> ```
>
> ```
> n = 20, Method: gaussian_elimination_full_pivot, Time: 0.0003 seconds
> 警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  1.276108e-19。
> Error norm compared to left division: 3.9477e+12
> ```
>
> ```
> n = 20, Method: balanced_weight_gaussian_elimination, Time: 0.0003 seconds
> 警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  1.276108e-19。
> Error norm compared to left division: 6.8510e+11
> ```
>
> ```
> n = 20, Method: matlab_backslash, Time: 0.0002 seconds
> 
> 警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  1.276108e-19。
> Error norm compared to left division: 0.0000e+00
> ```
>
> 
>
> ```
> n = 40, Method: gaussian_elimination, Time: 0.0003 seconds
> 警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  7.520353e-20。
> Error norm compared to left division: 6.8838e+12
> ```
>
> ```
> n = 40, Method: gaussian_elimination_pivot, Time: 0.0006 seconds
> 警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  7.520353e-20。
> Error norm compared to left division: 7.3147e+12
> ```
>
> ```
> n = 40, Method: gaussian_elimination_full_pivot, Time: 0.0007 seconds
> 警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  7.520353e-20。
> Error norm compared to left division: 5.3101e+12
> ```
>
> ```
> n = 40, Method: balanced_weight_gaussian_elimination, Time: 0.0010 seconds
> 警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  7.520353e-20。
> Error norm compared to left division: 4.7199e+12
> ```
>
> ```
> n = 40, Method: matlab_backslash, Time: 0.0003 seconds
> 警告: 矩阵接近奇异值，或者缩放不良。结果可能不准确。RCOND =  7.520353e-20。
> Error norm compared to left division: 0.0000e+00
> ```
>
> 
>
> ```
> 误差对比：
> ```
>
> ```
>    1.0e+12 *
> 
>     0.0000    0.0000    0.0000    0.0000         0
>     0.0000    0.0000    0.0000    0.0000         0
>     1.4478    0.8537    3.9477    0.6851         0
>     6.8838    7.3147    5.3101    4.7199         0
> ```
>
> ```
> 时间消耗对比：
> ```
>
> ```
>     0.0028    0.0036    0.0207    0.0040    0.0001
>     0.0013    0.0021    0.0038    0.0022    0.0000
>     0.0001    0.0004    0.0003    0.0003    0.0002
>     0.0003    0.0006    0.0007    0.0010    0.0003
> ```

### 题目3

```matlab
function solve_Ax_b(n)
    % 设置矩阵 B
    e = ones(n, 1);
    B = spdiags([-e -4*e -e], -1:1, n, n);  % 构造五对角矩阵 B
    
    % 单位矩阵 I
    I = speye(n);  % 稀疏单位矩阵
    
    % 构造矩阵 A
    A = kron(speye(n), B) + kron(spdiags(-e, -1, n, n), I) + kron(spdiags(-e, 1, n, n), I);
    
    % 构造向量 b
    b = (0:n^2-1)';
    
    % 求解方程 Ax = b
    x = A \ b;
    
    % 显示结果
    disp('The solution x is:');
    disp(x);
end

% 程序求解
n = 80;
solve_Ax_b(80);  
```



### 习题二 第2题

```matlab
初始增广矩阵:
     1    13    -2   -34    13
     2     6    -7   -10   -22
   -10    -1     5     9    14
    -3    -5     0    15   -36


第 1 列主元消元:
交换第 1 行和第 3 行:
   -10    -1     5     9    14
     2     6    -7   -10   -22
     1    13    -2   -34    13
    -3    -5     0    15   -36

消去第 2 行第 1 列，消元因子为 -0.2000:
  -10.0000   -1.0000    5.0000    9.0000   14.0000
         0    5.8000   -6.0000   -8.2000  -19.2000
    1.0000   13.0000   -2.0000  -34.0000   13.0000
   -3.0000   -5.0000         0   15.0000  -36.0000

消去第 3 行第 1 列，消元因子为 -0.1000:
  -10.0000   -1.0000    5.0000    9.0000   14.0000
         0    5.8000   -6.0000   -8.2000  -19.2000
         0   12.9000   -1.5000  -33.1000   14.4000
   -3.0000   -5.0000         0   15.0000  -36.0000

消去第 4 行第 1 列，消元因子为 0.3000:
  -10.0000   -1.0000    5.0000    9.0000   14.0000
         0    5.8000   -6.0000   -8.2000  -19.2000
         0   12.9000   -1.5000  -33.1000   14.4000
         0   -4.7000   -1.5000   12.3000  -40.2000


第 2 列主元消元:
交换第 2 行和第 3 行:
  -10.0000   -1.0000    5.0000    9.0000   14.0000
         0   12.9000   -1.5000  -33.1000   14.4000
         0    5.8000   -6.0000   -8.2000  -19.2000
         0   -4.7000   -1.5000   12.3000  -40.2000

消去第 3 行第 2 列，消元因子为 0.4496:
  -10.0000   -1.0000    5.0000    9.0000   14.0000
         0   12.9000   -1.5000  -33.1000   14.4000
         0         0   -5.3256    6.6822  -25.6744
         0   -4.7000   -1.5000   12.3000  -40.2000

消去第 4 行第 2 列，消元因子为 -0.3643:
  -10.0000   -1.0000    5.0000    9.0000   14.0000
         0   12.9000   -1.5000  -33.1000   14.4000
         0         0   -5.3256    6.6822  -25.6744
         0         0   -2.0465    0.2403  -34.9535


第 3 列主元消元:
消去第 4 行第 3 列，消元因子为 0.3843:
  -10.0000   -1.0000    5.0000    9.0000   14.0000
         0   12.9000   -1.5000  -33.1000   14.4000
         0         0   -5.3256    6.6822  -25.6744
         0         0         0   -2.3275  -25.0873


经过高斯消元后的上三角增广矩阵:
  -10.0000   -1.0000    5.0000    9.0000   14.0000
         0   12.9000   -1.5000  -33.1000   14.4000
         0         0   -5.3256    6.6822  -25.6744
         0         0         0   -2.3275  -25.0873


回代后的解向量:
   14.3827
   30.9062
   18.3452
   10.7786
```





### 习题二    第3题

```matlab
初始矩阵 A:
    15     7     0    10
     6    18    15     9
     0    10    28     7
     5     0     6    35


计算 U 的第 1 行:
    15     7     0    10
     0     0     0     0
     0     0     0     0
     0     0     0     0


计算 L 的第 1 列:
    1.0000         0         0         0
    0.4000    1.0000         0         0
         0         0    1.0000         0
    0.3333         0         0    1.0000


计算 U 的第 2 行:
   15.0000    7.0000         0   10.0000
         0   15.2000   15.0000    5.0000
         0         0         0         0
         0         0         0         0


计算 L 的第 2 列:
    1.0000         0         0         0
    0.4000    1.0000         0         0
         0    0.6579    1.0000         0
    0.3333   -0.1535         0    1.0000


计算 U 的第 3 行:
   15.0000    7.0000         0   10.0000
         0   15.2000   15.0000    5.0000
         0         0   18.1316    3.7105
         0         0         0         0


计算 L 的第 3 列:
    1.0000         0         0         0
    0.4000    1.0000         0         0
         0    0.6579    1.0000         0
    0.3333   -0.1535    0.4579    1.0000


计算 U 的第 4 行:
   15.0000    7.0000         0   10.0000
         0   15.2000   15.0000    5.0000
         0         0   18.1316    3.7105
         0         0         0   30.7351


计算 L 的第 4 列:
    1.0000         0         0         0
    0.4000    1.0000         0         0
         0    0.6579    1.0000         0
    0.3333   -0.1535    0.4579    1.0000


最终的 L 矩阵:
    1.0000         0         0         0
    0.4000    1.0000         0         0
         0    0.6579    1.0000         0
    0.3333   -0.1535    0.4579    1.0000

最终的 U 矩阵:
   15.0000    7.0000         0   10.0000
         0   15.2000   15.0000    5.0000
         0         0   18.1316    3.7105
         0         0         0   30.7351


前代法解 Ly = b 后的 y 向量:
    8.0000
    2.8000
    2.1579
   -1.2250


回代法解 Ux = y 后的 x 向量 (方程组解):
    0.5264
    0.0718
    0.1272
   -0.0399
```



### 习题二  第4题

```matlab
初始矩阵 A:
     4    -1     0     0
    -1     4    -1     0
     0    -1     4    -1
     0     0    -1     4


计算 L(1, 1): 2.0000
     2     0     0     0
     0     0     0     0
     0     0     0     0
     0     0     0     0

计算 L(2, 1): -0.5000
    2.0000         0         0         0
   -0.5000         0         0         0
         0         0         0         0
         0         0         0         0

计算 L(3, 1): 0.0000
    2.0000         0         0         0
   -0.5000         0         0         0
         0         0         0         0
         0         0         0         0

计算 L(4, 1): 0.0000
    2.0000         0         0         0
   -0.5000         0         0         0
         0         0         0         0
         0         0         0         0


计算 L(2, 2): 1.9365
    2.0000         0         0         0
   -0.5000    1.9365         0         0
         0         0         0         0
         0         0         0         0

计算 L(3, 2): -0.5164
    2.0000         0         0         0
   -0.5000    1.9365         0         0
         0   -0.5164         0         0
         0         0         0         0

计算 L(4, 2): 0.0000
    2.0000         0         0         0
   -0.5000    1.9365         0         0
         0   -0.5164         0         0
         0         0         0         0


计算 L(3, 3): 1.9322
    2.0000         0         0         0
   -0.5000    1.9365         0         0
         0   -0.5164    1.9322         0
         0         0         0         0

计算 L(4, 3): -0.5175
    2.0000         0         0         0
   -0.5000    1.9365         0         0
         0   -0.5164    1.9322         0
         0         0   -0.5175         0


计算 L(4, 4): 1.9319
    2.0000         0         0         0
   -0.5000    1.9365         0         0
         0   -0.5164    1.9322         0
         0         0   -0.5175    1.9319


最终的 L 矩阵:
    2.0000         0         0         0
   -0.5000    1.9365         0         0
         0   -0.5164    1.9322         0
         0         0   -0.5175    1.9319


前代法解 Ly = b 后的 y 向量:
    1.0000
    2.3238
    6.3141
   -1.9319


回代法解 L'x = y 后的 x 向量 (方程组解):
    1.0000
    2.0000
    3.0000
   -1.0000
```



### 数值实验  第1题

```matlab
function x = tridiagonal_solve(a, b, c, d)
    % 输入: 
    % a 是下对角线向量 (长度 n-1), 其中 a(1) = 0
    % b 是主对角线向量 (长度 n)
    % c 是上对角线向量 (长度 n-1), 其中 c(n) = 0
    % d 是常数向量 (长度 n)
    % 输出: 方程组的解 x (长度 n)

    n = length(b);  % 方程组大小
    x = zeros(n, 1);  % 初始化解向量
    % 前向消去
    for i = 2:n
        w = a(i-1) / b(i-1);
        b(i) = b(i) - w * c(i-1);
        d(i) = d(i) - w * d(i-1);
    end
    
    % 回代
    x(n) = d(n) / b(n);
    for i = n-1:-1:1
        x(i) = (d(i) - c(i) * x(i+1)) / b(i);
    end

    % 输出结果
    fprintf('解向量 x:\n');
    disp(x);
end

% 程序求解
n = 101;  

a = ones(n-1, 1);  
b = 12 * ones(n, 1);  
c = ones(n-1, 1);  
d = [11; 10 * ones(n-2, 1); 11];  

x = tridiagonal_solve(a, b, c, d);
```

