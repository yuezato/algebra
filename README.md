# algebra
Rustで簡単な代数のライブラリを作る

現在実装しているのは次の通り
+ 一般の体
+ 有限体（ガロア体, GF(2), GF(2^16))
+ 係数が体となる1変数多項式環（univariate polynomial)
+ 体を要素とする行列

未実装だがやることは次の通り
+ Vandermonde行列ベースのErasure Coding
+ https://dl.acm.org/doi/10.14778/2535573.2488339
+ https://ieeexplore.ieee.org/abstract/document/8340062
+ https://dl.acm.org/doi/abs/10.1145/2897518.2897525
