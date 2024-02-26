# Model Formulation

$$
q=
\left[
\matrix{
_Ir_{IB}\\
q_{IB}\\
q_j
}
\right] \in SE(3) \cross \mathbb{R}^{n_j},
\space\space\space
u=
\left[
\matrix{
_IV_{IB}\\
_B\omega_{IB}\\
\dot{q}_j
}
\right
] \in \mathbb{R}^{n_u}
$$

​        其中，$_Ir_{IB}$是机身坐标系相对于世界系的位置向量在世界坐标系下的表示，$q_{IB}$是将机身坐标系中的向量映射到世界系下的单位四元数，$q_j$是关节位置向量，$_IV_{IB}$是机身速度在世界坐标系下的表示，$_B\omega_{IB}$是机身坐标系相对世界坐标系的旋转角速度在机身坐标系下的表示。$q$和$u$分别是足式机器人的广义位置向量和广义速度向量。

​       足式机器人系统动力学可以写成如下方程：
$$
M(q)\dot{u}+h(q,u)=S^T\tau+J_s^T\lambda
$$
其中，$M(q)$是系统广义质量矩阵，$h(q,u)$包含了科里奥利力，离心力和重力项。$S=[0_{n_{\tau}\cross(n_{u}-n_{\tau})}\space \space \space \mathbb{I}_{n_{\tau}\cross n_{\tau}}]$是选择矩阵，主要选择哪些自由度是驱动状态（浮动基机器人前六个自由度是虚拟关节）。$J_s \in \mathbb{R}^{3n_c \cross n_u}$是支撑雅可比，可以由触地足的接触雅可比堆叠组成，$J_s=[J_{C_{1}}^T \cdots J_{C_{n_{c}}}^T]^T$，其中$n_c$是触地足的数量，$\lambda$是约束力向量（触地足的所受的外力组成的向量）。如果接触足是点接触，那么每个接触约束都可由三个方程描述：
$$
_Ir_{IC}(t)=const\\
_I\dot{r}_{IC} = J_Cu=0\\
_I\ddot{r}_{IC}=J_C\dot{u}+\dot{J}_Cu=0
$$


