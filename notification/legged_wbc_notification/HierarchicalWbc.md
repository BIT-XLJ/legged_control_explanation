# 结构：

* **命名空间 legged**
  * **类 HierarchicalWbc**：继承自基类WbcBase
    * **update函数：**返回值为vector_t类型的向量，重新实现了WbcBase基类中的虚函数update
    * 部分展开了WbcBase类，暂不清楚目的

# 实现：

* **vector_t HierarchicalWbc::update**函数

  先运行基类WbcBase中的虚函数update（注意，不是纯虚函数，WbcBase::update函数非空），然后整合三个任务。                                                                             

  * Task0：浮动基座Eom任务+满足力矩限制方程的任务+满足摩擦锥方程的任务+满足非接触运动方程的任务

  * Task1：满足基座加速度方程的任务+满足摆动足方程的任务

  * Task2：满足接触力方程的任务

  * 数字越低，优先级越高

  最后按照任务优先级顺序调用HoQp::hoQp函数解算分级WBC控制的输出，并返回解向量

  