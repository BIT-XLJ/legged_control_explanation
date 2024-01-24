//
// Created by qiayuan on 2022/6/28.
//
//
// Ref: https://github.com/bernhardpg/quadruped_locomotion
//

#pragma once

#include "legged_wbc/Task.h"

#include <memory>

namespace legged {
// Hierarchical Optimization Quadratic Program
class HoQp {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  using HoQpPtr = std::shared_ptr<HoQp>;      //这里的using类似typedef,取个别名

  explicit HoQp(const Task& task) : HoQp(task, nullptr){};    //explicit防止隐式转换，但是这里是什么意思呢？看不懂。如果HoQp构造函数的higherProblem是个空指针，则调用HoQp(const Task& task)进行构造，且不能隐式转换

  HoQp(Task task, HoQpPtr higherProblem);       //构造函数

  matrix_t getStackedZMatrix() const { return stackedZ_; }

  Task getStackedTasks() const { return stackedTasks_; }

  vector_t getStackedSlackSolutions() const { return stackedSlackVars_; }

  vector_t getSolutions() const {
    vector_t x = xPrev_ + stackedZPrev_ * decisionVarsSolutions_;
    return x;
  }

  size_t getSlackedNumVars() const { return stackedTasks_.d_.rows(); }

 private:
  void initVars();
  void formulateProblem();
  void solveProblem();

  void buildHMatrix();
  void buildCVector();
  void buildDMatrix();
  void buildFVector();

  void buildZMatrix();
  void stackSlackSolutions();

  Task task_, stackedTasksPrev_, stackedTasks_;
  HoQpPtr higherProblem_;

  bool hasEqConstraints_{}, hasIneqConstraints_{};
  size_t numSlackVars_{}, numDecisionVars_{};
  matrix_t stackedZPrev_, stackedZ_;
  vector_t stackedSlackSolutionsPrev_, xPrev_;
  size_t numPrevSlackVars_{};

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> h_, d_;
  vector_t c_, f_;
  vector_t stackedSlackVars_, slackVarsSolutions_, decisionVarsSolutions_;

  // Convenience matrices that are used multiple times
  matrix_t eyeNvNv_;
  matrix_t zeroNvNx_;
};

}  // namespace legged
