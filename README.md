# Hamming distance-based Heuristic for Flexible Job-shop Scheduling Problems

This repository provides the implementation of a heuristic optimization method for the **Flexible Job-Shop Scheduling Problem (FJSSP)** based on **Hamming distance**. The method is built upon a **Mixed-Integer Linear Programming (MILP)** formulation, where job routing and machine sequencing decisions are encoded using binary variables.

## 🔍 Overview

FJSSP is an NP-hard extension of the classical Job Shop Scheduling Problem, where each operation can be processed on multiple machines. This flexibility increases solution space complexity and motivates the use of heuristics for large-scale instances.

This research proposes a hybrid approach:

- A MILP model to represent job routing and scheduling.
- A Hamming-distance-based heuristic to explore the neighborhood of feasible solutions.

The heuristic starts from an initial feasible schedule and iteratively explores neighboring solutions by modifying binary variables within a predefined Hamming radius.

## 📈 Key Features

- MILP formulation of FJSSP using binary decision variables.
- Hamming distance parameters (`H_min`, `H_max`) to control exploration granularity.
- Efficient solution improvement through controlled binary perturbation.
- Tested on benchmark instances from literature.

## 🏆 Results

- Achieved **ε-optimal solutions** with up to **22% reduction** in computation time.
- Provided **near-optimal intermediate solutions** up to **65% faster** than baseline methods.

## 🔧 Future Work

- Perform **sensitivity analysis** on Hamming distance parameters to optimize the trade-off between local search precision and global exploration.
- Develop **adaptive radius tuning** strategies to improve convergence and robustness.
- Extend the method by integrating with ongoing research on:
  - **Dynamic scheduling**
  - **Predictive maintenance**

## 📜 Citation
If you use this work, please cite the corresponding conference paper (link to be added upon publication).
