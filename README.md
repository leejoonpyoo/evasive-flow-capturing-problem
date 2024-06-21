# Evasive Flow Capturing Problem

This project implements three models (MRS, SM, and PM) to solve the Evasive Flow Capturing Problem (EFCP) using the formulations provided in the paper by Okan Arslan et al. The primary objective is to strategically place law enforcement facilities in a 25-node network to intercept unlawful vehicle flows while minimizing the total setup and damage costs.

## Project Structure

- `functions.jl`: Contains functions to generate the 25-node network, find paths within tolerance, calculate path costs, and other utility functions.
- `model.jl`: Implements the three models (MRS, SM, and PM) using JuMP and Gurobi.
- `main.jl`: Main script to run each model and output the results.

## Getting Started

### Prerequisites

- Julia programming language
- JuMP package
- Gurobi package
- DataStructures package

### Installation

1. Install Julia from the [official website](https://julialang.org/).
2. Install the required Julia packages by running the following commands in the Julia REPL:

    ```julia
    using Pkg
    Pkg.add("JuMP")
    Pkg.add("Gurobi")
    Pkg.add("DataStructures")
    Pkg.add("Graphs")
    ```

3. Ensure you have Gurobi installed and properly configured. Follow the instructions on the [Gurobi website](https://www.gurobi.com/documentation/9.1/quickstart_mac/the_gurobi_command_line.html) for installation and licensing.

## Running the Models

1. Include the function and model files:

    ```julia
    include("functions.jl")
    include("model.jl")
    ```

2. Execute the `main.jl` script:

    ```julia
    include("main.jl")
    ```

## Data

The 25-node network data used in this project is presented in the paper by Okan Arslan et al. This network consists of 25 nodes and 86 arcs. The numbers on each arc represent distances, which are used to calculate costs and determine the allowable range for drivers seeking alternative paths.

## Explanation of Models

### MRS Model

The MRS (Markov, Rabani, and Schoenfeld) model addresses the EFCP by requiring the enumeration of all paths within a specified deviation tolerance value. The model successively solves \(k\)-shortest path problems, ensuring optimal placement of law enforcement facilities to minimize setup cost and damage from non-intercepted flows.

### SM Model

The SM (Single-stage Mixed-Integer Linear Programming) model is derived from the bilevel model by converting it into a single-stage equivalent using duality theory. The primary objective is to minimize the total cost of setting up new stations and the damage costs due to non-intercepted flows.

### PM Model

The PM (Projection Model) uses projection inequalities derived from the dual formulation of the SM model, transforming the problem into a more computationally efficient form. The objective remains to minimize the total setup and damage costs, ensuring flow interception or following the shortest non-intercepted path.

## Experimental Results

The experiments replicated the conditions described in the original paper, varying the deviation tolerance (\(\lambda\)), WIM installation cost (\(w_e\)), and cost per unit distance (\(c_f\)). The results demonstrated the effectiveness of the PM model in achieving significantly lower computation times compared to the MRS and SM models.

## Conclusion

- The PM model is the most efficient and practical approach for solving the EFCP in large-scale networks.
- Advanced optimization techniques and efficient solvers like Gurobi enhance computational performance.
- Selecting appropriate models and methods is crucial for addressing complex optimization problems effectively.

## Reference

This project is based on the paper: "Exact Solution of the Evasive Flow Capturing Problem" by Okan Arslan, Ola Jabali, and Gilbert Laporte, published in INFORMS in 2018. You can find the paper [here](https://pubsonline.informs.org/doi/10.1287/opre.2018.1756).

<!-- ## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. -->

## Acknowledgments

- Okan Arslan, Ola Jabali, and Gilbert Laporte for the foundational paper on the Evasive Flow Capturing Problem.
- The Julia and JuMP communities for their support and development of the necessary packages.

