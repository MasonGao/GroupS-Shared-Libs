<div align="center">

#  C++ Shared Library 

**The goal of this project is to manage generic C++ code of Group'S.**

<br />

</div>




## Introduction

This project aims to provide a comprehensive solution for managing and organizing generic C++ code in HEP calculation. It is designed to simplify the calculation process, improve code maintainability, and enhance collaboration among students in Group'S. 

## Features

### 1. Code Organization

The project offers a structured approach to organizing C++ code. It provides a directory structure that separates different components of the program, such as **constants**, **anomalous dimension**, and **hard functions**. This organization system helps in locating and managing code files efficiently. Such as QCD Beta function:

$$
\beta(\alpha_{s})=-2\alpha_{s}\sum_{n=0}^{\infty}\beta_{n}\left(\frac{\alpha_{s}}{4\pi}\right)^{n+1},
$$

```
double Beta[4] = {(11 * CA) / 3. - (4 * nf * TF) / 3., (34 * pow(CA, 2)) / 3. - (20 * CA * nf * TF) / 3. - 4 * CF *nf *TF, (2857 * pow(CA, 3)) / 54. + ((-1415 * pow(CA, 2)) / 27. - (205 * CA * CF) / 9. + 2 * pow(CF, 2)) * nf *TF + ((158 * CA) / 27. + (44 * CF) / 9.) * pow(nf, 2) * pow(TF, 2), 149753 / 6. + (1093 * pow(nf, 3)) / 729. + 3564 * riemann_zeta(3) + pow(nf, 2) * (50065 / 162. + (6472 * riemann_zeta(3)) / 81.) - nf * (1078361 / 162. + (6508 * riemann_zeta(3)) / 27.)};

```



### 2. Build System

The project includes a reliable build system that automates the compilation and linking process. It supports various build configurations, allowing developers to build the code for different platforms and environments effortlessly.


### 3. Code Formatting and Style Guidelines

Consistency in coding style is crucial for readability and maintainability. The project enforces code formatting rules and adheres to widely accepted style guidelines, making the codebase consistent and easier to understand.
1. Indentation: Use two spaces for indentation instead of tabs.
2. Line Length: Keep each line of code within 80 characters or fewer, and break lines that exceed this limit.
3. Braces: Place the opening brace on its own line, aligned with the control structure, and align the closing brace with the ending position of the control structure.
4. Variable Naming: Use lowercase letters and underscores for variable names, e.g., "my_variable".
5. Type Naming: Use camel case for type names, e.g., "MyClass".
6. Constant Naming: Use uppercase letters and underscores for constant names, e.g., "CONSTANT_VALUE".
7. Function Naming: Use camel case for function names, starting with a lowercase letter, e.g., "myFunction".
8. Comments: Use C++-style double-slash comments "//" and follow specific comment formatting guidelines.
9. File Header Comments: Each source file should include copyright information and a brief description in a file header comment.
10. Operator Spacing: Leave appropriate spaces around operators to improve readability.

### 4. Dependency Management

This project relies on the following C++ libraries or programs: 


| [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)  | [LHAPDF](https://lhapdf.hepforge.org/)  | [CUBA](http://www.feynarts.de/cuba/)   |
| ------- | ------- |------- |




## Call to Action

Whether you have a particular passion about a specific language or framework, want to help in creating a more robust codeblock, or generally have interesting ideas on how to make this better, we'd love to have you!




