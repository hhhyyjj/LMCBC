# LMCBC Boundary Condition for OpenFOAM

This repository provides the source code for the **LMCBC** implemented as an independent OpenFOAM boundary condition module. The LMCBC boundary condition is designed for **APE-2-based solvers** to improve the non-reflecting performance of outlet boundaries in computational aeroacoustics simulations.

## 🔨 Compilation

You can compile the boundary condition using the standard OpenFOAM `wmake` utility:

```
wmake
```

## 🚀 Usage

An example of how to configure the LMCBC boundary condition is provided in the `tutorials/`.

**Note:** This boundary condition is designed to be used with APE-2-based solvers.
 The full solver used in the corresponding study is not open-sourced at this stage. The provided configuration template is intended to help users integrate the LMCBC boundary condition into their own solvers.

## 📚 Reference

If you use this boundary condition, please cite the corresponding paper (to be added upon publication).