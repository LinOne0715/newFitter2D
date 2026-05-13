# C++/HLS Implementation and Verification of a 2D Track Fitting Algorithm

## Overview

This project implements and verifies a hardware-oriented 2D track fitting algorithm using C++ and Vitis HLS.

The goal is to translate an algorithmic C++ reference model into an HLS-compatible design while maintaining numerical accuracy and improving hardware performance. The design focuses on fixed-point arithmetic, lookup-table based computation, pipelining, and dataflow optimization.

Python scripts are used to compare the HLS output with reference results and analyze the numerical resolution of key track parameters such as `phi0` and `omega`.

This project was developed as part of my study on hardware acceleration, FPGA-based design flow, and numerical verification for high-throughput trigger algorithms.

## Motivation

Modern ASIC and FPGA design flows require not only correct RTL or HLS implementation, but also careful analysis of latency, initiation interval, resource usage, timing behavior, and numerical accuracy.

This project explores how a C++-based tracking algorithm can be converted into a hardware-oriented HLS implementation. The work includes algorithm restructuring, fixed-point conversion, LUT-based approximation, pipeline optimization, and verification against reference outputs.

The project is also a foundation for future AI-assisted EDA flow optimization, such as automated design-space exploration, performance prediction, and verification data analysis.

## Features

- C++ reference implementation of a 2D track fitting algorithm
- Vitis HLS implementation for hardware-oriented synthesis
- Fixed-point arithmetic design using `ap_fixed`
- Lookup-table based approximation for special functions
- Pipeline and dataflow optimization for improved throughput
- Comparison between C++ reference, HLS output, and verification data
- Python scripts for numerical error analysis and resolution plotting
- Evaluation of latency, initiation interval, DSP usage, LUT usage, and accuracy trade-offs

