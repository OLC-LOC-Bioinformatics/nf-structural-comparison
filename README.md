# nf-structural-comparison

*This workflow is currently under development and has not yet been validated.*

A Nextflow DSL2 workflow for the structural comparison of bacterial genomic assemblies.
This workflow uses [SyRI](https://github.com/schneebergerlab/syri) to compare alignments between genomic assemblies to identify structural differences and SNVs.

## Overview

Presently, only the left side of the overview shown below is available for use (comparing a scaffold-level genome to a reference genome):

![Overview of nf-structural-comparison workflow](assets/nf-structural-comparison-workflow.png)

Future updates will allow users to choose between analyzing a query genome that is at scaffold-level (low contiguity) or analyzing a query genome that is at chromosome-level (high contiguity).