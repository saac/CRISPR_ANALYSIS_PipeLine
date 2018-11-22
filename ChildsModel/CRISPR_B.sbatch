#!/bin/bash

#SBATCH --job-name=Smu1e-7
#SBATCH --output=Smu1e-7.out
#SBATCH --error=Smu1e-7.err
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

#./Crispr Dp mu T Pt Sp Rep
./Crispr $1 $2 $3 $4 $5 $6


# Dp = Initial Phage Strains
# mu = Phage Mutation Rate
# T = running simulation Time
# Pt = Number of protospacer.
# Sp = Number of spacers.
# Rep = Replicate.
