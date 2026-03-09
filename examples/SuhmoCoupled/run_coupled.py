#!/usr/bin/env python3
"""
Script-based orchestration for BISICLES-SUHMO coupled simulation.

Each model is restarted from checkpoint at each coupling interval.
The workflow is:

  BISICLES (N steps) → write H, u for SUHMO → write checkpoint
  SUHMO (to steady state) → read H, u → write Pw → write checkpoint
  BISICLES (N steps, restart from chk) → read Pw → write H, u → ...

Usage:
    python run_coupled.py [--nprocs 4] [--coupling-steps 10] [--total-steps 100]
"""
import subprocess
import os
import glob
import argparse
import sys


def find_latest_checkpoint(prefix, dim="2d"):
    """Find the highest-numbered checkpoint file matching prefix.NNNNNN.2d.hdf5"""
    pattern = f"{prefix}*.{dim}.hdf5"
    files = sorted(glob.glob(pattern))
    if not files:
        return None
    return files[-1]


def find_latest_plotfile(prefix, dim="2d"):
    """Find the highest-numbered plot file matching prefix.NNNNNN.2d.hdf5"""
    pattern = f"{prefix}*.{dim}.hdf5"
    files = sorted(glob.glob(pattern))
    if not files:
        return None
    return files[-1]


def run_cmd(cmd, label):
    """Run a shell command, printing label and checking return code."""
    print(f"\n{'='*60}")
    print(f"[COUPLED] {label}")
    print(f"  cmd: {cmd}")
    print(f"{'='*60}\n")
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        print(f"[COUPLED] ERROR: {label} exited with code {result.returncode}")
        sys.exit(result.returncode)


def run_bisicles(args, coupling_iter, restart_chk=None):
    """
    Run BISICLES for coupling_steps timesteps.
    On first call, starts fresh. On subsequent calls, restarts from checkpoint.
    Always writes a coupling output file for SUHMO and a checkpoint for itself.
    """
    input_file = args.bisicles_input

    # Build command with overrides via ParmParse command-line syntax
    cmd = f"mpirun -n {args.nprocs} {args.bisicles_exe} {input_file}"

    # Override max steps for this coupling window
    cmd += f" amr.maxStep={args.coupling_steps * (coupling_iter + 1)}"

    # Restart from checkpoint if available
    if restart_chk is not None:
        cmd += f" amr.restart_file={restart_chk}"

    # Tell BISICLES to write coupling output
    cmd += f" amr.writeSUHMOCoupling=true"
    cmd += f" amr.suhmoOutputFile={args.b2s_file}"

    # Force a checkpoint at the end
    cmd += f" amr.check_interval={args.coupling_steps}"

    # Point to SUHMO's water pressure file for reading N
    cmd += f" SUHMOCoupling.waterPressureFile={args.s2b_file}"

    run_cmd(cmd, f"BISICLES coupling iteration {coupling_iter}")


def run_suhmo(args, coupling_iter):
    """
    Run SUHMO to steady state, reading thickness/velocity from the
    BISICLES coupling file. Always starts fresh (no restart) since
    we want a steady-state solve for the current geometry.
    """
    input_file = args.suhmo_input

    cmd = f"mpirun -n {args.nprocs} {args.suhmo_exe} {input_file}"

    # Override to read from BISICLES file
    cmd += f" suhmo.readFromBISICLES=true"
    cmd += f" suhmo.bisiclesFile={args.b2s_file}"

    run_cmd(cmd, f"SUHMO steady-state solve, iteration {coupling_iter}")

    # After SUHMO finishes, its final plot file contains Pw.
    # Find it and copy/rename it to the standard coupling filename.
    latest_plot = find_latest_plotfile("plot")
    if latest_plot is None:
        print("[COUPLED] ERROR: no SUHMO plot file found after run")
        sys.exit(1)

    # Create a symlink or copy so BISICLES can find it at a fixed name
    if os.path.exists(args.s2b_file):
        os.remove(args.s2b_file)
    os.symlink(latest_plot, args.s2b_file)
    print(f"[COUPLED] Linked {latest_plot} -> {args.s2b_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Orchestrate BISICLES-SUHMO coupled simulation")

    parser.add_argument("--nprocs", type=int, default=4,
                        help="Number of MPI processes")
    parser.add_argument("--coupling-steps", type=int, default=10,
                        help="BISICLES timesteps per coupling interval")
    parser.add_argument("--num-coupling-iters", type=int, default=10,
                        help="Total number of coupling iterations")

    parser.add_argument("--bisicles-exe", type=str,
                        default="./driver2d.Linux.64.mpicxx.gfortran.OPT.MPI.ex",
                        help="Path to BISICLES executable")
    parser.add_argument("--bisicles-input", type=str,
                        default="inputs.bisicles",
                        help="BISICLES input file")

    parser.add_argument("--suhmo-exe", type=str,
                        default="./Suhmo2d.Linux.64.mpicxx.gfortran.OPT.MPI.ex",
                        help="Path to SUHMO executable")
    parser.add_argument("--suhmo-input", type=str,
                        default="input.hydro",
                        help="SUHMO input file")

    parser.add_argument("--b2s-file", type=str,
                        default="bisicles_for_suhmo.2d.hdf5",
                        help="BISICLES→SUHMO coupling file")
    parser.add_argument("--s2b-file", type=str,
                        default="suhmo_for_bisicles.2d.hdf5",
                        help="SUHMO→BISICLES coupling file")

    args = parser.parse_args()

    bisicles_chk = None

    for iteration in range(args.num_coupling_iters):
        print(f"\n{'#'*60}")
        print(f"# COUPLING ITERATION {iteration}")
        print(f"{'#'*60}")

        # ─── Run BISICLES ───
        run_bisicles(args, iteration, restart_chk=bisicles_chk)

        # Find the checkpoint BISICLES just wrote for next restart
        bisicles_chk = find_latest_checkpoint("chk")
        if bisicles_chk is None and iteration < args.num_coupling_iters - 1:
            print("[COUPLED] WARNING: no BISICLES checkpoint found for restart")

        # ─── Run SUHMO to steady state ───
        run_suhmo(args, iteration)

    print(f"\n{'#'*60}")
    print(f"# COUPLED SIMULATION COMPLETE")
    print(f"{'#'*60}")


if __name__ == "__main__":
    main()