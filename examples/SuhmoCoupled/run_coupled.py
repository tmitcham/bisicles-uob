#!/usr/bin/env python3
"""
Orchestration script for iteratively coupled BISICLES–SUHMO simulation.

Coupling workflow per iteration
────────────────────────────────
  1. Run BISICLES for `coupling_steps` timesteps (restart from checkpoint
     on iterations > 0).  BISICLES writes regular plot files containing
     ice thickness and velocity (thickness, xVel, yVel).

  2. Symlink the latest BISICLES plot file to a fixed name so SUHMO can
     find it at a known path.

  3. Run SUHMO to steady state.  SUHMO reads the BISICLES geometry file,
     solves subglacial hydrology, and writes effective pressure N to a
     fixed output file.

  4. BISICLES on the next iteration reads N from that file via
     EffectivePressureFromFile.file (overridden on the command line).

Data exchange files
───────────────────
  bisicles → SUHMO:  <b2s_file>  (symlink to latest BISICLES plot file)
                      contains: thickness, xVel, yVel
  SUHMO → bisicles:  <s2b_file>  (written directly by SUHMO)
                      contains: effectivePressure

Usage
─────
    python run_coupled.py [options]

    # Quick test (serial, 5 coupling iters of 5 BISICLES steps each):
    python run_coupled.py --nprocs 1 --coupling-steps 5 --num-coupling-iters 5

    # Parallel run with custom executables:
    python run_coupled.py \\
        --nprocs 8 \\
        --bisicles-exe ./driver2d.Linux.64.mpicxx.gfortran.OPT.MPI.ex \\
        --suhmo-exe    ./Suhmo2d.Linux.64.mpicxx.gfortran.OPT.MPI.ex  \\
        --bisicles-input inputs.bisicles \\
        --suhmo-input    inputs.suhmo   \\
        --coupling-steps 10 --num-coupling-iters 20
"""
import subprocess
import os
import glob
import argparse
import sys
import shutil


# ── helpers ──────────────────────────────────────────────────────────────────

def find_latest_file(prefix, dim="2d"):
    """Return the lexicographically last file matching <prefix>*.2d.hdf5."""
    pattern = f"{prefix}*.{dim}.hdf5"
    files = sorted(glob.glob(pattern))
    return files[-1] if files else None


def symlink_force(src, dst):
    """Create (or replace) a symlink dst → src."""
    if os.path.islink(dst) or os.path.exists(dst):
        os.remove(dst)
    os.symlink(src, dst)
    print(f"[COUPLED] linked {src} -> {dst}")


def run_cmd(cmd, label):
    """Print and execute a shell command; abort on non-zero exit."""
    print(f"\n{'='*60}")
    print(f"[COUPLED] {label}")
    print(f"  cmd: {cmd}")
    print(f"{'='*60}\n")
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        print(f"[COUPLED] ERROR: '{label}' exited with code {result.returncode}")
        sys.exit(result.returncode)


# ── model runners ─────────────────────────────────────────────────────────────

def write_bisicles_input(args, coupling_iter, restart_chk):
    """
    Write a per-iteration BISICLES input file by copying the base input file
    and appending the coupling overrides at the end.  Returns the path to the
    temporary file (caller is responsible for deleting it after the run).
    """
    tmp_input = f"inputs.bisicles.iter{coupling_iter:04d}"
    shutil.copy2(args.bisicles_input, tmp_input)

    with open(tmp_input, "a") as f:
        f.write(f"\n# --- coupling iteration {coupling_iter} overrides ---\n")

        # Advance up to the end of this coupling window
        max_step = args.coupling_steps * (coupling_iter + 1)
        f.write(f"main.maxStep = {max_step}\n")

        # Write a checkpoint at the end of each window
        f.write(f"amr.check_interval = {args.coupling_steps}\n")

        # Restart from previous checkpoint (not on the very first iteration)
        if restart_chk is not None:
            f.write(f"amr.restart_file = {restart_chk}\n")

        # Supply effective pressure from SUHMO once available (iteration > 0).
        # On the first iteration the input file defaults to hydrostatic N.
        if os.path.exists(args.s2b_file):
            f.write(f"main.effectivePressure = LevelData\n")
            f.write(f"LevelDataEffectivePressure.file = {args.s2b_file}\n")

        # Per-iteration pout name so logs are not overwritten between restarts.
        f.write(f"main.poutBaseName = pout.bisicles.iter{coupling_iter:04d}\n")

    return tmp_input


def run_bisicles(args, coupling_iter, restart_chk=None):
    """
    Run BISICLES for one coupling window.

    On the first call (coupling_iter == 0, restart_chk is None) BISICLES
    starts from scratch.  On subsequent calls it restarts from the checkpoint
    written at the end of the previous window.

    Coupling overrides (maxStep, restart_file, effectivePressure) are written
    into a temporary copy of the input file so that ParmParse picks them up
    reliably.
    """
    tmp_input = write_bisicles_input(args, coupling_iter, restart_chk)

    cmd = f"mpirun -n {args.nprocs} {args.bisicles_exe} {tmp_input}"
    print(f"[COUPLED] Running BISICLES with input file {tmp_input} "
            f"for iteration {coupling_iter}")
    run_cmd(cmd, f"BISICLES iteration {coupling_iter}")

    # Symlink the latest BISICLES plot file to the agreed coupling filename
    latest_plot = find_latest_file(args.bisicles_plot_prefix)
    if latest_plot is None:
        print(f"[COUPLED] ERROR: no BISICLES plot file found "
              f"(prefix='{args.bisicles_plot_prefix}')")
        sys.exit(1)
    symlink_force(latest_plot, args.b2s_file)


def run_suhmo(args, coupling_iter):
    """
    Run SUHMO to steady state for the current ice geometry.

    SUHMO reads ice thickness and velocity from the BISICLES coupling file
    and writes the effective pressure field to a fixed output file.

    Coupling overrides are written into a temporary copy of the SUHMO input
    file so that ParmParse picks them up reliably:
      suhmo.coupled_to_bisicles = true
      suhmo.bisicles_input_file = <b2s_file>
      suhmo.output_N_file       = <s2b_file>
    """
    tmp_input = f"/home/tm17544/BISICLES/SUHMO-BISICLES-sep-exec/exec/BisiclesCoupled/inputs.suhmo.iter{coupling_iter:04d}"
    shutil.copy2(args.suhmo_input, tmp_input)

    with open(tmp_input, "a") as f:
        f.write(f"\n# --- coupling iteration {coupling_iter} overrides ---\n")
        f.write(f"suhmo.coupled_to_bisicles = true\n")
        f.write(f"suhmo.bisicles_input_file = {args.b2s_file}\n")
        f.write(f"suhmo.output_N_file = {args.s2b_file}\n")
        # Per-iteration prefixes so pout/plot/chk files are not overwritten.
        f.write(f"main.poutBaseName = pout.suhmo.iter{coupling_iter:04d}\n")
        f.write(f"AmrHydro.plot_prefix = plot.suhmo.iter{coupling_iter:04d}.\n")
        f.write(f"AmrHydro.check_prefix = chk.suhmo.iter{coupling_iter:04d}.\n")

    cmd = f"mpirun -n {args.nprocs} {args.suhmo_exe} {tmp_input}"
    run_cmd(cmd, f"SUHMO steady-state solve, iteration {coupling_iter}")

    if not os.path.exists(args.s2b_file):
        print(f"[COUPLED] ERROR: SUHMO did not produce '{args.s2b_file}'")
        sys.exit(1)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Orchestrate iteratively coupled BISICLES–SUHMO simulation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--nprocs", type=int, default=4,
                        help="MPI process count for both models")

    parser.add_argument("--coupling-steps", type=int, default=5,
                        help="BISICLES timesteps per coupling interval")
    parser.add_argument("--num-coupling-iters", type=int, default=5,
                        help="Total number of coupling iterations")

    # Executables
    parser.add_argument("--bisicles-exe", type=str,
                        default="/home/tm17544/BISICLES/bisicles-uob-SUHMO-sep-exec/code/exec2D/driver2d.Linux.64.mpiCC.gfortran.DEBUG.OPT.MPI.ex",
                        help="Path to BISICLES driver executable")
    parser.add_argument("--suhmo-exe", type=str,
                        default="/home/tm17544/BISICLES/SUHMO-BISICLES-sep-exec/exec/BisiclesCoupled/Suhmo2d.Linux.64.mpiCC.gfortran.DEBUG.OPT.MPI.ex",
                        help="Path to SUHMO executable")

    # Input files
    parser.add_argument("--bisicles-input", type=str,
                        default="inputs.bisicles",
                        help="BISICLES ParmParse input file")
    parser.add_argument("--suhmo-input", type=str,
                        default="/home/tm17544/BISICLES/SUHMO-BISICLES-sep-exec/exec/BisiclesCoupled/inputs.suhmo",
                        help="SUHMO ParmParse input file")

    # Plot prefix used to locate BISICLES output (must match amr.plot_prefix)
    parser.add_argument("--bisicles-plot-prefix", type=str,
                        default="plot.bisicles.",
                        help="BISICLES amr.plot_prefix (used to find latest plot file)")

    # Coupling exchange filenames
    parser.add_argument("--b2s-file", type=str,
                        default="bisicles_for_suhmo.2d.hdf5",
                        help="BISICLES→SUHMO coupling file (symlink to latest plot)")
    parser.add_argument("--s2b-file", type=str,
                        default="suhmo_for_bisicles.2d.hdf5",
                        help="SUHMO→BISICLES effective pressure file")

    args = parser.parse_args()

    bisicles_chk = None

    for iteration in range(args.num_coupling_iters):
        print(f"\n{'#'*60}")
        print(f"# COUPLING ITERATION {iteration} / {args.num_coupling_iters - 1}")
        print(f"{'#'*60}")

        # ── BISICLES ─────────────────────────────────────────────────────
        run_bisicles(args, iteration, restart_chk=bisicles_chk)

        # Locate checkpoint for next restart
        bisicles_chk = find_latest_file("chk.bisicles.")
        if bisicles_chk is None and iteration < args.num_coupling_iters - 1:
            print("[COUPLED] WARNING: no BISICLES checkpoint found for restart")

        # ── SUHMO ────────────────────────────────────────────────────────
        run_suhmo(args, iteration)

    print(f"\n{'#'*60}")
    print(f"# COUPLED SIMULATION COMPLETE "
          f"({args.num_coupling_iters} iterations, "
          f"{args.num_coupling_iters * args.coupling_steps} BISICLES steps)")
    print(f"{'#'*60}")


if __name__ == "__main__":
    main()
