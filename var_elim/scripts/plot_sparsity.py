#  ___________________________________________________________________________
#
#  Variable Elimination: Research code for variable elimination in NLPs
#
#  Copyright (c) 2023. Triad National Security, LLC. All rights reserved.
#
#  This program was produced under U.S. Government contract 89233218CNA000001
#  for Los Alamos National Laboratory (LANL), which is operated by Triad
#  National Security, LLC for the U.S. Department of Energy/National Nuclear
#  Security Administration. All rights in the program are reserved by Triad
#  National Security, LLC, and the U.S. Department of Energy/National Nuclear
#  Security Administration. The Government is granted for itself and others
#  acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license
#  in this material to reproduce, prepare derivative works, distribute copies
#  to the public, perform publicly and display publicly, and to permit others
#  to do so.
#
#  This software is distributed under the 3-clause BSD license.
#  ___________________________________________________________________________

import os
import matplotlib.pyplot as plt
from pyomo.contrib.incidence_analysis.visualize import spy_dulmage_mendelsohn

import var_elim.scripts.config as config
from var_elim.elimination_callbacks import matching_elim_callback


plt.rcParams["font.size"] = 14
plt.rcParams["font.family"] = "serif"


def plot_sparsity(m_before, m_after, title=None):
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    spy_dulmage_mendelsohn(
        m_before,
        ax=axes[0],
        highlight_coarse=False,
        highlight_fine=False,
    )
    axes[0].set_title("Original model")
    axes[0].set_xlabel("Variables")
    axes[0].set_ylabel("Constraints")

    spy_dulmage_mendelsohn(
        m_after,
        ax=axes[1],
        highlight_coarse=False,
        highlight_fine=False,
    )
    axes[1].set_title("Linear matching")
    axes[1].set_xlabel("Variables")
    axes[1].set_ylabel("Constraints")

    if title is not None:
        fig.suptitle(title)

    fig.tight_layout()
    return fig, axes


def main(args):
    model = config.CONSTRUCTOR_LOOKUP["mb-steady"]()
    # Copy model for before/after without rebuilding
    before_model = model.clone()
    after_model = model
    matching_elim_callback(after_model)

    fig, _ = plot_sparsity(before_model, after_model, title=args.title)

    if not args.no_save:
        suff = "" if args.suffix is None else f"-{args.suffix}"
        fname = f"mb-steady-matching-sparsity{suff}.pdf"
        image_dir = config.validate_dir(args.image_dir)
        fpath = os.path.join(image_dir, fname)
        print(f"Saving sparsity plot to {fpath}")
        fig.savefig(fpath, transparent=not args.opaque)

    if args.show:
        plt.show()


if __name__ == "__main__":
    argparser = config.get_plot_argparser()
    args = argparser.parse_args()
    if args.model is None:
        args.model = "mb-steady"
    main(args)
