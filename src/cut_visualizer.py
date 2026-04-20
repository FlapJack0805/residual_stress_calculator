#!/usr/bin/env python3
"""
Simple cut-by-cut visualizer for the residual stress pipeline.

Features
- Loads mesh_BASE.vtk and step_###/mesh_*.vtk from a case directory
- Optional overlay of toolpath STL files
- Keyboard controls:
    n / Right  -> next step
    p / Left   -> previous step
    t          -> toggle toolpath overlay
    w          -> toggle wireframe for active mesh
    o          -> reset camera
    q          -> quit
- Optional slider widget for jumping to a specific step

Recommended install:
    python3 -m pip install pyvista pyvistaqt vtk

Usage:
    python3 cut_visualizer.py --case /path/to/case_000
    python3 cut_visualizer.py --case /path/to/case_000 --toolpaths /path/to/toolpath_stls
"""

from __future__ import annotations

import argparse
import glob
import os
import re
from dataclasses import dataclass
from typing import List, Optional

import pyvista as pv


@dataclass
class StepData:
    label: str
    mesh_path: str
    toolpath_path: Optional[str] = None


STEP_DIR_RE = re.compile(r"step_(\d+)$")
CUT_MESH_RE = re.compile(r"mesh_.*\.vtk$", re.IGNORECASE)
TOOLPATH_RE = re.compile(r"(cut|toolpath).*?(\d+).*\.stl$", re.IGNORECASE)


def natural_key(path: str):
    parts = re.split(r"(\d+)", os.path.basename(path))
    out = []
    for part in parts:
        out.append(int(part) if part.isdigit() else part.lower())
    return out


def find_base_mesh(case_dir: str) -> str:
    candidates = [
        os.path.join(case_dir, "mesh_BASE.vtk"),
        *glob.glob(os.path.join(case_dir, "mesh_*.vtk")),
    ]
    for path in sorted(candidates, key=natural_key):
        if os.path.isfile(path):
            return path
    raise FileNotFoundError(f"Could not find a base mesh in {case_dir}")


def find_step_meshes(case_dir: str) -> List[StepData]:
    step_entries: List[StepData] = []
    for entry in sorted(os.listdir(case_dir)):
        full = os.path.join(case_dir, entry)
        if not os.path.isdir(full):
            continue
        m = STEP_DIR_RE.match(entry)
        if not m:
            continue

        vtk_candidates = sorted(
            glob.glob(os.path.join(full, "mesh_*.vtk")), key=natural_key
        )
        if not vtk_candidates:
            continue

        step_idx = int(m.group(1))
        label = f"CUT {step_idx:03d}"
        step_entries.append(StepData(label=label, mesh_path=vtk_candidates[0]))

    return step_entries


def assign_toolpaths(steps: List[StepData], toolpath_dir: Optional[str]) -> None:
    if not toolpath_dir or not os.path.isdir(toolpath_dir):
        return

    candidates = sorted(glob.glob(os.path.join(toolpath_dir, "*.stl")), key=natural_key)
    step_to_tool = {}
    for path in candidates:
        name = os.path.basename(path)
        m = TOOLPATH_RE.search(name)
        if m:
            step_to_tool[int(m.group(2))] = path

    # Fallback: if counts match, assign in order.
    unassigned = [i for i, s in enumerate(steps, start=1) if i not in step_to_tool]
    if unassigned and len(candidates) == len(steps):
        for idx, path in enumerate(candidates, start=1):
            step_to_tool.setdefault(idx, path)

    for idx, step in enumerate(steps, start=1):
        step.toolpath_path = step_to_tool.get(idx)


class CutViewer:
    def __init__(self, steps: List[StepData], show_edges: bool = False) -> None:
        self.steps = steps
        self.index = 0
        self.show_toolpath = True
        self.show_edges = show_edges

        self.plotter = pv.Plotter(window_size=(1400, 900))
        self.plotter.set_background("white")

        self.mesh_actor = None
        self.tool_actor = None
        self.text_actor = None
        self.bounds_actor = None

        self._setup_ui()
        self._render_current(reset_camera=True)

    def _setup_ui(self) -> None:
        self.plotter.add_axes()
        self.plotter.add_key_event("n", self.next_step)
        self.plotter.add_key_event("Right", self.next_step)
        self.plotter.add_key_event("p", self.prev_step)
        self.plotter.add_key_event("Left", self.prev_step)
        self.plotter.add_key_event("t", self.toggle_toolpath)
        self.plotter.add_key_event("w", self.toggle_wireframe)
        self.plotter.add_key_event("o", self.reset_camera)
        self.plotter.add_key_event("q", self.plotter.close)

        if len(self.steps) > 1:
            self.plotter.add_slider_widget(
                callback=self._slider_callback,
                rng=[0, len(self.steps) - 1],
                value=0,
                title="Step",
                pointa=(0.025, 0.08),
                pointb=(0.32, 0.08),
                style="modern",
                event_type="always",
            )

        help_text = (
            "n/Right: next | p/Left: previous | t: toolpath | "
            "w: wireframe | o: reset view | q: quit"
        )
        self.plotter.add_text(help_text, position="lower_left", font_size=10, color="black")

    def _slider_callback(self, value: float) -> None:
        new_index = int(round(value))
        if new_index != self.index:
            self.index = max(0, min(len(self.steps) - 1, new_index))
            self._render_current(reset_camera=False)

    def next_step(self) -> None:
        if self.index < len(self.steps) - 1:
            self.index += 1
            self._render_current(reset_camera=False)

    def prev_step(self) -> None:
        if self.index > 0:
            self.index -= 1
            self._render_current(reset_camera=False)

    def toggle_toolpath(self) -> None:
        self.show_toolpath = not self.show_toolpath
        self._render_current(reset_camera=False)

    def toggle_wireframe(self) -> None:
        self.show_edges = not self.show_edges
        self._render_current(reset_camera=False)

    def reset_camera(self) -> None:
        self.plotter.reset_camera()
        self.plotter.render()

    def _clear_actors(self) -> None:
        for actor_name in ["mesh_actor", "tool_actor", "text_actor", "bounds_actor"]:
            actor = getattr(self, actor_name)
            if actor is not None:
                try:
                    self.plotter.remove_actor(actor)
                except Exception:
                    pass
                setattr(self, actor_name, None)

    def _load_mesh(self, path: str):
        mesh = pv.read(path)
        if mesh.n_points == 0:
            raise RuntimeError(f"Empty mesh: {path}")
        return mesh

    def _render_current(self, reset_camera: bool) -> None:
        self._clear_actors()
        step = self.steps[self.index]

        mesh = self._load_mesh(step.mesh_path)
        scalars = "elem_label" if "elem_label" in mesh.array_names else None

        self.mesh_actor = self.plotter.add_mesh(
            mesh,
            scalars=scalars,
            show_edges=self.show_edges,
            opacity=0.95,
            cmap="viridis",
            scalar_bar_args={"title": "elem_label"} if scalars else None,
            color=None if scalars else "lightgray",
        )

        if self.show_toolpath and step.toolpath_path and os.path.isfile(step.toolpath_path):
            tool = self._load_mesh(step.toolpath_path)
            self.tool_actor = self.plotter.add_mesh(
                tool,
                style="wireframe",
                line_width=2,
                color="red",
                opacity=0.8,
            )

        bounds = mesh.bounds
        info = [
            f"Step {self.index}/{len(self.steps) - 1}: {step.label}",
            f"Mesh: {os.path.basename(step.mesh_path)}",
            f"Points: {mesh.n_points}  Cells: {mesh.n_cells}",
            (
                f"Bounds: X[{bounds[0]:.3f}, {bounds[1]:.3f}]  "
                f"Y[{bounds[2]:.3f}, {bounds[3]:.3f}]  "
                f"Z[{bounds[4]:.3f}, {bounds[5]:.3f}]"
            ),
            (
                f"Toolpath: {os.path.basename(step.toolpath_path)}"
                if step.toolpath_path else "Toolpath: none"
            ),
        ]
        self.text_actor = self.plotter.add_text(
            "\n".join(info), position="upper_left", font_size=11, color="black"
        )

        self.bounds_actor = self.plotter.show_bounds(
            grid="back",
            location="outer",
            all_edges=True,
            ticks="outside",
            xtitle="X",
            ytitle="Y",
            ztitle="Z",
            font_size=10,
            use_2d=False,
        )

        if reset_camera:
            self.plotter.reset_camera()
            self.plotter.camera_position = "iso"
        self.plotter.render()

    def show(self) -> None:
        self.plotter.show()


def build_steps(case_dir: str, toolpath_dir: Optional[str]) -> List[StepData]:
    base_mesh = find_base_mesh(case_dir)
    steps = [StepData(label="BASE", mesh_path=base_mesh)]
    steps.extend(find_step_meshes(case_dir))
    assign_toolpaths(steps[1:], toolpath_dir)
    return steps


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Visualize per-cut VTK meshes and optional toolpath STLs.")
    ap.add_argument("--case", required=True, help="Path to case directory containing mesh_BASE.vtk and step_### folders")
    ap.add_argument("--toolpaths", default=None, help="Optional directory containing toolpath STL files")
    ap.add_argument("--edges", action="store_true", help="Show mesh edges initially")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    steps = build_steps(args.case, args.toolpaths)
    viewer = CutViewer(steps, show_edges=args.edges)
    viewer.show()


if __name__ == "__main__":
    main()
