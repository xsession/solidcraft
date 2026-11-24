# Making this site

The pages are Markdown in `docs/multiphysics-site/docs`, built with
[MkDocs](https://www.mkdocs.org/) and the
[Material for MkDocs](https://squidfunk.github.io/mkdocs-material/) theme. Both are build
tools only (BSD-2-Clause and MIT); nothing from them is part of SolidCraft.

## Build the site

```sh
python3 -m venv .venv
.venv/bin/pip install mkdocs mkdocs-material
.venv/bin/mkdocs build --strict -f docs/multiphysics-site/mkdocs.yml
```

The site is written to `target/multiphysics-site`. `mkdocs serve -f
docs/multiphysics-site/mkdocs.yml` serves it at <http://127.0.0.1:8000> while you edit. The
site fetches no web fonts or scripts from elsewhere, so it reads offline.

## The example design

`docs/examples/board-in-duct.sc` holds the tutorials' design and two real runs. To make it
again, start SolidCraft with the solver companion and a control port, then run the script:

```sh
SOLIDCRAFT_CFD_SOLVE=/path/to/cfd-solve solidcraft --control 7979 &
cd docs/multiphysics-site/tools
SOLIDCRAFT_CONTROL_TOKEN=$(cat ~/.local/share/solidcraft/control/7979.token) \
    python3 make_study.py ../docs/examples/board-in-duct.sc
```

The script sends commands only (no clicks), so it runs the same on any desktop.

## The pictures

Every picture is drawn by SolidCraft from the example design, on the GPU, by a test of the
Multiphysics workspace that opens the design, sets each scene up through the workspace's own
panels and saves the frame:

```sh
cargo test -p solidcraft-multiphysics-ui --lib documentation_pictures -- --ignored
```

It writes `docs/images/*.png` and the batch report in `docs/examples/report`. Without a GPU
adapter it skips. Run it again after a change to the workspace, and the pictures follow.
