#!/usr/bin/env bash
# Fetch mathlib, download its prebuilt cache, build the smoke test. Run from repo root.
set -e
export PATH="$HOME/.elan/bin:$PATH"
cd lean
lake update
lake exe cache get
lake build
