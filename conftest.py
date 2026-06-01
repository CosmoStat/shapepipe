"""Shared pytest configuration for the ShapePipe test suite."""

import os

from hypothesis import settings


settings.register_profile("ci", derandomize=True, max_examples=50)
settings.register_profile("dev", max_examples=200)
settings.load_profile(os.environ.get("HYPOTHESIS_PROFILE", "ci"))
