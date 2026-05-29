"""Shared pytest configuration."""

from hypothesis import HealthCheck, settings

settings.register_profile(
    "ci",
    derandomize=True,
    max_examples=50,
    suppress_health_check=[HealthCheck.function_scoped_fixture],
)
settings.load_profile("ci")
