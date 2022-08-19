"""Automated testing linting and formatting apparatus."""

import nox
from nox.sessions import Session

nox.options.sessions = "lint", "tests", "cover", "mypy"  # default session
locations = "sotb_wrapper", "test", "noxfile.py"  # Linting locations
locations_mypy = "sotb_wrapper/interface.py"  # mypy locations
pyversions = ["3.8", "3.9", "3.10"]


# Testing
@nox.session(python=pyversions)
def tests(session: Session) -> None:
    """Run tests."""
    pkg_path = f"{session.virtualenv.location}/lib/python{session.python}/site-packages"
    strg = [
        "--cov",
        f"{pkg_path}/sotb_wrapper",
        "-s",
        "--import-mode=importlib",
        "--cov-report=",
    ]
    args = session.posargs + strg
    session.install("pytest", "pytest-cov")
    session.install(".")
    session.run("pytest", *args, env={"COVERAGE_FILE": f".coverage.{session.python}"})
    session.notify("cover")
    
@nox.session
def cover(session: Session) -> None:
    """Coverage analysis."""
    session.install("coverage[toml]")
    session.run("coverage", "combine")
    session.run("coverage", "report", "--show-missing")
    session.run("coverage", "erase")

# Linting
@nox.session(python="3.9")
def lint(session: Session) -> None:
    """Lint code."""
    args = session.posargs or locations
    session.install(
        "flake8",
        "flake8-black",
        "flake8-bugbear",
        "flake8-import-order",
        "flake8-annotations",
        "flake8-docstrings",
        "darglint",
    )
    session.run("flake8", *args)


# Code formatting
@nox.session(python="3.9")
def black(session: Session) -> None:
    """Format code."""
    args = session.posargs or locations
    session.install("black")
    session.run("black", *args)


# Static typing
@nox.session(python="3.9")
def mypy(session: Session) -> None:
    """Run the static type checker."""
    args = session.posargs or locations_mypy
    session.install("mypy")
    session.install("types-pkg_resources")
    session.run("mypy", args)
