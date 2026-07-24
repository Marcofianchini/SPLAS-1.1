PYTHON ?= python3
FC = gfortran
export FC
F90FLAGS ?= -O2 -fPIC -ffree-form -ffree-line-length-none
DEBUG_F90FLAGS ?= -O0 -g -fPIC -ffree-form -ffree-line-length-none -Wall -Wextra -Wimplicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
PACKAGE_DIR := src/splas
CORE_SOURCE := src/fortran/splas_core.f90

.PHONY: all build build-debug install test test-python smoke reproduce-thesis clean audit

all: build

build:
	@command -v $(FC) >/dev/null || { echo "ERROR: Fortran compiler '$(FC)' not found."; exit 2; }
	$(PYTHON) -m numpy.f2py -c -m _core $(CORE_SOURCE) --backend meson --f90flags="$(F90FLAGS)" --build-dir build/f2py
	$(PYTHON) scripts/install_extension.py --module "_core*.so" --destination "$(PACKAGE_DIR)"
	PYTHONPATH=src $(PYTHON) -c "from splas.backend import require_core; require_core(); print('SPLAS Fortran core import: OK')"

build-debug:
	@command -v $(FC) >/dev/null || { echo "ERROR: Fortran compiler '$(FC)' not found."; exit 2; }
	$(PYTHON) -m numpy.f2py -c -m _core $(CORE_SOURCE) --backend meson --f90flags="$(DEBUG_F90FLAGS)" --build-dir build/f2py-debug
	$(PYTHON) scripts/install_extension.py --module "_core*.so" --destination "$(PACKAGE_DIR)"
	PYTHONPATH=src $(PYTHON) -c "from splas.backend import require_core; require_core(); print('SPLAS debug core import: OK')"

install:
	$(PYTHON) -m pip install -e .

test: build
	PYTHONPATH=src MPLCONFIGDIR=/tmp/splas-mpl SPLAS_REQUIRE_CORE=1 $(PYTHON) -m unittest discover -s tests -v

test-python:
	PYTHONPATH=src MPLCONFIGDIR=/tmp/splas-mpl $(PYTHON) -m unittest discover -s tests -v

smoke: build
	PYTHONPATH=src $(PYTHON) -m splas run configs/reference/smoke.json --output reproduction/results/smoke
	PYTHONPATH=src $(PYTHON) -m splas analyze reproduction/results/smoke --output reproduction/tables/smoke
	PYTHONPATH=src MPLCONFIGDIR=/tmp/splas-mpl $(PYTHON) -m splas figures reproduction/results/smoke --output reproduction/figures/smoke

reproduce-thesis: build
	PYTHONPATH=src MPLCONFIGDIR=/tmp/splas-mpl $(PYTHON) -m splas reproduce-thesis --output reproduction

audit:
	PYTHONPATH=src $(PYTHON) scripts/baseline_audit.py

clean:
	$(PYTHON) scripts/clean_build.py
