"""Minimal self-describing NetCDF3 I/O based on SciPy."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.io import netcdf_file


VARIABLES: dict[str, tuple[tuple[str, ...], str, str]] = {
    "P": (("time", "depth", "phyto_size"), "uM N", "phytoplankton nitrogen biomass"),
    "Z": (("time", "depth", "zoo_size"), "uM N", "zooplankton nitrogen biomass"),
    "quota_P": (("time", "depth", "phyto_size"), "uM N cell-1", "phytoplankton cellular nitrogen quota"),
    "cell_number_P": (("time", "depth", "phyto_size"), "cell m-3, historical scaling", "phytoplankton cell number derived as P / quota"),
    "nutrient_pool": (("time", "depth", "nutrient"), "uM", "dissolved pools; N is prognostic and P/C are diagnostic Redfield equivalents"),
    "detritus": (("time", "depth"), "uM N", "detrital nitrogen"),
    "surface_irradiance": (("time",), "W m-2 PAR", "surface light forcing before attenuation"),
    "irradiance": (("time", "depth"), "W m-2 PAR", "attenuated in-water irradiance"),
    "nutrient_limitation": (("time", "depth", "phyto_size"), "1", "fractional nutrient limitation"),
    "grazing_limitation": (("time", "depth", "phyto_size"), "1", "grazing loss normalized by potential growth"),
    "light_limitation": (("time", "depth", "phyto_size"), "1", "fraction of potential growth retained under light"),
    "phyto_growth_rate": (("time", "depth", "phyto_size"), "d-1", "Droop growth rate before light limitation"),
    "grazing_loss": (("time", "depth", "phyto_size"), "uM N d-1", "total grazing loss by phytoplankton class"),
    "net_phyto_growth_rate": (("time", "depth", "phyto_size"), "d-1", "light-limited Droop growth minus grazing, before mortality and transport"),
    "total_growth_reduction": (("time", "depth", "phyto_size"), "1", "reduction relative to mu_max from quota, light, and grazing; may exceed one"),
    "phyto_biomass": (("time", "depth"), "uM N", "total phytoplankton biomass"),
    "zoo_biomass": (("time", "depth"), "uM N", "total zooplankton biomass"),
    "shannon_entropy_P": (("time", "depth"), "nat", "Shannon entropy over phytoplankton size classes"),
    "shannon_evenness_P": (("time", "depth"), "1", "Shannon entropy divided by log of present classes"),
    "diversity_valid_P": (("time", "depth"), "1", "one where biomass is sufficient for diversity interpretation"),
    "carbon_per_cell_P": (("phyto_size",), "pg C cell-1", "Montagnes et al. size-to-carbon relationship"),
    "carbon_biomass_P": (("time", "depth", "phyto_size"), "mmol C m-3", "phytoplankton carbon inferred from cell number"),
    "export_phyto": (("time",), "mg C m-2 d-1", "phytoplankton sinking export at lower boundary"),
    "export_detritus": (("time",), "mg C m-2 d-1", "detrital sinking export at lower boundary"),
    "export_total": (("time",), "mg C m-2 d-1", "total lower-boundary carbon export"),
    "mass_balance_error": (("time",), "uM N m", "closed-system nitrogen inventory change from initial"),
    "nonnegative_correction": (("time",), "uM N m", "cumulative mass introduced by roundoff-level nonnegative correction"),
}


def write_dataset(
    path: str | Path,
    coords: dict[str, np.ndarray],
    data: dict[str, np.ndarray],
    config: dict[str, Any],
    metadata: dict[str, Any],
) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with netcdf_file(destination, "w", version=2) as dataset:
        for name, values in coords.items():
            array = np.asarray(values)
            dataset.createDimension(name, int(array.size))
            variable = dataset.createVariable(name, "f8", (name,))
            variable[:] = array.astype(np.float64)
        dataset.variables["time"].units = "days since simulation start"
        dataset.variables["depth"].units = "m; positive downward; layer centers"
        dataset.variables["phyto_size"].units = "um ESD"
        dataset.variables["zoo_size"].units = "um ESD"
        dataset.variables["nutrient"].labels = "N,P_equivalent,C_equivalent"
        for name, values in data.items():
            if name not in VARIABLES:
                raise KeyError(f"No NetCDF schema registered for variable {name}")
            dims, units, long_name = VARIABLES[name]
            array = np.asarray(values)
            expected = tuple(len(coords[dim]) for dim in dims)
            if array.shape != expected:
                raise ValueError(f"{name} has shape {array.shape}; expected {expected}")
            variable = dataset.createVariable(name, "f8", dims)
            variable[:] = array.astype(np.float64)
            variable.units = units
            variable.long_name = long_name
        dataset.title = "SPLAS 1.1 reproducible simulation"
        dataset.Conventions = "CF-1.8-inspired"
        dataset.config_json = json.dumps(config, sort_keys=True)
        dataset.metadata_json = json.dumps(metadata, sort_keys=True)
    return destination


def read_dataset(path: str | Path) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    source = Path(path)
    with netcdf_file(source, "r", mmap=False) as dataset:
        arrays = {name: np.array(var[:], copy=True) for name, var in dataset.variables.items()}
        attrs = {
            "config": json.loads(_decode(dataset.config_json)),
            "metadata": json.loads(_decode(dataset.metadata_json)),
        }
    return arrays, attrs


def _decode(value: Any) -> str:
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)
