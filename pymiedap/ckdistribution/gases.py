"""
Registry of HITRAN molecules supported by the ckdistribution package.

Each entry contains:
    hitran_id      : HITRAN molecule number (M field in .par files)
    molar_mass     : molar mass in g/mol (main isotopologue)
    global_iso_ids : list of HITRAN *global* isotopologue IDs for fetch_by_ids()
    local_iso_ids  : list of HITRAN *local* isotopologue IDs (1-based, per molecule)
    name           : human-readable name

Global isotopologue IDs follow the 2004+ HITRAN numbering scheme used by HAPI's
``fetch_by_ids``.  Local IDs follow the 1-based per-molecule ordering used by
HAPI's ``fetch``.

Reference: https://hitran.org/docs/iso-meta/
"""

# ---------------------------------------------------------------------------
# HITRAN molecule registry
# ---------------------------------------------------------------------------
GASES = {
    # Molecule  hitran_id  molar_mass(g/mol)  global_iso_ids         local_iso_ids  name
    'H2O': {
        'hitran_id':    1,
        'molar_mass':  18.01528,
        'global_iso_ids': [1, 2, 3, 4, 5, 6],
        'local_iso_ids':  [1, 2, 3, 4, 5, 6],
        'name': 'Water vapour',
    },
    'CO2': {
        'hitran_id':    2,
        'molar_mass':  44.01000,
        'global_iso_ids': [7, 8, 9, 10, 11, 12, 13, 14, 121, 15, 120, 122],
        'local_iso_ids':  [1, 2, 3, 4, 5,  6,  7,  8,   9, 10,  11,  12],
        'name': 'Carbon dioxide',
    },
    'O3': {
        'hitran_id':    3,
        'molar_mass':  47.99820,
        'global_iso_ids': [16, 17, 18],
        'local_iso_ids':  [1,  2,  3],
        'name': 'Ozone',
    },
    'N2O': {
        'hitran_id':    4,
        'molar_mass':  44.01280,
        'global_iso_ids': [19, 20, 21, 22, 23],
        'local_iso_ids':  [1,  2,  3,  4,  5],
        'name': 'Nitrous oxide',
    },
    'CO': {
        'hitran_id':    5,
        'molar_mass':  27.99491,
        'global_iso_ids': [24, 25, 26, 27, 28, 29],
        'local_iso_ids':  [1,  2,  3,  4,  5,  6],
        'name': 'Carbon monoxide',
    },
    'CH4': {
        'hitran_id':    6,
        'molar_mass':  16.04303,
        'global_iso_ids': [30, 31, 32, 33],
        'local_iso_ids':  [1,  2,  3,  4],
        'name': 'Methane',
    },
    'O2': {
        'hitran_id':    7,
        'molar_mass':  31.99880,
        'global_iso_ids': [34, 35, 36],
        'local_iso_ids':  [1,  2,  3],
        'name': 'Oxygen',
    },
    'NO': {
        'hitran_id':    8,
        'molar_mass':  29.99799,
        'global_iso_ids': [37, 38, 39],
        'local_iso_ids':  [1,  2,  3],
        'name': 'Nitric oxide',
    },
    'SO2': {
        'hitran_id':    9,
        'molar_mass':  63.96200,
        'global_iso_ids': [40, 41],
        'local_iso_ids':  [1,  2],
        'name': 'Sulphur dioxide',
    },
    'NO2': {
        'hitran_id':   10,
        'molar_mass':  45.99291,
        'global_iso_ids': [42],
        'local_iso_ids':  [1],
        'name': 'Nitrogen dioxide',
    },
    'NH3': {
        'hitran_id':   11,
        'molar_mass':  17.02655,
        'global_iso_ids': [43, 44],
        'local_iso_ids':  [1,  2],
        'name': 'Ammonia',
    },
    'HNO3': {
        'hitran_id':   12,
        'molar_mass':  62.99564,
        'global_iso_ids': [45],
        'local_iso_ids':  [1],
        'name': 'Nitric acid',
    },
    'OH': {
        'hitran_id':   13,
        'molar_mass':  17.00274,
        'global_iso_ids': [46, 47, 48],
        'local_iso_ids':  [1,  2,  3],
        'name': 'Hydroxyl radical',
    },
    'HF': {
        'hitran_id':   14,
        'molar_mass':  20.00623,
        'global_iso_ids': [49, 110],
        'local_iso_ids':  [1,   2],
        'name': 'Hydrogen fluoride',
    },
    'HCl': {
        'hitran_id':   15,
        'molar_mass':  35.97668,
        'global_iso_ids': [50, 51, 107, 108],
        'local_iso_ids':  [1,  2,   3,   4],
        'name': 'Hydrogen chloride',
    },
    'HBr': {
        'hitran_id':   16,
        'molar_mass':  79.90439,
        'global_iso_ids': [52, 53, 111, 112],
        'local_iso_ids':  [1,  2,   3,   4],
        'name': 'Hydrogen bromide',
    },
    'HI': {
        'hitran_id':   17,
        'molar_mass': 127.91239,
        'global_iso_ids': [54],
        'local_iso_ids':  [1],
        'name': 'Hydrogen iodide',
    },
    'ClO': {
        'hitran_id':   18,
        'molar_mass':  50.96371,
        'global_iso_ids': [55, 56],
        'local_iso_ids':  [1,  2],
        'name': 'Chlorine monoxide',
    },
    'OCS': {
        'hitran_id':   19,
        'molar_mass':  59.96685,
        'global_iso_ids': [57, 58, 59, 60, 61],
        'local_iso_ids':  [1,  2,  3,  4,  5],
        'name': 'Carbonyl sulphide',
    },
    'H2CO': {
        'hitran_id':   20,
        'molar_mass':  30.01057,
        'global_iso_ids': [62, 63, 64],
        'local_iso_ids':  [1,  2,  3],
        'name': 'Formaldehyde',
    },
    'HOCl': {
        'hitran_id':   21,
        'molar_mass':  51.97181,
        'global_iso_ids': [65, 66],
        'local_iso_ids':  [1,  2],
        'name': 'Hypochlorous acid',
    },
    'N2': {
        'hitran_id':   22,
        'molar_mass':  28.01348,
        'global_iso_ids': [67],
        'local_iso_ids':  [1],
        'name': 'Molecular nitrogen',
    },
    'HCN': {
        'hitran_id':   23,
        'molar_mass':  27.01090,
        'global_iso_ids': [68, 69, 70],
        'local_iso_ids':  [1,  2,  3],
        'name': 'Hydrogen cyanide',
    },
    'CH3Cl': {
        'hitran_id':   24,
        'molar_mass':  49.99232,
        'global_iso_ids': [71, 72],
        'local_iso_ids':  [1,  2],
        'name': 'Methyl chloride',
    },
    'H2O2': {
        'hitran_id':   25,
        'molar_mass':  34.00548,
        'global_iso_ids': [73],
        'local_iso_ids':  [1],
        'name': 'Hydrogen peroxide',
    },
    'C2H2': {
        'hitran_id':   26,
        'molar_mass':  26.01565,
        'global_iso_ids': [74, 75, 105],
        'local_iso_ids':  [1,  2,   3],
        'name': 'Acetylene',
    },
    'C2H6': {
        'hitran_id':   27,
        'molar_mass':  30.04695,
        'global_iso_ids': [76, 77],
        'local_iso_ids':  [1,  2],
        'name': 'Ethane',
    },
    'PH3': {
        'hitran_id':   28,
        'molar_mass':  33.99738,
        'global_iso_ids': [78],
        'local_iso_ids':  [1],
        'name': 'Phosphine',
    },
    'COF2': {
        'hitran_id':   29,
        'molar_mass':  65.99156,
        'global_iso_ids': [79],
        'local_iso_ids':  [1],
        'name': 'Carbonyl fluoride',
    },
    'SF6': {
        'hitran_id':   30,
        'molar_mass': 145.95362,
        'global_iso_ids': [80],
        'local_iso_ids':  [1],
        'name': 'Sulphur hexafluoride',
    },
    'H2S': {
        'hitran_id':   31,
        'molar_mass':  33.98772,
        'global_iso_ids': [81, 82, 83],
        'local_iso_ids':  [1,  2,  3],
        'name': 'Hydrogen sulphide',
    },
    'HCOOH': {
        'hitran_id':   32,
        'molar_mass':  46.00548,
        'global_iso_ids': [84],
        'local_iso_ids':  [1],
        'name': 'Formic acid',
    },
    'HO2': {
        'hitran_id':   33,
        'molar_mass':  32.99746,
        'global_iso_ids': [85],
        'local_iso_ids':  [1],
        'name': 'Hydroperoxyl radical',
    },
    'C2H4': {
        'hitran_id':   38,
        'molar_mass':  28.03130,
        'global_iso_ids': [92, 93],
        'local_iso_ids':  [1,  2],
        'name': 'Ethylene',
    },
    'CH3OH': {
        'hitran_id':   39,
        'molar_mass':  32.04187,
        'global_iso_ids': [94],
        'local_iso_ids':  [1],
        'name': 'Methanol',
    },
}

# ---------------------------------------------------------------------------
# Convenience look-up helpers
# ---------------------------------------------------------------------------

def get_gas(name: str) -> dict:
    """Return the registry entry for *name* (case-insensitive).

    Parameters
    ----------
    name : str
        Molecule name, e.g. ``'CO2'``, ``'ch4'``.

    Returns
    -------
    dict
        Registry entry with keys ``hitran_id``, ``molar_mass``,
        ``global_iso_ids``, ``local_iso_ids``, ``name``.

    Raises
    ------
    KeyError
        If the molecule is not in the registry.
    """
    key = name.upper()
    if key not in GASES:
        # Try a case-insensitive scan
        for k in GASES:
            if k.upper() == key:
                return GASES[k]
        raise KeyError(
            f"Molecule '{name}' not found in gases registry. "
            f"Available: {list(GASES.keys())}"
        )
    return GASES[key]


def molar_mass(name: str) -> float:
    """Return molar mass (g/mol) for *name*."""
    return get_gas(name)['molar_mass']


def hitran_id(name: str) -> int:
    """Return the HITRAN molecule number for *name*."""
    return get_gas(name)['hitran_id']


def list_gases() -> list:
    """Return a sorted list of all registered molecule names."""
    return sorted(GASES.keys())
