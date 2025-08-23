## Neutron Decay Event Generator [![license](https://img.shields.io/badge/License-GPLv3-blue.svg)](COPYRIGHT)
[![CI-CD](https://github.com/jose-luis-rs/R3B_Neutron_Generator/actions/workflows/main.yml/badge.svg)](https://github.com/jose-luis-rs/R3B_Neutron_Generator/actions/workflows/main.yml)

This file is used to generate input data for the [R3BRoot](https://github.com/R3BRootGroup/R3BRoot) simulation framework. It supports the generation of 1n and 2n decay events, including multiple decay mechanisms for the 2n decay mode.

### License
This software is distributed under the terms of the GNU Lesser General Public Licence version 3 ([LGPLv3](https://github.com/jose-luis-rs/R3B_Neutron_Generator/blob/main/LICENSE)).

### Decay Mode Options
The type of decay can be configured using the `decay_opt` parameter, which accepts the following values:

- `decay_opt = 0` — **Basic parameters for 1n decay**
- `decay_opt = 1` — **N–N correlations for 2n** (based on Lednicky’s formalism)
- `decay_opt = 2` — **Sequential decay for a 2n decay**
- `decay_opt = 3` — **Dineutron decay for 2n**

### How to Use

Clone the repository:

```bash
git clone https://github.com/jose-luis-rs/R3B_Neutron_Generator.git
cd R3B_Neutron_Generator
```

Run the generator using ROOT:

```bash
root -l NeutronDecayGenerator.C
or
root -l 'NeutronDecayGenerator.C("fileName", nb_events, A, Z, Ekin, E_BW, W_BW, nb_n, decay_opt)'
```

| Parameter   | Description                              |
| ----------- | ---------------------------------------- |
| `fileName`  | Output ROOT file name (string)           |
| `nb_events` | Number of events to generate (int)       |
| `A`         | Mass number of the fragment              |
| `Z`         | Atomic number of the fragment            |
| `Ekin`      | Kinetic energy in MeV (lab frame)        |
| `E_BW`      | Average energy for the Breit-Wigner      |
| `W_BW`      | Sigma of the Breit-Wigner                |
| `nb_n`      | Number of decay neutrons                 |
| `decay_opt` | Decay mode selection (see above options) |

### Nuclear Mass Calculation

The nuclear masses used in this generator are computed using the formula:

$$
M_\text{nucleus} = Z \cdot m_p + N \cdot m_n - B
$$

Where:

- `m_p` is the proton mass (938.272088 MeV, [PDG 2022](https://pdg.lbl.gov/2025/listings/rpp2025-list-p.pdf))
- `m_n` is the neutron mass (939.565420 MeV, [PDG 2022](https://pdg.lbl.gov/2025/listings/rpp2025-list-n.pdf))
- `B` is the nuclear binding energy from the [AME2020 atomic mass evaluation](https://www-nds.iaea.org/amdc/) 

For nuclei not listed or with missing binding energy data, the generator returns `-1`.

