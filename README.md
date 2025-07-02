## Neutron Decay Event Generator

This file is used to generate input data for the **R3BRoot** simulation framework. It supports the generation of 1n and 2n decay events, including multiple decay mechanisms for the 2n decay mode.

### Decay Mode Options

The type of decay can be configured using the `decay_opt` parameter, which accepts the following values:

- `decay_opt = 0` — **Phase space decay**
- `decay_opt = 1` — **N–N correlations** (based on Lednicky’s formalism)
- `decay_opt = 2` — **Sequential decay**
- `decay_opt = 3` — **Dineutron decay**

### How to Use

Clone the repository:

```bash
git clone https://github.com/jose-luis-rs/R3B_Neutron_Generator.git
cd R3B_Neutron_Generator

Run the generator using ROOT:

root -l
.L NeutronDecayGenerator.cxx++
EventGenerator_Ndecay("fileName", nb_events, A, Z, Ekin, decay_opt)


| Parameter   | Description                              |
| ----------- | ---------------------------------------- |
| `fileName`  | Output ROOT file name (string)           |
| `nb_events` | Number of events to generate (int)       |
| `A`         | Mass number of the fragment              |
| `Z`         | Atomic number of the fragment            |
| `Ekin`      | Kinetic energy in MeV (lab frame)        |
| `decay_opt` | Decay mode selection (see above options) |

