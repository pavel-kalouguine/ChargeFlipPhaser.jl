# ChargeFlipPhaser.jl

**ChargeFlipPhaser.jl** is a pure Julia implementation of a **generalized charge-flipping phasing algorithm** for resolving crystal structures with **arbitrary space groups in any dimension**. 

Designed as a flexible and customizable alternative to [Superflip](http://superflip.fzu.cz/), it features a modern architecture that prioritizes extensibility and user control.

---

## Features

- **Arbitrary space groups and dimensionality**: Leverages [SpaceGroups.jl](https://github.com/pkfrance/SpaceGroups.jl) for universal symmetry support.
- **Generalized charge-flipping algorithm**: Robust implementation for versatile phasing applications.
- **Symmetry-aware sampling**: Explicitly breaks structural symmetry to eliminate redundancy (unlike Superflip's P1 workaround).
- **Modular plugin system**:
  - Customizable phasing process control.
  - Real-time visualization of intermediate results.
  - Support for user-defined algorithm extensions.
- **Pure Julia codebase**: Ensures transparency and cross-platform compatibility.
- **Research-focused**: Ideal for advanced use cases requiring full algorithmic flexibility.

---

## Comparison with Superflip

| Feature               | ChargeFlipPhaser.jl                          | Superflip                          |
|-----------------------|---------------------------------------------|-----------------------------------|
| Language              | Pure Julia                                  | Fortran 90                        |
| Space Group Support   | Arbitrary (via SpaceGroups.jl)              | Built-in list or manual input     |
| Redundancy Handling   | Symmetry-breaking sampling                 | P1 workaround                    |
| Extensibility         | Plugin system, user extensions              | Not extensible                   |
| Visualization         | Customizable, interactive                  | Static output                    |

---

## Installation

To install the unregistered package, run the following in the Julia `Pkg` REPL (accessed by pressing `]` in the Julia REPL):

```julia
add https://github.com/pavel-kalouguine/ChargeFlipPhaser.jl
```

## Documentation
*Work in progress*. For now, please refer to the source code and docstrings.

## Contributing
Contributions, bug reports, and feature suggestions are welcome! Feel free to open an issue or a pull request.