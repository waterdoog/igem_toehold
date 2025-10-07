# 🧬 iGEM Toehold Switches Toolkit

A comprehensive toolkit for designing, analyzing, and visualizing toehold riboswitches for synthetic biology applications in iGEM competitions.

## 📋 Overview

This toolkit provides computational tools for:
- **Trigger RNA Design**: Generate optimal trigger RNA sequences that activate toehold switches
- **Toehold Switch Design**: Design complete toehold switch systems with integrated toxin genes
- **RNA Structure Visualization**: Plot RNA secondary structures with annotated functional regions
- **Structure Validation**: Verify toehold switch functionality through computational analysis

## 🔬 Biological Background

**Toehold switches** are RNA-based regulatory devices that control gene expression in response to specific trigger RNAs. They consist of:
- A **toehold region** that binds the trigger RNA
- A **switching stem** that unfolds upon trigger binding
- A **ribosome binding site (RBS)** and **start codon** that become accessible
- A **universal linker** for consistent behavior
- A **downstream gene** (e.g., toxin) under translational control

When the correct trigger RNA binds, it exposes the RBS and AUG, enabling translation of the controlled gene.

## 📁 File Structure

```
igem_toehold/
├── trigger_mrna_generator.py    # Core trigger RNA design algorithm
├── design_toehold.py           # Complete toehold switch design pipeline
├── plot_nucleotide.py          # RNA secondary structure visualization
├── test_trigger.py             # Simple trigger RNA testing script
└── README.md                   # This documentation
```

## 🚀 Quick Start

### Prerequisites
```bash
# Install required packages
pip install vienna-rna biopython matplotlib numpy

# Or on Ubuntu/Debian:
sudo apt install vienna-rna

# Or on macOS:
brew install vienna-rna
```

### Basic Usage

#### 1. Generate Trigger RNA
```python
from trigger_mrna_generator import trigger_candidate

# Find optimal trigger RNA for a toehold switch
result = trigger_candidate(36, "GGAAGGAGGUAACAAUG")
if result:
    score, sequence, mfe, gc, mfe_freq, ensemble_div, no_repeats, rbs_exposed = result
    print(f"Best trigger: {sequence}")
    print(f"Score: {score:.2f}")
```

#### 2. Design Complete Toehold Switch
```python
from design_toehold import design_seriesB_toehold

# Design toehold switch with toxin gene
toxin_sequence = "AUGCGCGAAUUCGCG..."  # Your toxin gene sequence
designs = design_seriesB_toehold(toxin_sequence, trigger_length=25)

for design in designs:
    print(f"Trigger RNA: {design['trigger_mRNA']}")
    print(f"Complete sequence: {design['toehold_sequence']}")
```

#### 3. Visualize RNA Structure
```python
# Use plot_nucleotide.py to visualize toehold structures
# (Requires manual execution due to matplotlib GUI requirements)
python plot_nucleotide.py
```

## 🔧 Detailed Usage

### trigger_mrna_generator.py

**Core algorithm for designing trigger RNAs with optimal properties:**

- **GC content**: 40-55% for proper stability
- **No long repeats**: Avoid 4+ consecutive identical nucleotides
- **Thermodynamic stability**: MFE between -2.0 and 0 kcal/mol
- **Structural ensemble**: >30% MFE structure frequency
- **RBS accessibility**: Ensures RBS and AUG are exposed when trigger binds

**Command Line Interface:**
```bash
python trigger_mrna_generator.py 36 --switch GGAAGGAGGUAACAAUG
```

**Three Design Modes:**
- `main()`: Generate 10 top candidates (command-line)
- `trigger_candidate()`: Return single best candidate
- `trigger_candidate_fast()`: Return first valid candidate

### design_toehold.py

**Complete toehold switch design pipeline:**

1. **Trigger RNA Generation**: Uses `trigger_candidate_fast()` for speed
2. **Toehold Construction**: Combines trigger complement + conserved sequence + linker + toxin
3. **Structure Validation**: Verifies conserved sequence and linker regions are single-stranded
4. **Interaction Check**: Ensures no unwanted base-pairing between toehold and toxin

**Key Components:**
- **Conserved sequence**: `GGACUUUAGAACAGAGGAGAUAAAGAUG` (25 nt)
- **Universal linker**: `AACCUGGCGGCAGCGCAAAAG` (22 nt)
- **Toehold region**: Reverse complement of trigger RNA

### plot_nucleotide.py

**RNA secondary structure visualization:**

- **Stem-loop structures**: Displays RNA folding patterns
- **Annotated regions**: Labels toehold, linker, and conserved sequences
- **Nucleotide numbering**: Shows position information
- **Interactive display**: Matplotlib-based visualization

## 📊 Output Interpretation

### Trigger RNA Scores
Scores range from 0-20+, with components:
- **MFE stability**: Optimal around -1 kcal/mol (40% weight)
- **GC content**: Optimal at 47.5% (20% weight)
- **Ensemble diversity**: Normalized by sequence length (20% weight)
- **MFE frequency**: >30% in thermodynamic ensemble (20% weight)

### Structure Validation
- **Conserved sequence**: Must be single-stranded (all `.` in dot-bracket)
- **Universal linker**: Must be single-stranded for consistent behavior
- **Toxin independence**: No base-pairing between toehold and toxin regions

## 🔬 Biological Constraints

### RNA Design Rules
- **GC balance**: 40-55% for proper folding without over-stability
- **No homopolymers**: Avoid 4+ consecutive identical nucleotides
- **Thermodynamic window**: MFE -2.0 to 0 kcal/mol for switchable structures
- **Structural diversity**: Minimum ensemble diversity for robust behavior

### Sequence Elements
- **RBS pattern**: `AGGAGG` (Shine-Dalgarno sequence)
- **Start codon**: `AUG` for translation initiation
- **Toehold length**: Typically 25-36 nucleotides
- **Conserved region**: 25 nt for switch mechanism
- **Linker length**: 22 nt for standardization

## 🛠️ Customization

### Modifying Design Parameters
Edit constants in `trigger_mrna_generator.py`:
```python
GC_MIN = 0.40      # Minimum GC content
GC_MAX = 0.55      # Maximum GC content
DG_MIN = -2.0      # Minimum free energy
DG_MAX = 0.0       # Maximum free energy
```

### Adding New Toxin Genes
Modify `design_seriesB_toehold()` in `design_toehold.py`:
- Update toxin sequence validation rules
- Adjust structure requirements if needed
- Modify output file paths as required

## 📈 Performance Tips

- **Speed vs Quality**: Use `trigger_candidate_fast()` for quick results, `trigger_candidate()` for optimal designs
- **Batch Processing**: Generate multiple designs and screen for best performance
- **Parallelization**: Consider parallel processing for large-scale design campaigns
- **Structure Validation**: Always validate final designs experimentally

## 🔍 Troubleshooting

### Common Issues
- **RNAfold not found**: Install ViennaRNA package
- **No valid candidates**: Try different switch sequences or adjust length
- **Poor scores**: Check GC content and sequence complexity
- **Structure validation fails**: Verify toxin sequence compatibility

### Debug Mode
Enable detailed output by modifying print statements in the source files.

## 📚 References

- Green, A. A., et al. "Toehold switches: De-novo-designed regulators of gene expression." *Cell* 159.4 (2014): 925-939.
- iGEM Synthetic Biology resources and protocols
- ViennaRNA package documentation

## 👥 Contributing

This toolkit is designed for iGEM teams working with toehold switches. Contributions welcome for:
- New design algorithms
- Additional validation methods
- Visualization improvements
- Performance optimizations

## 📄 License

Open source for educational and research use in synthetic biology applications.

---

**Happy designing! 🧬🔬**
