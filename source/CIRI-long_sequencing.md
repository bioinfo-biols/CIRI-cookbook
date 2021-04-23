# CIRI-long Nanopore Library Preparation

## 1. Total RNA Extraction & Ribosomal RNA Depletion

- Total RNA is isolated using [TRIzol (Invitrogen)](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/trizol_reagent.pdf)

- [RiboErase kit (human/mouse/rat, KAPA Biosystems)](https://rochesequencingstore.com/wp-content/uploads/2017/10/KAPA-RiboErase-KitHMR_KR1142-%E2%80%93-v4.19.pdf) is used to remove rRNA from **1ug of total RNA**.

- Elute rRNA-depleted RNA from beads with 17uL nuclease free water.

## 2. Poly(A) Tailing & RNase R Treatment

Then, additional poly(A) tails are added to the linear transcripts to increase RNase R digestion efficiency.

### 2.1 Add Poly(A) Tails To Linear RNAs

E-PAP treatment is used to add poly(A) tails to the 3' end of linear RNAs, which can increase the RNase R digestion ability to RNAs with secondary structures.

- Add the following components in the order specified:

Component | Volume
-|-
rRNA-depleted RNA | 15 uL
10X *E.coli* Poly(A) Polymerase Reaction Buffer | 2uL
ATP (10mM) | 2 uL
*E.coli* Poly(A) Polymerase (5 U/uL) | 1 uL
**Total Volume:** | **20 uL**

- Ribosomal-depleted total RNA is incubated at 37ºC with 1uL of [*E.coli* Poly(A) Polymerase (NEBNext)](https://international.neb.com/protocols/2014/08/13/poly-a-tailing-of-rna-using-e-coli-poly-a-polymerase-neb-m0276) for 30min.
- Stop the reaction by proceeding to the cleanup step.

### 2.2 Purification After Poly(A) Treatment

AMPure XP is used to remove contamination after poly(A) treatment.

Component | Volume
-|-
AMPure XP | 44 uL
Polyadenylated RNA | 20 uL
**Total Volume:** | **64 uL**

- Add 2.2uL of [Agencourt AMPure XP magnetic beads (Beckman)](https://www.beckmancoulter.com/wsrportal/techdocs?docname=B37419) per 1.0 uL of sample.
- Wash beads + RNA fragments twice with 75% Ethanol to remove contaminants.
- Elute purified RNA from beads with 20uL H2O.

### 2.3 RNase R Treatment To Effectively Digest Linear RNAs

Polyadenylted RNA is treated using RNase R to remove linear RNAs.

Component | Volume
-|-
Polyadenylated RNA  | 17.5 uL
RNase R Buffer | 2 uL
RNase R (20 U/uL) | 0.5 uL
**Total** | **20 uL**

- Polyadenylated RNA was incubated with [RNase R (Epicentre)](https://www.lucigen.com/docs/manuals/MA266E-RNase-R.pdf) at 37ºC for 15 min.

### 2.4 Purification after RNase R Treatment

2.2x bead-based cleanup is used to remove contamination after RNase R treatment as described above (See 2.2).

- Add 2.2uL of AMPure XP per 1.0 uL of sample.
- Wash beads + RNA fragments twice with 75% Ethanol to remove contaminants.
- Elute purified RNA with 5uL H2O. 

## 3. SMARTer Reverse Transcription

Then, RNase R-treated RNA is reverse transcribed using random hexamers and [SMARTer cDNA synthesis kit (Takara Bio)](https://www.takarabio.com/documents/User%20Manual/SMARTer%20PCR%20cDNA%20Synthesis%20Kit%20User%20Manual%20%28PT4097/SMARTer%20PCR%20cDNA%20Synthesis%20Kit%20User%20Manual%20%28PT4097-1%29_040114.pdf) according to the manufacturer's instructions. The 3' SMART CDS Primer II A `5'-AAGCAGTGGTATCAACGCAGAGTACT(30)N-1N-3'` was replaced with `5'-AAGCAGTGGTATCAACGCAGAGTACNNNNNN-3'` to amplify circular RNAs without poly(A) sequences.

- Prepare reaction as follows:

Component | Volume
-|-
RNase R treated RNA | 3.5 uL
SMARTer CDS random primer (12 uM) | 1 uL
**Total Volume:** | **4.5 uL**

- Incubated at 72ºC for 3min, 25ºC for 10min, hold at 42ºC.

- Add the following mixture:

Component | Volume
-|-
5x First Strand Buffer (RNase-Free) | 2 uL
Dithiothreitol (DTT; 100 mM) | 0.25 uL
dNTP Mix (10 mM) | 1 uL
SMARTer II A Oligonucleotide (12 uM) | 1 uL
RNase Inhibitor (40 U/uL) | 0.25 uL
SMARTScribe Reverse Transcriptase (100 U/uL) | 1 uL
**Total Volume:** | **5.5 uL**

- Incubation at 42ºC for 90 min.
- Denatured at 70ºC for 10 min.

## 4. cDNA PCR Amplification

To obtain sufficient cDNA products for sequencing, PCR amplification is performed using **2uL** of cDNA with [NEBNext LongAmp Taq DNA Polymerase](https://international.neb.com/protocols/2012/10/15/m0323-longamp-taq-dna-polymerase-protocol) and **SMARTer primers** under the following conditions:

STEP | TEMP | TIME
-|-|-
Initial Denaturation | 95ºC | 30 s
19-20 Cycles | 95ºC<br>62ºC<br>65ºC | 15 s<br>15 s<br>2 min 
Final Extension | 65ºC | 2 min
Hold | 4-10ºC

## 5. Fragment Size Selection

Afterward, 0.5x AMPure XP is used for size selection of the cDNA fragments:

- Add 0.5uL of AMPure XP per 1.0 uL of sample.
- Wash beads + cDNA fragments twice with 75% Ethanol to remove contaminants.
- Elute purified cDNA with 20-30uL H2O.

## 6. Nanopore Sequencing

Finally, cDNA libraries are prepared according to the ONT protocol `SQK-LSK109` and barcoded with `EXP-NBD104` / `EXP-NBD114` kits, and nanopore sequencing is performed using the MinION (MN26543) platform with a `FLOW-MIN106` flow cell. Please refer to the [Nanopore Community](https://nanoporetech.com/) for detailed instructions.
