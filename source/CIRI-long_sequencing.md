# Nanopore library preparation for CIRI-long method

## 1. Total RNA Extraction & Ribosomal RNA Depletion

In our manuscript:

- Total RNA is isolated using TRIzol (Invitrogen)
- **RiboErase kit (human/mouse/rat, KAPA Biosystems)** is used to remove rRNA from **1ug of total RNA**.

## 2. Poly(A) Tailing & RNase R Treatment

Then, additional poly(A) tails are added to the linear transcripts to increase RNase R digestion efficiency.

### 2.1 Add poly(A) tails to linear RNAs

E-PAP treatment is used to add poly(A) tails to the 3' end of linear RNAs, which can increase the RNase R digestion ability to RNAs with secondary structures.

Component | Volume
-|-
RNA | 17 uL
PolyA Buffer | 2 uL
NEB *E.coli* Poly(A) Polymerase (5 U/uL) | 1 uL
**Total Volume:** | **20 uL**

- Ribosomal-depleted total RNA is incubated at 37ºC with 1uL of Poly(A) Polymerase for 30min.

### 2.2 Purification after poly(A) treatment

AMPure XP is used to remove contamination after poly(A) treatment.

Component | Volume
-|-
Agencourt AMPure XP magnetic beads (Beckman) | 44 uL
Polyadenylated RNA | 20 uL
**Total Volume:** | **64 uL**

- Add 2.2uL of AMPure XP per 1.0 uL of sample.
- Wash beads + RNA fragments twice with 75% Ethanol to remove contaminants.
- Elute purified RNA with 20uL H2O.

### 2.3 RNase R Treatment to effectively digest linear RNAs

Polyadenylted RNA is treated using RNase R to remove linear RNAs.

Component | Volume
-|-
Polyadenylated RNA  | 17.5 uL
RNase R Buffer | 2 uL
Epicentre RNase R (20 U/uL) | 0.5 uL
**Total** | **20 uL**

- Polyadenylated RNA was incubated with RNase R at 37ºC for 15 min.

### 2.4 Purification after RNase R Treatment

AMPure XP is used to remove contamination after RNase R treatment.

Component | Volume
-|-
AMPure XP | 44 uL
Polyadenylated RNA | 20 uL
**Total Volume:** | **64 uL**

- Add 2.2uL of AMPure XP per 1.0 uL of sample.
- Wash beads + RNA fragments twice with 75% Ethanol to remove contaminants.
- Elute purified RNA with 5uL H2O. 

## 3. SMARTer Reverse Transcription

Then, RNase R-treated RNA is reverse transcribed using random hexamers and SMARTer cDNA synthesis kit (Takara Bio) according to the manufacturer's instructions. The 3' SMART CDS Primer II A `5'-AAGCAGTGGTATCAACGCAGAGTACT(30)N-1N-3'` was replaced with `5'-AAGCAGTGGTATCAACGCAGAGTACNNNNNN-3'` to amplify circular RNAs without poly(A) sequences.

- Prepare reaction as follow:

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

To obtain sufficient cDNA products for sequencing, PCR amplification is performed using **2uL** of cDNA with **NEBNext LongAmp Taq DNA Polymerase** and **SMARTer primers** under the following conditions:

STEP | TEMP | TIME
-|-|-
Initial Denaturation | 95ºC | 30 s
19-20 Cycles | 95ºC<br>62ºC<br>65ºC | 15 s<br>15 s<br>2 min 
Final Extension | 65ºC | 2 min
Hold | 4-10ºC

## 5. Fragment Size selection

Afterward, AMPure XP is used for size selection of the cDNA fragments:

- Add 0.5 of AMPure XP per 1.0 uL of sample.
- Wash beads + cDNA fragments twice with 75% Ethanol to remove contaminants.
- Elute purified cDNA

## 6. Nanopore Sequencing

Finally, cDNA libraries are prepared according to the ONT protocol `SQK-LSK109` and barcoded with `EXP-NBD104` / `EXP-NBD114` kits, and nanopore sequencing is performed using the MinION (MN26543) platform with a `FLOW-MIN106` flow cell.
