---
title: ' `medRCT`: Causal mediation analysis estimating interventional effects mapped
  to a target trial in R'
tags:
- R
- causal inference
- mediation analysis
- interventional effect
- target trial
date: "11 January 2025"
affiliations:
- name: Clinical Epidemiology and Biostatistics Unit, Murdoch Children’s Research
    Institute, Australia
  index: 1
- name: Melbourne Dental School, University of Melbourne, Australia
  index: 2
- name: Department of Paediatrics, University of Melbourne, Australia
  index: 3
authors:
- name: Tong Chen
  orcid: "0000-0003-2751-2412"
  affiliation: 1, 2
  corresponding: true
- name: S. Ghazaleh Dashti
  affiliation: 1, 3
- name: "Margarita Moreno-Betancur"
  affiliation: 1, 3
bibliography: ../inst/references.bib
---

# Summary
Causal mediation analysis is often used to explore the causal pathways linking 
exposures to outcomes across various research domains, including epidemiology, 
public health, and social sciences. It quantifies the extent to which the causal effect of an 
exposure on an outcome is mediated through one or more intermediate variables, known as mediators. 
The goal of studies conducting mediation analyses, particularly in health research, is often to understand 
how hypothetical interventions targeting mediators might counter exposure effects in order to guide 
future treatment and policy decisions as well as intervention development.

The interventional effects [@Geneletti2007Identifying; @vanderweele2014effect] are causal mediation estimands that overcome certain
identifiability challenges associated with traditional estimands, known as the natural effects [@Robins1992]. 
Unlike the natural effects, the interventional effects implicitly emulate target trials [@hernan2016using] 
that assess the effects of distributional shifts in mediators [@MorenoBetancur2018]. As a result, interventional effects are 
more closely aligned with the goal of examining the potential impacts of hypothetical interventions that shift mediators. 
However, the original definitions of interventional effects do not necessarily assess shifts that are relevant to 
policy and practice. To bridge this gap, @Moreno2021Mediation proposed definitions of interventional effects that map 
explicitly to a target trial assessing policy-relevant hypothetical mediator interventions.

`medRCT` provides a user-friendly interface for estimating policy-relevant interventional effects as
defined in @Moreno2021Mediation using a Monte Carlo simulation-based g-computation approach. 
Key features of `medRCT` include support for multiple mediators, 
intermediate confounders (exposure-induced mediator-outcome confounders), and a Shiny 
application [@chang2024shiny] for comprehensive model assessment. Detailed methods and estimands 
implemented in `medRCT` are described in @Moreno2021Mediation.

# Statement of need

Several software packages are available for causal mediation analysis. Existing tools, 
such as `mediation` [@Tingley2014mediation] and `medflex` [@Steen2017medflex], focus primarily on 
natural effects estimands, which are prone to identification challenges and rely on assumptions that 
are difficult to validate in real-world studies [@Vansteelandt2017]. While software like 
`medoutcon` [@hejazi2022medoutcon] leverages causal machine learning methods to estimate 
interventional effects, they do not consider estimands defined explicitly in terms of a target 
trial examining policy-relevant interventions. This limits their direct practical relevance in 
informing decision-making. 

The `medRCT` package addresses these gaps by estimating policy-relevant interventional effects, explicitly mapped to a 
target trial. There are several published applications of this methodology in real-world studies, including @Dashti2022, 
@Goldfeld2023, and @Afshar2024, demonstrating its utility across diverse research contexts. For instance, @Dashti2022 studied 
the mediating roles of C-reactive protein, leptin, fasting insulin, and estradiol in the effect 
of adiposity on cancer risk among postmenopausal women. @Goldfeld2023 assessed the extent to which 
inequities could be mitigated by improving disadvantaged children’s parental mental health and 
preschool attendance. Before the development of `medRCT`, these studies had to rely on user-written code, which was prone 
to errors and limited both accessibility and reproducibility. Several ongoing studies at 
the Murdoch Children’s Research Institute are now using `medRCT`. For example, one study aims to investigate 
the impact on cardiovascular outcomes of a hypothetical intervention that shifts the distribution of inflammatory markers 
in adolescents from high-income households to the levels in those from low-income households, 
using data from multiple longitudinal cohort studies. 

`medRCT` has also been used in education and training. It has been used in workshops and 
training sessions, such as the [ViCBiostat Summer School](https://www.vicbiostat.org.au/event/summer-school-2024-causal-mediation-analysis) 
and workshop at the Society for Epidemiologic Research 2024 meeting.

The package has been extensively documented on its [website](https://t0ngchen.github.io/medRCT/). Future developments for `medRCT` 
include extending its functionality to handle survey design weights to broaden its applicability to real-world studies involving
complex survey designs.




# Acknowledgements

This work was supported by an Investigator Grant fellowships to MMB [grant ID: 2009572] 
and SGD [grant ID: 2027171] from the National Health and Medical Research Council. 
The Murdoch Children’s Research Institute is supported by the Victorian Government's 
Operational Infrastructure Support Program.

# References
