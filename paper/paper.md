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
Causal mediation analysis is often used to explore the causal pathways linking exposures to outcomes across various research domains, including epidemiology, public health, and social sciences. It quantifies the extent to which the causal effect of an exposure on an outcome is mediated through one or more intermediate variables, known as mediators. The goal of studies conducting mediation analyses, particularly in health research, is often to understand how hypothetical interventions targeting mediators might counter exposure effects in order to guide future treatment and policy decisions as well as intervention development. 


The modern causal inference literature initially defined the natural (in)direct effects [@Robins1992; @Pearl2001] as the estimands of interest in causal mediation analysis. Grounded in the potential outcomes framework, these effects are defined based on cross-world counterfactuals [@Robins2011], and their identifiability relies on a cross-world independence assumption. This assumption, however, can never be guaranteed to hold, not even with an experiment [@Robins2011; @Didelez2006]. Further, this assumption is not met in the presence of a mediator-outcome confounder that is itself affected by the exposure [@Avin2005; @vanderweele2014effect; @Vansteelandt2012Natural], and more generally, in the presence of multiple mediators of interest [@VanderWeele2014]. 


Interventional effects [@Geneletti2007Identifying; @vanderweele2014effect] have been proposed as an alternative to address these limitations. In particular, these effects have been shown to implicitly emulate target trials [@hernan2016using] that assess the impact of hypothetical interventions shifting the distribution of the mediators [@MorenoBetancur2018]. It has thus been proposed that interventional effects be explicitly defined by mapping to a hypothetical randomised trial (a target trial) that assesses the hypothetical interventions of interest [@Moreno2021Mediation]. Specifying the target trial clarifies the causal question and makes the interventional effects directly policy-relevant and practically meaningful.


@Moreno2021Mediation proposed three interventional effects, with each one examining a distinct policy-relevant question about how hypothetical interventions that would shift mediator distributions individually, together, or in sequence might impact the outcome. For example, they examined the extent to which the increased risk of financial hardship in mid-adulthood (outcome) among adolescent self-harmers (exposed group) could be countered by a hypothetical intervention that would shift the distribution of weekly cannabis use in young adulthood (mediator) among adolescent self-harmers to the levels of those who did not self-harm (unexposed group), either treating other mediators as independent or allowing for flow-on effects on its causal descendants, such as education and employment. They also considered an interventional effect capturing the effect of an intervention that shifts the joint distribution of all mediators (in their example, depression or anxiety, cannabis use in young adulthood, education, and employment status) among adolescent self-harmers to the levels in those who did not self-harm.

The `medRCT` package provides a user-friendly interface for estimating policy-relevant interventional effects as defined in @Moreno2021Mediation using a Monte Carlo simulation-based g-computation approach. Key features of `medRCT` include support for multiple mediators, intermediate confounders (exposure-induced mediator-outcome confounders), and a Shiny application [@chang2024shiny] for comprehensive model assessment. Detailed definition of the mediation estimands that can be estimated using the `medRCT` package, as well as the target trial to which they map, their identifiability and estimation procedures are described in @Moreno2021Mediation.

# Statement of need

Several R [@R] software packages are available for causal mediation analysis. The `mediation` package [@Tingley2014mediation] implements the estimation of natural effects using a generic approach based on Monte-Carlo integration methods as described by @Imai2010. The package can also estimate path-specific effects when there are multiple mediators, and allows researchers to conduct sensitivity analyses to evaluate the robustness of their results to potential unmeasured confounding. The `medflex` package [@Steen2017medflex] implements the methods proposed by @Vansteelandt2012Imputation, which directly model the natural effects based on a class of mean models for nested counterfactuals. The `medoutcon` package [@hejazi2022medoutcon] implements asymptotically efficient causal machine learning-based estimators of both the natural and interventional direct and indirect effects [@Díaz2020]. However, the interventional effects considered are not defined explicitly in terms of a target trial examining policy-relevant interventions, which may limit their direct practical relevance in informing decision-making. 


The `medRCT` package addresses these gaps by estimating policy-relevant interventional effects, explicitly mapped to a target trial. However, as `medRCT` uses a Monte Carlo simulation-based g-computation approach, it requires a large number of replications to produce stable and reliable inference, which can be computationally intensive. Additionally, the approach is sensitive to model misspecification, as all nuisance parameters can only be estimated via restrictive parametric modelling. Furthermore, the method is not suited for settings with high-dimensional mediators or intermediate confounders, where causal machine learning-based estimators are required [@Díaz2020; @Rudolph2024; @liu2024general; @chen2025causal].


There are several published applications of this methodology in real-world studies, including @Dashti2022, @Goldfeld2023, and @Afshar2024, demonstrating its utility across diverse research contexts. For instance, @Dashti2022 studied the mediating roles of C-reactive protein, leptin, fasting insulin, and estradiol in the effect of adiposity on cancer risk among postmenopausal women. @Goldfeld2023 assessed the extent to which inequities could be mitigated by improving disadvantaged children’s parental mental health and preschool attendance. Before the development of `medRCT`, these studies had to rely on user-written code, which is prone to errors and limits both accessibility and reproducibility. Several ongoing studies at the Murdoch Children’s Research Institute are now using `medRCT`. For example, one study aims to investigate the impact on cardiovascular outcomes of a hypothetical intervention that shifts the distribution of inflammatory markers in adolescents from high-income households to the levels in those from low-income households, using data from multiple longitudinal cohort studies. 

`medRCT` has also been used in education and training. It has been used in workshops and training sessions, such as the [ViCBiostat Summer School](https://www.vicbiostat.org.au/event/summer-school-2024-causal-mediation-analysis) and workshop at the Society for Epidemiologic Research 2024 meeting.


The package has been extensively documented on its [website](https://t0ngchen.github.io/medRCT/). Future developments for `medRCT` include extending its functionality to handle survey design weights to broaden its applicability to real-world studies involving complex survey designs.


# Acknowledgements


This work was supported by an Investigator Grant fellowships to MMB [grant ID: 2009572] and SGD [grant ID: 2027171] from the National Health and Medical Research Council. The Murdoch Children’s Research Institute is supported by the Victorian Government's Operational Infrastructure Support Program.


# References
