# Pop Bio Global Health: Daily Schedule

**Monday, January 8: Introduction to population biology of infectious diseases**


* 9:30-9:45am: Introductions 
* 9:45-9:50am: Overview of [Syllabus](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/course-info/Syllabus-PopBio-ID-GlobalHealth-Jan2024.pdf)
* 9:50-10:50am: [Lecture: Introduction to the population biology of infectious diseases](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/pdf-lectures/IntroPopBioGlobalHealth.pdf)
  * Includes:
      * overview of causative agents in infectious disease
      * emphasis on diseases of importance in global health
      * differences between population biology approach and classic epidemiology
      * what is a model?
      * introduction to the R programming language
* 10:50-11:10am: Activity: Formulating research questions for infectious disease modeling - part 1
  * **For HW, start to brainstorm a disease of interest and a question of your own for an in-class activity on Wednesday.**
* 11:10-11:45am (continue for homework until you are comfortable): [Computer Tutorial:  Intro to R and R Studio](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/tutorials/Intro-R.zip)
  * Includes:
      * Knowing your working environment
      * Assigning variables
      * Basic arithmetic
      * Running a script
* [Here](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/assignment-templates/Reading-Digest-Template.pdf) are your instructions for Reading Digests due at the beginning of each class day. **All assignments should be uploaded to Canvas by 9:00am on the due date specified.**
* [Here](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/tree/main/papers) is a link to the folder with all papers for the course.
  

**Tuesday, January 9: Understanding compartmental models of infectious diseases**
* 9:30-9:50am: Reading recap
* 9:50-11:00am: [Lecture: Understanding compartmental models of infectious diseases](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/pdf-lectures/Intro_CompartmentalModels.pdf)
  * Includes:
    * discrete vs. continuous time
    * deterministic vs. stochastic models
    * review of basic differential equations and application to disease modeling
    * R0, RE
    * critical community size, herd immunity, critical vaccination threshold
* 11:00-11:45am: [Epidemic Cards](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/activities/Epidemic_Cards_Activity/Directions-Epidemic_Card-Game.pdf) modeling of an epidemic curve
* [Homework due on Thursday, January 11](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/assignment-templates/HW-Jan11-Disease-Q-States-Processes.pdf): Choose a disease of interest, and formulate a research question about it that can be addressed with a dynamical (compartmental) model. List the ‘states’ and ‘processes’ associated with your research question. **Have a draft of your disease and question ready for an in-class activity tomorrow (Wednesday)!**

**Wednesday, January 10: Thresholds to persistence in infectious diseases**

* 9:30-9:50am: Reading recap
* 9:50-10:50am: [Computer Tutorial: Visualizing and modeling data from Epidemic Cards](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/activities/Epidemic_Cards_Activity/Epidemic_Cards_Modeling.zip)
* 11:00-11:45am: Activity: Formulating research questions for infectious disease modeling - part 2
* [Homework due on Thursday, January 11](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/assignment-templates/HW-Jan11-Disease-Q-States-Processes.pdf): Choose a disease of interest, and formulate a research question about it that can be addressed with a dynamical (compartmental) model. List the ‘states’ and ‘processes’ associated with your research question.


**Thursday, January 11: Modeling interventions in infectious disease control**

* 9:30-9:40am: Reading recap
* 9:40-10:50am: Activity: Refining your research questions - part 3
* 10:50-10:55am: [Homework due on Monday, January 15](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/assignment-templates/HW-Jan15-Model-Diagram.pdf): Build a model diagram for your disease of interest and define its states and processes.
* 10:55-11:45am: [Activity: Dynamical Fever](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/activities/Dynamical_Fever/Dynamical_Fever_Download.zip) (group exercise and discussion)

**Friday, January 12: Course field trip to AIGHD and Museum Vrolik**

**Monday, January 15: Transmission dynamics and interventions for HIV**

* 9:30-9:50am: Reading recap
* 9:50-10:50am: Model Telephone
* 10:50-11:10am:  [Lecture: Transmission dynamics and interventions for HIV](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/tree/main/pdf-lectures)
	* Includes:
	    * history of HIV, with emphasis on past and present dynamics in LMICs, particularly Sub-Saharan Africa
	    * application of transmission models to understanding and intervening in HIV spread
	    * coupled dynamics of HIV and TB in LMICs
* 11:10-11:45am: Computer Tutorial: Compartmental modeling of HIV in Harare via Shiny App
  * You will need the most recent version of the ICI3D R package for this tutorial. If you are using your laptop for the tutorials, please install the package before you begin by running the command in R Studio:
  ```
  install.packages('devtools')
  devtools::install_github('ICI3D/ici3d-pkg')
  ```
  * To run the tutorial (after installing the package), type ```ICI3D::hivTutorial()```
* [Homework due on Tuesday, January 16](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/assignment-templates/HW-Jan16-Updated-Model-Diagram.pdf): Refine your model diagram after today’s ‘Model Telephone’ activity.


**Tuesday, January 16: Transmission dynamics and interventions for vector-borne diseases**

* 9:30-9:50am: Reading recap
* 9:50-10:10am:  Computer Tutorial: Compartmental modeling of HIV in Harare via Shiny App
  * You will need the most recent version of the ICI3D R package for this tutorial. If you are using your laptop for the tutorials, please install the package before you begin by running the command in R Studio:
  ```
  install.packages('devtools')
  devtools::install_github('ICI3D/ici3d-pkg')
  ```
* 10:10-10:50am:  [Lecture: Transmission dynamics and interventions for vector-borne diseases](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/pdf-lectures/IntroVBD.pdf)
	* Includes:
	    * biology of malaria and impact on human evolution
	    * models of malaria and challenges to control
	    * biology of dengue and challenges to control
	    * brief overview of other VBDs
	    * impacts of climate change on VBD transmission
* 10:50-11:45am: Activity: Writing equations for a model world
  * Take 15 minutes on your computer or iPad, and write out the equations that you think best represent your compartmental model.
  * Once complete, we will go over a subset of these as examples in class.
* [Homework due on Wednesday, January 17](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/assignment-templates/HW-Jan17-Model-Equations.pdf): Write out equations for your model.


**Wednesday, January 17: Transmission dynamics and interventions to control zoonotic diseases (case study: rabies)**

* 9:30-9:50am: Reading recap
* 9:50-10:35am:  [Lecture: Transmission dynamics and interventions to control zoonotic diseases](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/pdf-lectures/IntroZoonoticDiseases.pdf)
	* Includes:
	  * key terms in understanding zoonotic diseases: reservoir vs. spillover hosts, maintenance vs. target populations 
	  * reservoir culling and vaccination 
	  * population biology of rabies persistence and elimination
* 10:35-11:45am: [Computer Tutorial: Building compartmental models for infectious diseases in R](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/tutorials/CompartmentalModels.zip)
  * We will start with "TutorialSIRModels.R" as this teaches you how to build a disease model in R. This is similar to what you need to turn in for your homework Friday. Please note that your model does not need to be working by Friday--you just need to demonstrate an attempt to write the code to get credit for the assignment. We will help you edit it until it is working.
  * If you have extra time or are simply interested, work through the second script "BonusTutorialCompartmentalModelsExpanded.R" at your leisure.
* [Homework due on Friday, January 19](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/assignment-templates/HW-Jan19-Draft-R-Script.pdf): Draft R script of your infectious disease model


**Thursday, January 18: Transmission dynamics and interventions for diseases with complex life cycles (case studies: schistosomiasis, cholera)**

* 9:30-9:50am: Reading recap
* 9:50-10:20am:  [Brief Lecture: Nutrition, immunity, and poverty in the control of infectious diseases](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/pdf-lectures/IntroComplexLifeCycles.pdf)
	* Includes:
	  * transmission and control of schistosomiasis
	  * hygiene hypothesis
	  * vaccination, immunity, and nutrition
	  * poverty traps and coupled social-ecological models
* 10:20-10:40am: [Building an example ODE for our classroom example (all together)](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/tutorials/ExampleModel.zip)
* 10:40-11:15am: [Computer Tutorial: Model evaluation and comparison](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/tutorials/Age-Prevalence-Model-Comp.zip)
* 11:15-11:45am: Mentored work time on your own projects. WHen you finish above, use the available class time to start coding your own model, and seek help as needed from Gwen and Cara.
* [Homework due on Friday, January 19](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/assignment-templates/HW-Jan19-Draft-R-Script.pdf): Draft R script of your infectious disease model

**Friday, January 19: Course field trip to Edible Insect Farming Working Group at the Université de Strasbourg and vist to the 'In the Days of AIDS' exhibit at the Museum of Contemporary and Modern Art in Strasbourg.**


**Monday, January 22: Introduction to phylogenetic models**

* 9:30-9:50am: Reading recap
* 9:50-10:50am:  Lecture: Introduction to phylogenetics (Gwen)
* 10:50-11:45am: Computer Tutorial: [Building and interpreting simple maximum likelihood phylogenies using MEGA and RAxML](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/tutorials/Intro_phylogenetics.zip) (Gwen)
* Programs to download before class: 
  * [MEGA](https://www.megasoftware.net)
  * [RAxML](https://antonellilab.github.io/raxmlGUI/)
  * [FigTree](https://github.com/rambaut/figtree/releases) - this is not required before class but it could be useful to use on your own
  * In R please install the following packages
     ```
      install.packages("BiocManager")
      install.packages('ggtree')
      install.packages('tidyverse')
      install.packages('ggplot2')
      install.packages('cowplot')
      install.packages('ape')
      install.packages('ggnewscale')

     ```
* Homework due on Thursday, January 25: [Final R script + term paper](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/assignment-templates/Final-Paper-Model.pdf) + [final presentations](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/assignment-templates/Final-Presentation.pdf)


**Tuesday, January 23: Molecular epidemiology of SARS-CoV-2**

* 9:30-9:50am: Reading recap
* 9:50-10:50am:  [Lecture: Introduction to phylodynamics and molecular epidemiology](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/pdf-lectures/IntroPhylodynamics.pdf)
	* Includes:
	  * history of phylodynamics as a field
	  * overview of viral genomics: structure and functions
	  * zoonotic origins of SARS-CoV-2
* 10:50-11:45am: [Computer Tutorial: Building and interpreting TimeTrees using BEAST](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/tutorials/BEAST-SIR.zip)
  * Prior to class, please download [BEAST2](https://www.beast2.org) and download and install [Tracer](https://github.com/beast-dev/tracer/releases/tag/v1.7.2) on your computer
* Homework due on Thursday, January 25: [Final R script + term paper](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/assignment-templates/Final-Paper-Model.pdf) + [final presentations](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/assignment-templates/Final-Presentation.pdf)


**Wednesday, January 24: Next Generation Sequencing and Emerging Infectious Diseases**

* 9:00-9:20am: Reading recap
* 9:20-9:50am:  [Lecture: Next Generation Sequencing and Emerging Infectious Diseases](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/pdf-lectures/IntroSequencing.pdf)
	* Includes:
	  * Sanger vs. Illumina vs. nanopore sequencing
	  * introduction to Nextstrain and GISAID
	  * global expansion of pathogen genomic sequencing in response to COVID-19
	  * realtime outbreak response from NGS data
zoonotic origins of SARS-CoV-2
* 9:50-11:45am: [Activity: Investigating novel viruses from Next Generation Sequencing data](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/activities/NGS-Outbreak-Investigation.pdf)
  * [index case sequencing](https://artic-network.github.io/artic-live/gc/index_case_title.html)
  * heavy metal associated cases [one](https://artic-network.github.io/artic-live/gc/case_2_title.html), [two](https://artic-network.github.io/artic-live/gc/case_3_title.html), [three](https://artic-network.github.io/artic-live/gc/case_4_title.html)
  * heavy metal associated controls [one](https://artic-network.github.io/artic-live/gc/control_1_title.html), [two](https://artic-network.github.io/artic-live/gc/control_2_title.html), [three](https://artic-network.github.io/artic-live/gc/control_3_title.html)
  * [pet shop environmental sampling](https://artic-network.github.io/artic-live/gc/pet_shop_env_title.html)
  * [wastewater sampling](https://artic-network.github.io/artic-live/gc/wastewater_seq_title.html)
  * [Genome Segment S](https://artic-network.github.io/artic-live/gc/segment_S.fasta.txt)
  * [Genome Segment L](https://artic-network.github.io/artic-live/gc/segment_L.fasta.txt)
  * [Segment L Phylogeny](https://nextstrain.org/community/emmahodcroft/GC/arenavirus/L)
  * [Segment S Phylogeny](https://nextstrain.org/community/emmahodcroft/GC/arenavirus/S)
  * [Tanglegram](https://nextstrain.org/community/emmahodcroft/GC/arenavirus/S:community/emmahodcroft/GC/arenavirus/L)
  * [Metal cluster sequencing](https://nextstrain.org/community/emmahodcroft/GC/HMFV/FC)
  * [Phylogenetic analysis of rock concert sequences](https://nextstrain.org/community/emmahodcroft/GC/HMFV/FC?p=grid)
  * [Final Nextstrain phylogeny](https://nextstrain.org/community/emmahodcroft/GC/HMFV/SC?p=grid)
* Homework due on Thursday, January 25: [Final R script + term paper](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/assignment-templates/Final-Paper-Model.pdf) + [final presentations](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/assignment-templates/Final-Presentation.pdf)

**Thursday, January 25: Final Presentations**

*  [Final R script + term paper](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/assignment-templates/Final-Paper-Model.pdf) due by 9:30 am
* 9:30-11:45am: [Final Presentations](https://github.com/brooklabteam/BIOS27816-pop-bio-global-health/blob/main/assignment-templates/Final-Presentation.pdf) (10 minutes per student + 5 minutes for questions) 
