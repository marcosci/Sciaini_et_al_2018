;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                                                                                                                                                                                                           ;;
;;  SwiFCoIBM - Indivudal-Based Swine Fever Community Model (based on model versions published in Kramer-Schadt et al. 2009 and Lange et al. 2012)                                                                                                           ;;
;;                                                                                                                                                                                                                                                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                                                                                                                                                                                                           ;;
;;  The model is composed of two major components:                                                                                                                                                                                                           ;;
;;    (1) a wild boar demography model considering seasonal reproduction, herd splitting (dispersal) and mortality and                                                                                                                                       ;;
;;    (2) a Classical Swine Fever (CSF) virus model operating on the emerging wild boar population.                                                                                                                                                          ;;
;;                                                                                                                                                                                                                                                           ;;
;;  Wild boar population density and structure are affected by the disease via virus-induced mortality and litter size depression.                                                                                                                           ;;
;;  The crucial model entity is the wil boar individual, characterised by age in weeks (one week represents the approx. CSF incubation time; Artois et al. 2002, Moenning et al. 2003),                                                                      ;;
;;    resulting in age classes piglets (age <= 8 months), yearlings (subadults; 8 months < age <= 2 years) and adults (age > 2 years).                                                                                                                       ;;
;;  Each host has a location, which denotes its home range cell and the individual's family group (HerdID).                                                                                                                                                  ;;
;;                                                                                                                                                                                                                                                           ;;
;;  Movement: Adult movement is restrcited to the 'activity zone' which equals their home range cell (staying strategy). Movement is similar similar for individuals of the sounder (social unit organised around adult females and their offspring). ;;
;;            Subadult females split in specified weeks of the year from their herds searching for an unoccupied home range cell (see Procedure "HerdSplit").                                                                                                ;;
;;            Subadult males disperse from natal family groups in subgroups at an age the males begin to seek mates (~5 months), joining other groups in the surrounding area.                                                                               ;;                                                                                                                                    ;;
;;                                                                                                                                                                                                                                                           ;;
;;  Landscape: The model landscape is represented by a grid of 4 km^2 square cells which each encompass a wild boar family group's home range (Leaper et al. 1999).                                                                                          ;;
;;             Each habitat cell is characterised by a breeding capacity (Quality), denoting habitat quality by the number of female boars that are allowed to breed, representing density regulation in the model.                                          ;;
;;                                                                                                                                                                                                                                                           ;;
;;                                                                                                                                                                                                                                                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                                                                                                                                                                                                           ;;
;;  STATE VARIABLES                                                                                                                                                                                                                                          ;;
;;                                                                                                                                                                                                                                                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



; NAME                            TYPE/RANGE     DEFAULT        DESCRIPTION

globals
[
  #Landscape                  ;   [string]       -              way of landscape generation or import
  Habitat                     ;   [set]          -              habitat cells (cells with Quality > 0)
  Matrix                      ;   [set]          -              matrix cells (cells with Quality = 0)
  QualityRange                ;   [#]            9              maximum value for habitat quality to define range of habitat quality (quality between 0 and QualityRange)
  YearRelease                 ;   [#]            6              year in which an randomly individual becomes infected (range of WeekRelease)
  WeekRelease                 ;   [#]            261-312        random week in YearRelease of virus release following an expontential distribution with lambda = 1.5
  MaxAge                      ;   [#]            572 weeks      maximum age in weeks
  Longevity                   ;   [#]            11 years       maximum age in years (MaxAge / 52)
  AgeBlur                     ;   [#]            6 weeks        stochastic deviation from age class thresholds in weeks (in both directions (age class threshold +/- AgeBlur) following a normal distribution
  DispDist                    ;   [#]            3 cells        dispersal radius for subadults
  SurvivalProb                ;   [%]            0.65           annual survival rate of adults and subadults (mean = 0.65, min = 0.4)
  SurvivalProbPig             ;   [%]            0.5            annual survival rate of piglets (mean = 0.5, min = 0.1)
  currHerdID                  ;   [#]            1              current maximum herd ID
  CaseFatality                ;   [%]            0.8            general risk for getting lethally infected by a disease transmission
  BetaWithin                  ;   [#, 0-1]       0.0208         transmission coefficient within the herd
  BetaBetween                 ;   [#, 0-1]       0.00208        transmission coefficient between herds
  PreNatInf                   ;   [#, 0-1]       0.5            probability of pre-natal infection
  FertRedInf                  ;   [#, 0-1]       0.625          reduction in fertility of infected females
  mue                         ;   [#]            3 weeks        mean of exponential survival time distribution for lethally infected individuals
  T_trans                     ;   [#]            1 week         infectious period of transiently infected individuals
  T_anti                      ;   [#]            15 weeks       maximum persistence of maternal antibodies (Depner et al. 2000)
  ProbRepro                   ;   [#]                           monthly reproduction probability
  week                        ;   [#, 1-52]      0              current week of the year
  SimulatedYears              ;   [#]            12             number of simulated years during one run
  ProbFem                     ;   [%]            0.5            sex ratio (0: only males; 1: only females)
  NewLeth                     ;   [#]            0              number of individuals which got lethally infected during this timestep
  NewTrans                    ;   [#]            0              number of individuals which got transient infected during this timestep
  NmaxList                    ;   [list]         []             input data: sums of the number of patches adjacent to each patch in a circle of the same area as the infected/infectious area (see Report_fragmentation)
  ReproList                   ;   [list]         []             input data: cumulative reproduction probabilities for each month
  DONE                        ;   [boolean]      0              1 if stop conditions are TRUE

  ;; OUTPUT
  Time                        ;   [min]          0              counting minutes since starting setup
  Stopped                     ;   [string]       -              reason for stop of siumulation
  InfDist                     ;   [cells]        0              maximum spatial distance of disease transmission
  InfX                        ;   [#]            -              coordinate x of the pathogen-introduced cell
  InfY                        ;   [#]            -              coordinate y of the pathogen-introduced cell
  MaxDist?                    ;   [#]            -              week when infection reached other side of the tubus
  WeekLastInf                 ;   [#]            -              week in which the last infectious indiviual was noted (excl. No-Shedders)
  WeekLast                    ;   [#]            -              week in which the last infected indiviual was noted (incl. No-Shedders)
  F_infected                  ;   [0-1]          -              fragmentation index for infectious cells (patches which contain shedding individuals)
  F_infectious                ;   [0-1]          -              fragmentation index for infected cells (patches which were infected at least one week during the simulation)
  PopStructList               ;   [list]         []             current age class distribution of all individuals in years
  InfStructList               ;   [list]         []             current age class distribution of infected individuals (EpiStat = esLeth or esTrans) in years
  InfTransStructList          ;   [list]         []             current age class distribution of transient infected individuals (EpiStat = esTrans) in years
  InfLethStructList           ;   [list]         []             current age class distribution of lethally infected individuals (EpiStat = esLeth) in years
  CurrInfPeriodList           ;   [list]         []             current infectious periods of lethally infected individuals (EpiStat = esLeth) in weeks
  LethPeriodList              ;   [list]         []             all infectious periods of lethally infected individuals (EpiStat = esLeth) in weeks
  InfPeriodList               ;   [list]         []             all infectious periods of lethally and transient infected individuals (EpiStat = esLeth and esTrans) in weeks
]

turtles-own
[
  Age                         ;   [#, 0-MaxAge]  -              individual age in weeks
  is-female?                  ;   [boolean]      -              individual sex: female or male
  HerdID                      ;   [#]            -              herd membership
  IndCaseFatality             ;   [%]            ind.           individual mortality rate (CaseFatality ^ 2 for adults, sqrt CaseFatality for Juveniles)
  AgeGroup                    ;   [string]       -              age group turtle belongs to: piglet, subadult, adult; differs between sexes
  DemStat                     ;   [string]       dsResident     demographic status: dsResident, dsOffspring, dsDispF, dsDispM
  EpiStat                     ;   [string]       esSusc         epidemiological status: esSusc, esTrans, esLeth, esImm, esImmMat
  shedding?                   ;   [boolean]      0              1 if the individual is infected (EpiStat = esTrans or esLeth
  TickInf                     ;   [#]            -              week of infection
  TickLeth                    ;   [#]            -              week of death after infection (for esLeth)
  Family                      ;   [#]            -              contains state variable [who] of the mother
  ChanceBreeding              ;   [#, 0-1]       -              random number for stochastic effects in monthly reproduction for the given year
  TickBirth                   ;   [#]            -              week of birth
  matAB?                      ;   [boolean]      0              1 if mother had anti bodies (= was immune)
  mother?                     ;   [boolean]      0              1 if female had offspring in current year
  breeder?                    ;   [boolean]      0              1 if female is allowed to breed in current year
]

patches-own
[
  occupied?                   ;   [boolean]      0              1 if a patch contains a herd of wild boars
  Quality                     ;   [#]            [0-9]          habitat quality in terms of breeding capacity = number of females alllowed to breed
  habitat?                    ;   [boolean]      0              true if Quality > 0
  Qual_around                 ;   [#]            [0-9]          summed habitat quality of surrounding cells (to decide if movement is possible)
  Capacity                    ;   [#]            [0-40]         overall capacity = Quality * (mean number of 20 boars per cell / mean Quality) = Quality * (20 / 4.5) = 4.4445 boars per quality unit
  NoInfected                  ;   [#]            0              number of infectious indiviuals of this habitat patch
  InfPress                    ;   [#]            0              infection pressure based on infected indivudals in and around this patch
  infectious?                 ;   [#]            0              1 if patch is infected at the current step
  infected?                   ;   [#]            0              1 if patch has been infected during simulation
  NoRepro                     ;   [#]            0              number of reproductive females (subadult + adult) per cell
  NoDispF                     ;   [#]            0              number of female dispersers per cell
  DispCells                   ;   [set]          -              cells in dispersal distance of the focal cell
]

extensions [ profiler ]






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                                                                                                                                                                                                           ;;
;;  SETUP + GO PROCEDURES                                                                                                                                                                                                                                    ;;
;;                                                                                                                                                                                                                                                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;; SETUP ______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

to Setup

  ca
  reset-ticks
  reset-timer

  set DONE 0

  Init_Parameters
  Init_Landscape
  Init_Population

  View

end



;; GO _________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

to Go

  ;; * ticks = 0 -> week = 1 -> ticks = week - 1 and week = ticks + 1

  set week week + 1

  if ticks mod 52 = 0               [ set week 1 ]

  if ticks = WeekRelease - 1        [ Disease_Release ]  ;; - 1 to be synchronous*

  if ticks >= WeekRelease - 1       [ Disease_Transmission ]   ;; - 1 to be synchronous*

  if week = 17                      [ Movement_Dispersal-Males ]

  if week = 29                      [ Movement_Dispersal-Females ]

  Host_Reproduction

  Host_Death

  Host_Ageing

  if ticks >= WeekRelease - 1       [ Disease_Transition ]

  View

  Report_fragmentation

  if (ticks >= WeekRelease - 1) AND ((not any? turtles with [shedding? = 1]) OR (ticks = SimulatedYears * 52)) [ set DONE 1 ]

  ;; OUTPUT
  if any? turtles with [shedding? = 1] [ set WeekLastInf ticks + 1 ]   ;; + 1 to be synchronous*
  if MaxDist? = 0 AND (any? patches with [pycor = 0 AND infectious? = 1]) [ set MaxDist? ticks
  ]
  if DONE = 1 [ stop ]

  tick

end






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                                                                                                                                                                                                           ;;
;;  INITALIZATION PROCEDURES                                                                                                                                                                                                                                 ;;
;;                                                                                                                                                                                                                                                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;; SET PARAMETERS _____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; Sets all parameters needed for initalisation:
  ;; Current number of habitat patches and maximum number of habitat patches, thresholds for landscape pattern, case fatality, risk of getting lethally infected, week of virus release


to Init_Parameters

  set YearRelease 2
  set SimulatedYears 1000
  set week 0
  set WeekRelease ((YearRelease - 1) * 52) + random 52 + 1   ; YearRelase = 2 -> "during the second year" -> reduce number to year 1; + 1 to be between week 1 and 52
  set QualityRange 9
  set Longevity 11
  set AgeBlur 6
  set ProbFem 0.5
  set CaseFatality 0.5
  set mue 5
  set BetaWithin 0.0208
  set BetaBetween BetaWithin / 10
  set T_trans 1
  set T_anti 15
  set PreNatInf 0.5
  set FertRedInf 0.625
  set MaxAge Longevity * 52
  set InfDist 0
  set MaxDist? 0
  set ReproList [ 0.06    ;; January
                  0.16    ;; February
                  0.39    ;; March
                  0.73    ;; April
                  0.8     ;; May
                  0.88    ;; June
                  0.94    ;; July
                  0.97    ;; August
                  1.0     ;; September
                  0.0     ;; October
                  0.0     ;; November
                  0.0 ]   ;; December

  ;; INPUT DATA
  file-open "RadiusAll.txt"
  set NmaxList []
  while [ not file-at-end? ] [ set NmaxList lput file-read NmaxList ]
  file-close

  ;; OUTPUT
  set InfStructList []
  set InfTransStructList []
  set InfLethStructList []
  set CurrInfPeriodList []
  set LethPeriodList []
  set InfPeriodList []

end



;; LOAD LANDSCAPE _____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; Habitat patches varying in breeding capacity (Quality) are distributed randomly or homogeneous.
  ;; The darker the color is, the higher is the breeding capacity of the home range.


to Init_Landscape

  let x 0
  let y 0

  ;; file import
  file-open (word "nlms/" Filename ".txt")

  while [ not file-at-end? ]
  [
    while [ y < max-pycor + 1 ]
    [
      set x 0
      while [ x < max-pxcor + 1 ]
      [
        let next-value file-read ;; iterate over file entries seperated by white-space
        ask patch x y [ set Quality next-value ]
        set x x + 1
      ]
      set y y + 1
  ] ]
  file-close

  set Habitat patches with [Quality > 0]
  set Matrix patches with [Quality <= 0]
  ask Matrix [ set Quality 0 ]

  let Qual-MW mean [Quality] of patches

  ask Habitat
  [
    set Capacity (Quality * (20 / Qual-MW))   ; Quality * (mean number of 20 boars per cell / mean Quality) = 4.4445 boars per quality unit
    set NoRepro (count turtles-here with [is-female? = 1 AND mother? = 0 AND (AgeGroup = "Adult" OR AgeGroup = "Subdult")])
    set DispCells Habitat in-radius DispDist
    set Qual_around sum [Quality] of neighbors
    set habitat? 1
  ]
  ask Matrix [ set Capacity 0 ]

end



;; CREATE START POPULATION ____________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; One wild boar and its herd is allocated to each habitat patch, group size is initialised as three times breeding capacity (Quality) of each patch.
  ;; Initial age distributions are taken from the results of a 100 years model run conducted by Kramer-Schadt et al. 2009.


to Init_Population

  crt count Habitat
  [
     set Age 0
     set DemStat "dsStart"
     set EpiStat "esSusc"
     set shedding? 0
     set breeder? 0
     set mother? 0
  ]

  while [any? turtles with [DemStat = "dsStart"]]
  [
    ask turtles with [HerdID = 0]
    [
      move-to one-of Habitat with [occupied? = 0]
      set currHerdID currHerdID + 1
      set HerdID currHerdID
      set DemStat "dsResident"
      hatch (4.5 * Quality - 1)  ;; near carrying capacity
      set occupied? 1 ;patches-own
  ] ]

  while [any? turtles with [Age = 0]]
  [
    ask turtles with [Age = 0]
    [
      ifelse random-float 1 <= ProbFem [ set is-female? 1 ] [ set is-female? 0 ]
      let RandomNo random-float 1
      ; generate variation in age of start population
      let Var round random-normal 0 (AgeBlur / 2)
      if Var <= AgeBlur AND Var >= (- AgeBlur)
      [
        ;; intial age distribution based on Kramer-Schadt et al. 2009
        if RandomNo <= 0.38                     [set Age  52 + Var]
        if RandomNo > 0.38 AND RandomNo <= 0.62 [set Age 104 + Var]
        if RandomNo > 0.62 AND RandomNo <= 0.77 [set Age 156 + Var]
        if RandomNo > 0.77 AND RandomNo <= 0.86 [set Age 208 + Var]
        if RandomNo > 0.86 AND RandomNo <= 0.92 [set Age 260 + Var]
        if RandomNo > 0.92 AND RandomNo <= 0.95 [set Age 312 + Var]
        if RandomNo > 0.95 AND RandomNo <= 0.97 [set Age 364 + Var]
        if RandomNo > 0.97 AND RandomNo <= 0.98 [set Age 416 + Var]
        if RandomNo > 0.98 AND RandomNo <= 0.99 [set Age 468 + Var]
        if RandomNo > 0.99 AND RandomNo <= 1    [set Age 520 + Var]
      ]

    ;; age classes males
    if Age <= 21 AND is-female? = 0 [ set AgeGroup "Piglet" ]
    if Age > 21 AND Age <= 104 AND is-female? = 0 [ set AgeGroup "Subadult" ]
    if Age > 104 AND is-female? = 0 [ set AgeGroup "Adult" ]

    ;; age classes females
    if Age <= 34 AND is-female? = 1 [ set AgeGroup "Piglet" ]
    if Age > 34 AND Age <= 52 AND is-female? = 1 [ set AgeGroup "Subadult" ]
    if Age > 52 AND is-female? = 1 [ set AgeGroup "Adult" ]

  ]]

;; OUTPUT
  set PopStructList []
  ask turtles [ set PopStructList fput Age PopStructList ]
  set PopStructList map [? / 52] PopStructList

end






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                                                                                                                                                                                                           ;;
;;  GO PROCEDURES                                                                                                                                                                                                                                            ;;
;;                                                                                                                                                                                                                                                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  DISEASE PROCEDURES                                                                                                                                                                                                                                       ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;; DISEASE RELEASE ____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; One randomly selected wild boar or selected group (cell) in the middle of the top row gets lethally infected in a random week ("WeekRelease") of a given year (based on input "YearRelease").


to Disease_Release

  ask one-of patches with [pycor = max-pycor AND count turtles-here >= 1]
  [
    set InfX pxcor
    set InfY pycor
    set infectious? 1
    set infected? 1
    ask turtles-here [ Disease_Infection ]
  ]

end



;; DISEASE TRANSMISSION _______________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; First, infection pressure is calculated for all habitat cells based on infected indivudals in the herd and the eight herds surrounding the home range cell.
  ;; Afterwards, the model iterates over all individuals and stochastically sets susceptible indivudals to infected if a drawn random number (0-1) is smaller than the infection pressure of the cell.
  ;; The transmission parameter was reversly fitted to the estimated disease spread velocity of approx. 8 km per quarter (Rossi et al. 2010). The transmission probabilities were assigned constant as 0.0208 within a herd and 0.00208 between herds.
  ;; Based on the case fatality rate (CaseFatality), on infection the host is stochastically assigned either as lethally infected (esLeth) or transiently infected (esTrans).
  ;; CaseFatality for lethally infected individuals applies unchanged for yearlings (subadults), is decreased for adults (CaseFatality ^ 2) and increased for piglets (sqrt CaseFatality) to represent age-dependent disease outcomes (Dahle & Liess 1992).
  ;; Lethally infected hosts receive their individual infectious period / survival times (in weeks) which are drawn from an exponential distribution with mean mue (resulting in TickLeth).


to Disease_Transmission

  set NewLeth 0
  set NewTrans 0

  ;; update infection status of cells
  ask Habitat [ set infectious? 0 ]
  ask Habitat with [any? turtles-here with [shedding? = 1]]
  [
    set infectious? 1
    set infected? 1
  ]

  ask Habitat [ set NoInfected count turtles-here with [shedding? = 1] ]

  ask Habitat [ set InfPress (1 - (1 - BetaWithin) ^ NoInfected * (1 - BetaBetween) ^ sum [NoInfected] of neighbors) ]

  ask turtles with [EPiStat = "esSusc" ]
  [
    if random-float 1 < [InfPress] of patch-here [ Disease_Infection ]
  ]

  ;; OUTPUT
  set CurrInfPeriodList []
  ask turtles with [EpiStat = "esLeth"] [ set CurrInfPeriodList fput (TickLeth - ticks) CurrInfPeriodList ]

end



;; DISEASE STATE TRANSITION ___________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; Transient shedders (esTrans) are converted to immune (esImm) after T_trans.
  ;; Individuals protected by maternal antibodies (esImmMat) turn suscpetible (esSusc) after reaching the maximum duration of immunity of 12 weeks.


to Disease_Transition

  ;; transient infecteds
  ask turtles with [EpiStat = "esTrans" AND TickInf + T_trans <= ticks]
  [
    set EpiStat "esImm"
    set shedding? 0
  ]
  ;; maternally protected pigs
  ask turtles with [EpiStat = "esImmMat" AND Age >= T_anti]
  [
    set EpiStat "esSusc"
    set MatAB? 0
  ]

  ask turtles with [matAB? = 1 AND Age > T_anti - 4]   ; 1 month with varying immunity for piglets
  [
    if random-float 1 < 0.5
    [
      set EpiStat "esSusc"
      set matAB? 0
  ] ]

  ;; OUTPUT
  set InfStructList []
  ask turtles with [shedding? = 1]  [ set InfStructList fput Age InfStructList ]
  set InfStructList map [? / 52] InfStructList

  set InfTransStructList []
  ask turtles with [EpiStat = "esTrans"]  [ set InfTransStructList fput Age InfTransStructList ]
  set InfTransStructList map [? / 52] InfTransStructList

  set InfLethStructList []
  ask turtles with [EpiStat = "esLeth"]  [ set InfLethStructList fput Age InfLethStructList ]
  set InfLethStructList map [? / 52] InfLethStructList

end



;; INFECTIOUS PROCESS _________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; turtle procedure

  ;; Individual infection based on case fatality


to Disease_Infection

  if EpiStat = "esSusc" AND AgeGroup = "Piglet"   [ set IndCaseFatality sqrt CaseFatality ]
  if EpiStat = "esSusc" AND AgeGroup = "Subadult" [ set IndCaseFatality CaseFatality ]
  if EpiStat = "esSusc" AND AgeGroup = "Adult"    [ set IndCaseFatality CaseFatality ^ 2 ]

  ifelse random-float 1 <= IndCaseFatality
  [
    set EpiStat "esLeth"
    set shedding? 1
    set TickLeth ticks + 1 + floor (- (mue - 0.5) * ln random-float 1)
    set NewLeth NewLeth + 1
    set LethPeriodList fput (TickLeth - ticks) LethPeriodList  ; list for lethally infected individuals only
    set InfPeriodList fput (TickLeth - ticks) InfPeriodList    ; list for lethally and transient infected individuals
    set infected? 1 ;patches-own
    if distancexy InfX InfY > InfDist [ set InfDist distancexy InfX InfY ]  ;patches-own
  ][ ; else
    set EpiStat "esTrans"
    set shedding? 1
    set TickInf ticks + 1
    set NewTrans NewTrans + 1
    set InfPeriodList fput T_trans InfPeriodList
    set infected? 1  ;patches-own
    if distancexy InfX InfY > InfDist [ set InfDist distancexy InfX InfY ]  ;patches-own
  ]

end




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  HOST PROCEDURES                                                                                                                                                                                                                                          ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;; REPRODUCTION _______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; Adult and subadult females reproduce only once a year depending on their age class.
  ;; In the first week of the year females are checked wether they are able breed (age + no piglets) and if there is at least one reproducible male within their range (> 2 years).
  ;; Litter size is drawn from a truncated normal distribution and reduced to a constant fraction for infected individuals.
  ;; Depending on the disease state of breeding indiviuals its piglets disease states are adjusted.


to Host_Reproduction

  ask Habitat [ set NoRepro count turtles-here with [is-female? = 1 AND mother? = 0 AND Age > 34] ]

  ;; at the beginning of each year check...
  if week = 1
  [
    ;; ... for habitat patches with at least one male to reproduce ...
    ask Habitat
    [
      ;; ... and for all females which can breed
      let num-fem NoRepro
      ;; if number of potential breeders is higher than breeding capacity
      if (num-fem > Quality)
      [
        ;; … sort breeders descending by age …
        let sort-fem sort-on [(- Age)] turtles-here with [is-female? = 1 AND mother? = 0 AND Age > 34]
        set num-fem floor Quality
        if length sort-fem > 0
        [
          while [num-fem > 0]
          [
            ;; … and let num-fem of them (difference between capacity and number of breeders) reproduce
            ask first sort-fem
            [
              set breeder? 1
              set ChanceBreeding random-float 1
            ]
            set sort-fem remove-item 0 sort-fem
            set num-fem (num-fem - 1)
      ] ] ]

      ;; if number of breeders is equal or below than breeding capacity all breeders reproduce
      if (num-fem > 0 AND num-fem <= Quality)
      [
        ask turtles-here with [is-female? = 1 AND AgeGroup = "Adult" AND mother? = 0]
        [
          set breeder? 1
          set ChanceBreeding random-float 1
  ] ] ] ]

  ;; cumulative reproduction probabilities for each month (January to September)
  ifelse week <= 4 [ set ProbRepro item 0 ReproList  Host_prob-Reproduce ][
    ifelse week <= 8 [ set ProbRepro item 1 ReproList  Host_prob-Reproduce ][
      ifelse week <= 13 [ set ProbRepro item 2 ReproList  Host_prob-Reproduce ][
        ifelse week <= 17 [ set ProbRepro item 3 ReproList  Host_prob-Reproduce ][
          ifelse week <= 21 [ set ProbRepro item 4 ReproList  Host_prob-Reproduce ][
            ifelse week <= 26 [ set ProbRepro item 5 ReproList  Host_prob-Reproduce ][
              ifelse week <= 30 [ set ProbRepro item 6 ReproList  Host_prob-Reproduce ][
                ifelse week <= 34 [ set ProbRepro item 7 ReproList  Host_prob-Reproduce ][
                  ifelse week <= 39 [ set ProbRepro item 8 ReproList  Host_prob-Reproduce ][
                    ifelse week <= 43 [ set ProbRepro item 9 ReproList  Host_prob-Reproduce ][
                      ifelse week <= 47 [ set ProbRepro item 10 ReproList  Host_prob-Reproduce ][
                        if week <= 52 [ set ProbRepro item 11 ReproList  Host_prob-Reproduce
  ] ] ] ] ] ] ] ] ] ] ] ]

end

to Host_prob-Reproduce
  ;; observer procedure

    ask turtles with [breeder? = 1]
    [
      ;; stochastic effect in probabilities of reproduction
      if ChanceBreeding < ProbRepro
      [
        let RandomNo random-float 1

        ;; stochastic effect in number of piglets
        if (RandomNo > 0.01305859 AND RandomNo <= 0.082208805)
        [
          ;; if the breeder is ifencted number of piglets is reduced
          ifelse shedding? = 1
          [
            hatch floor (1 * FertRedInf)   ;; round to integer
            [
              set DemStat "dsOffspring"
              set family who
          ] ]
          [
          ;; if the breeder is not infected
            hatch 1
            [
              set DemStat "dsOffspring"
              set family who
        ] ] ]

        if (RandomNo > 0.082208805 AND RandomNo <= 0.24510931)
        [
          ;; if the breeder is ifencted number of piglets is reduced
          ifelse shedding? = 1
          [
            hatch floor (2 * FertRedInf)
            [
              set DemStat "dsOffspring"
              set family who
            ]
          ][
            ;; if the breeder is not infected
            hatch 2
            [
              set DemStat "dsOffspring"
              set family who
        ] ] ]

        if (RandomNo > 0.24510931 AND RandomNo <= 0.495050085)
        [
          ;; if the breeder is ifencted number of piglets is reduced
          ifelse shedding? = 1
          [
            hatch floor (3 * FertRedInf)
            [
              set DemStat "dsOffspring"
              set family who
            ]
          ][
            ;; if the breeder is not infected
            hatch 3
            [
              set DemStat "dsOffspring"
              set family who
        ] ] ]

        if (RandomNo > 0.495050085 AND RandomNo <= 0.744990859)
        [
          ;; if the breeder is ifencted number of piglets is reduced
          ifelse shedding? = 1
          [
            hatch floor (4 * FertRedInf)
            [
              set DemStat "dsOffspring"
              set family who
          ] ]
          [
            ;; if the breeder is not infected
            hatch 4
            [
              set DemStat "dsOffspring"
              set family who
        ] ] ]

        if (RandomNo > 0.744990859 AND RandomNo <= 0.907891364)
        [
          ;; if the breeder is ifencted number of piglets is reduced
          ifelse shedding? = 1
          [
            hatch floor (5 * FertRedInf)
            [
              set DemStat "dsOffspring"
              set family who
            ]
          ][
            ;; if the breeder is not infected
            hatch 5
            [
              set DemStat "dsOffspring"
              set family who
        ] ] ]

        if (RandomNo > 0.907891364 AND RandomNo <= 0.977041579)
        [
          ;; if the breeder is ifencted number of piglets is reduced
          ifelse shedding? = 1
          [
            hatch floor (6 * FertRedInf)
            [
              set DemStat "dsOffspring"
              set family who
            ]
          ][
            ;; if the breeder is not infected
            hatch 6
            [
              set DemStat "dsOffspring"
              set family who
        ] ] ]

        if (RandomNo > 0.977041579 AND RandomNo <= 0.996144138)
        [
          ;; if the breeder is ifencted number of piglets is reduced
          ifelse shedding? = 1
          [
            hatch floor (7 * FertRedInf)
            [
              set DemStat "dsOffspring"
              set family who
          ] ]
          [
            ;; if the breeder is not infected
            hatch 7
            [
              set DemStat "dsOffspring"
              set family who
        ] ] ]

        if (RandomNo > 0.996144138 AND RandomNo <= 0.999575749)
        [
          ;; if the breeder is ifencted number of piglets is reduced
          ifelse shedding? = 1
          [
            hatch floor (8 * FertRedInf)
            [
              set DemStat "dsOffspring"
              set family who
          ] ]
          [
            ;; if the breeder is not infected
            hatch 8
            [
              set DemStat "dsOffspring"
              set family who
        ] ] ]

        if (RandomNo > 0.999575749 AND RandomNo <= 0.99997575)
        [
          ;; if the breeder is ifencted number of piglets is reduced
          ifelse shedding? = 1
          [
            hatch floor (9 * FertRedInf)
            [
              set DemStat "dsOffspring"
              set family who
          ] ]
          [
            ;; if the breeder is not infected
            hatch 9
            [
              set DemStat "dsOffspring"
              set family who
        ] ] ]

        if (RandomNo > 0.99997575 AND RandomNo <= 1.0)
        [
          ;; if the breeder is ifencted number of piglets is reduced
          ifelse shedding? = 1
          [
            hatch floor (10 * FertRedInf)
            [
              set DemStat "dsOffspring"
              set family who
          ] ]
          [
            ;; if the breeder is not infected
            hatch 10
            [
              set DemStat "dsOffspring"
              set family who
        ] ] ]

        set TickBirth ticks
        set breeder? 0
        set mother? 1
  ] ]

  ask turtles with [DemStat = "dsOffspring"]
  [
    ifelse random-float 1 <= ProbFem [ set is-female? 1 ] [ set is-female? 0 ]
    set Age 0
    set AgeGroup "Piglet"
    set breeder? 0

    ;; Susceptibles
    if [EpiStat] of turtle family = "esSusc"
    [
      set EpiStat "esSusc"
      set shedding? 0
      set DemStat "dsResident"
    ]

    ;; Immunes
    if [EpiStat] of turtle family = "esImm"
    [
      set EpiStat "esImmMat"
      set shedding? 0
      set DemStat "dsResident"
      set matAB? 1
    ]

    ;; Lethaly Infecteds
    if [EpiStat] of turtle family = "esLeth"
    [
      set EpiStat "esLeth"
      set shedding? 1
      set DemStat "dsResident"
      set TickLeth ticks + 1 + floor (- (mue - 0.5) * ln random-float 1)
      set NewLeth NewLeth + 1
      set LethPeriodList fput (TickLeth - ticks) LethPeriodList
      set InfPeriodList fput (TickLeth - ticks) InfPeriodList
    ]

    ;; Transient Infecteds
    if [EpiStat] of turtle family = "esTrans"
    [
      ifelse random-float 1 < PreNatInf
      [
        set EpiStat "esLeth"
        set shedding? 1
        set DemStat "dsResident"
        set TickLeth ticks + 1 + floor (- (mue - 0.5) * ln random-float 1)
        set NewLeth NewLeth + 1
        set LethPeriodList fput (TickLeth - ticks) LethPeriodList
        set InfPeriodList fput (TickLeth - ticks) InfPeriodList
      ]
      [
        set EpiStat "esSusc"
        set shedding? 0
        set DemStat "dsResident"
        set TickInf ticks
  ] ] ]

end



;; DEATH ______________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; Iterating over the entire population, each individual either stochastically dies due to age class dependent mortality rates by drawing a random number,
  ;;    due to reaching a certain maximum age or due to a lethal infection after a certain ifnection time span (TickLeth).
  ;; Stochastic baseline mortality is age-dependent and are drawn from a Gaussian distribution which was adjusted to annual survival estimates found in the literature and is cutted symmetrically around the mean.
  ;; The stochastic effect resembles 'good' versus 'bad' years for wild boars (e.g. environmental noise).
  ;; Each time step the adjusted age-dependent mortality (MortAge) is applied to all individual.


to Host_Death

  if week = 1
  [
    let Noise 999
    while [ Noise > 1 OR Noise < -1 ]
    [
      set Noise random-normal 0 (1 / 2.576)   ; 99% between -1 and 1
      ;; survival probability for subadults and adults
      set SurvivalProb 0.65 + Noise * (0.65 - 0.4)   ; mean + Noise * (mean - min) to scale for age group
      ;; survival probability for piglets
      set SurvivalProbPig 0.5 + Noise * (0.5 - 0.1)
    ]
  ]

  ask turtles
  [
    ifelse Age > 34
    [
      if random-float 1 < (1 - (SurvivalProb ^ (1 / 52))) [ die ]
    ][  ;; mean survival probability for yearlings (Gaillard et al. 1987) = 0.65 = mean survival probability for adults (Focardi et al. 1996)
      if random-float 1 < (1 - (SurvivalProbPig ^ (1 / 52))) [ die ]
  ] ]

  ;; mortality due to a lethal infection
  ask turtles with [EpiStat = "esLeth"]
  [
    if TickLeth = ticks [ die ]
  ]

  ;; mortality due to reaching maximum age (more than 10 years)
  ask turtles with [Age > MaxAge] [ die ]

  ask Habitat with [not any? turtles-here AND occupied? = 1] [ set occupied? 0 ]

  ;;OUTPUT
  set PopStructList []
  ask turtles [ set PopStructList fput Age PopStructList ]
  set PopStructList map [? / 52] PopStructList

end



;; AGE ________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; For each individual age is incremented one week and disease state transitions are performed.
  ;; Females who gave birth more than 33 weeks before can reproduce again.

to Host_Ageing

  ask turtles
  [
    set Age Age + 1

    ifelse Age = 35
    [
      set AgeGroup "Subadult"
      if is-female? = 0 [ set DemStat "dsDispM" ]
    ]
    [
      ifelse Age = 53 AND is-female? = 1
      [ set AgeGroup "Adult" ]
      [
        if Age = 105 AND is-female? = 0
        [
          set AgeGroup "Adult"
          if Qual_around > 0 [ set DemStat "dsResident" ]  ;patches-own
    ] ] ]

    if mother? = 1 AND (TickBirth + 34) = ticks [ set mother? 0 ]
  ]

end




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  HOST MOVEMENT PROCEDURES                                                                                                                                                                                                                                 ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;; HERD SPLITTING _____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; During summer, subadult females either remain with their mothers or establish new territories nearby.
  ;; Herd splitting is performed in one week per year only where female subadults (8 months = 34 weeks < Age <= 2 years = 104 weeks) without offspring are able to move between herds.
  ;; All herds to split are extracted matching the conditions of containing a number of females exceeding the cells carrying capacity and containing at least 2 dispersers.
  ;; For each splittable herd an habitat cell without female groups and a carrying capacity above 0 within the dispersal distance of 6 km (3 cells) is selected randomly. All dispersers from the source cell find a new herd on that cell.


to Movement_Dispersal-Females

  let OverCap Habitat with [count turtles-here with [is-female? = 1 AND Age > 34] - Quality >= 2 AND Qual_around > 0]

  ask OverCap
  [
    ;; count females above breeding capacity (Quality)
    let num-overfem ceiling (count turtles-here with [is-female? = 1 AND Age > 34] - Quality)
    ;; ... assign those females a disperser-status...
    let num-disp count turtles-here with [AgeGroup = "Subadult" AND is-female? = 1 AND breeder? = 0]
    ;; ... which are over breeding capacity ...
    if num-disp > num-overfem [ set num-disp num-overfem ]
    ;; ... and if there are at least 2 dispersers ...
    if num-disp >= 2
    [ ask n-of num-disp turtles-here with [AgeGroup = "Subadult" AND is-female? = 1 AND breeder? = 0] [ set DemStat "dsDispF" ] ]
    set NoDispF 0
  ]

  ask OverCap with [any? turtles-here with [DemStat = "dsDispF"]]
  [

    let ListDispF []
    ask turtles-here with [DemStat = "dsDispF"] [ set ListDispF fput self ListDispF ]
    set NoDispF length ListDispF   ; number of dispersing females of this cell
    let tempDispCells DispCells with [not any? turtles-here with [is-female? = 1]]
    ifelse count tempDispCells > 0
    [
      let newHabitat one-of tempDispCells
      while [not empty? ListDispF]
      [
        ask first ListDispF
        [
          move-to newhabitat
          let cHerdID HerdID
          set HerdID currHerdID
          set DemStat "esResident"
          if shedding? = 1
          [
            set infected? 1 ;patches-own
            if distancexy InfX InfY > InfDist [ set InfDist distancexy InfX InfY ] ;patches-own
          ]
          set ListDispF remove-item 0 ListDispF
        ]
      ]
      ask newhabitat [ set occupied? 1 ]
      set ListDispF []
    ]
    ;; else = if there are no unoccupied habitat cells within their dispersal distance the dispersers emmigrate/die
    [
      ask turtles-here with [DemStat = "dsDispF"]
      [
        ;; version used in Fernandez et al. 2006: dispersers stay in their maternal cell if no empty habitat cell is available
        set DemStat "esResident"
    ] ]
    set currHerdID currHerdID + 1
  ]

end



;;; MALE DISPERSAL _____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
  ;; observer procedure

  ;; Males perform their natal dispersal in an age of 8 months (35 weeks) upwards. All dispersing males form a group of at least 3 subadults and move to a cells below capacity with the highest habitat quality in dispersal distance.

to Movement_Dispersal-Males

  let NatalMales 0

  ask Habitat with [any? turtles-here with [DemStat = "dsDispM"] AND Qual_around > 0]
  [
    set NatalMales turtles-here with [DemStat = "dsDispM"]
    if count NatalMales > 2
    [
      let NewHome DispCells with [count turtles-here < Capacity]
      if count NewHome > 0
      [
        set NewHome max-one-of NewHome [Quality]
        ask NatalMales
        [
          move-to NewHome
          set DemStat "dsResident"
  ] ] ] ]

end





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; HELP FUNCTIONS                                                                                                                                                                                                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;; FRAGMENTAITON INDEX by Mitchell & Powell 2004 ______________________________________________________________________________________________________________________________________________________________________________________________________________

to Report_fragmentation

  set F_infected 0
  set F_infectious 0

  let Nmax_infected count patches with [infected? = 1]
  let Nmax_infectious count patches with [any? turtles-here with [shedding? = 1]]

  set Nmax_infected item Nmax_infected NmaxList
  set Nmax_infectious item Nmax_infectious NmaxList

  if Nmax_infected > 0 [ set F_infected (1 - ((sum ([count neighbors with [infected? = 1]] of patches with [infected? = 1])) / Nmax_infected )) ]
  if Nmax_infectious > 0 [ set F_infectious (1 - ((sum ([count neighbors with [any? turtles-here with [shedding? = 1]]] of patches with [any? turtles-here with [shedding? = 1]])) / Nmax_infectious )) ]

end



;; UPDATE VIEW ________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________

to View

  ask turtles with [EpiStat = "esSusc"] [ set color black ]
  ask turtles with [EpiStat = "esImm" OR EpiStat = "esImmMat"] [ set color green ]
  ask turtles with [EpiStat = "esLeth"]
  [
    set color red
    set size 0.75
  ]
  ask turtles with [EpiStat = "esTrans"]
  [
    set color yellow
    set size 0.75
  ]
  ask patches [ set pcolor scale-color grey Quality 12 0 ]
  ask Habitat with [any? turtles-here with [shedding? = 1]] [ set pcolor scale-color orange Quality 12 0 ]
  if ticks >= WeekRelease - 1 [ ask patch InfX InfY [ set pcolor yellow ] ]

end







;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                                                                                                                                                                                                                           ;;
;; END                                                                                                                                                                                                                                                       ;;
;;                                                                                                                                                                                                                                                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@#$#@#$#@
GRAPHICS-WINDOW
223
10
650
876
-1
-1
16.7
1
10
1
1
1
0
0
0
1
0
24
0
49
1
1
1
ticks
30.0

BUTTON
7
10
62
43
Setup
Setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
64
10
138
43
Run Week
Go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
140
10
214
43
Run Year
Go\nif ticks mod 52 = 0 [stop]
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
123
45
214
78
Run Until End
if ticks <= (SimulatedYears * 52) [ Go ]
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
7
45
120
78
Run Until Infection
while [ ticks + 1 < (WeekRelease) ] [ Go ]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
10
158
60
203
Year
ceiling (ticks / 52)
0
1
11

TEXTBOX
188
548
338
566
NIL
11
0.0
1

MONITOR
11
252
85
293
New Transient
NewTrans
17
1
10

MONITOR
11
209
85
250
New Lethal
NewLeth
17
1
10

MONITOR
65
158
115
203
Week
week
17
1
11

INPUTBOX
7
82
119
142
Filename
frag_high_4
1
0
String

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 5.3
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="Experiment" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>count turtles with [is-female? = 0]</metric>
    <metric>count turtles with [is-female? = 1]</metric>
    <metric>count turtles with [EpiStat = "esTrans"]</metric>
    <metric>count turtles with [EpiStat = "esLeth"]</metric>
    <metric>count turtles with [EpiStat = "esImmu"]</metric>
    <metric>count turtles with [Age &lt;= 34]</metric>
    <metric>count turtles with [Age &lt;= 34 AND EpiStat = "esTrans"]</metric>
    <metric>count turtles with [Age &lt;= 34 AND EpiStat = "esLeth"]</metric>
    <metric>count turtles with [Age &lt;= 34 AND EpiStat = "esImmu"]</metric>
    <metric>count turtles with [Age &gt; 104]</metric>
    <metric>count turtles with [Age &gt; 104 AND EpiStat = "esTrans"]</metric>
    <metric>count turtles with [Age &gt; 104 AND EpiStat = "esLeth"]</metric>
    <metric>count turtles with [Age &gt; 104 AND EpiStat = "esImmu"]</metric>
    <metric>NewTrans</metric>
    <metric>NewLeth</metric>
    <metric>MaxLethPeriod</metric>
    <metric>MeanLethPeriod</metric>
    <metric>count patches with [infected? = 1]</metric>
    <metric>count patches with [infectious? = 1]</metric>
    <metric>F_infectious</metric>
    <metric>InfDist</metric>
    <metric>MaxDist?</metric>
    <metric>WeekRelease</metric>
    <metric>WeekLast</metric>
    <enumeratedValueSet variable="Filename">
      <value value="&quot;mpd-0.0&quot;"/>
      <value value="&quot;mpd-0.5&quot;"/>
      <value value="&quot;mpd-1.0&quot;"/>
      <value value="&quot;mpd-1.5&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seed?">
      <value value="&quot;BS&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Landscape">
      <value value="&quot;File-Import&quot;"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
0
@#$#@#$#@
