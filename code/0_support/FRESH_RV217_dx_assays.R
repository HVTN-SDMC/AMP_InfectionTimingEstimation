   all_assays <- list()
 
   # BASICs for UNIT TESTING
   all_assays[[tolower('step_unit_testing')]] <- list(
     class = 'unit testing',
     full_assayname = 'Unit Testing',
     short_assayname = 'Unit_Testing',
     form = 'step',
     source = 'made-up',
     fun = 'step_assay_dynamics',
     params = list(diagnostic_delay = 10)
   )
   all_assays[[tolower('linear_unit_testing')]] <- list(
     class = 'unit testing',
     full_assayname = 'Unit Testing',
     short_assayname = 'Unit_Testing',
     form = 'linear_abs_spread',
     source = 'made-up',
     fun = 'linear_assay_dynamics',
     params = list(diagnostic_delay = 10, abs_spread = 10)
   )
   all_assays[[tolower('weib3_unit_testing')]] <- list(
     class = 'unit testing',
     full_assayname = 'Unit Testing',
     short_assayname = 'Unit_Testing',
     form = 'weib3',
     source = 'aptima_in_delaney_2017',
     fun = 'weib3_assay_dynamics',
     params = list(location = 4.8, shape = 1.35, scale = 9)
   )
 

  # ITRIs for USE IN ESTIMATING DDI like IDT does. Note that you'll want to add ~4.2 to all of these (except Aptima), because this is the approximate delay between 1 copy / ml and the beginning of the ITRI, according to Facente et al; but also because this is about the average diff between the IDT estimates and our tsic2 estimates. This is the diagnostic delay for the aptima assay according to Facente/Grebe. It appears that when the assay is not aptima, it is relative to aptima, but aptima is relative to 1 copy / ml.
  all_assays[['aptima']] <- list(
    class = 'RNA',
    full_assayname = 'Aptima HIV-1 RNA Qualitative Assay',
    short_assayname = 'Aptima RNA',
    form = 'weib3',
    source = 'IDT',
    fun = 'weib3_assay_dynamics',
    params = list(location = 2.53, shape = 1, scale =1/1.417)
  )
  all_assays[['nuclisens']] <- list(
    class = 'RNA',
    full_assayname = 'Organon Teknika EQL Nuclisens System',
    short_assayname = 'Nuclisens',
    form = 'weib3',
    source = "manufacturer's details",
    fun = 'weib3_assay_dynamics',
    params = list(location = 0, shape = 3.93, scale = 5.02 )
  )
 all_assays[['determine']] <- list(
   class = 'Ag/Ab_Rapid',
   full_assayname = 'Determine HIV-1/2 Ag/Ab Combo',
   short_assayname = 'Determine',
   form = 'weib3',
   source = 'IDT',
   fun = 'weib3_assay_dynamics',
   params = list(location = 7.56, shape = 1, scale =1/0.474)
 )
 all_assays[['unigold']] <- list(
   class = 'Ag/Ab_Rapid',
   full_assayname = 'Trinity Biotech Unigold Rapid HIV Test',
   short_assayname = 'Unigold',
   form = 'weib3',
   source = 'IDT',
   fun = 'weib3_assay_dynamics',
   params = list(location = 15.06, shape = 1, scale =1/0.238)
 )
 all_assays[['cobas_combi']] <- list(
   class = 'Ag/Ab_Lab',
   full_assayname = ' Roche cobas Core Anti-HIV-1/HIV-2 DAGS',
   short_assayname = 'Cobas',
   form = 'weib3',
   source = 'IDT',
   fun = 'weib3_assay_dynamics',
   params = list(location = 10.02, shape = 1, scale =1/0.358)
 )
 all_assays[['centaur']] <- list(
   class = 'Ag/Ab_Lab',
   full_assayname = 'Siemens Advia Centaur HIV 1/O/2 Enhanced Assay',
   short_assayname = 'Centaur',
   form = 'weib3',
   source = 'IDT',
   fun = 'weib3_assay_dynamics',
   params = list(location = 9.12, shape = 1, scale =1/0.393)
 )
 all_assays[['wb_fr']] <- list(
   class = 'Western',
   full_assayname = 'BioRad GS HIV-1 Western blot Fully Reactive',
   short_assayname = 'Western_fr',
   form = 'weib3',
   source = 'IDT',
   fun = 'weib3_assay_dynamics',
   params = list(location = 17.76, shape = 1, scale =1/0.202)
 )
 all_assays[['wb_indet']] <- list(
   class = 'Western',
   full_assayname = 'BioRad GS HIV-1 Western blot Indeterminate',
   short_assayname = 'Western_indet',
   form = 'weib3',
   source = 'IDT',
   fun = 'weib3_assay_dynamics',
   params = list(location = 8.88, shape = 1, scale =1/0.404)
 )
 all_assays[['gs_eia']] <- list(
   class = 'IgG/IgM_Lab',
   full_assayname = 'GS HIV-1/HIV-2 PLUS O EIA',
   short_assayname = 'GS PLUS O EIA',
   form = 'weib3',
   source = 'IDT',
   fun = 'weib3_assay_dynamics',
   params = list(location = 10.86, shape = 1, scale = 1/0.33)
 )

all_assay_dynamics <- all_assays;
