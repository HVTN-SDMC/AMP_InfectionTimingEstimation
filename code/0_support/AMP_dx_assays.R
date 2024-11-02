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
 
  # ITRIs for USE IN ESTIMATING DDI like IDT does. Note that you'll want to add ~4.2 to all of these, because this is the approximate delay between 1 copy / ml and the beginning of the ITRI, according to Facente et al; but also because this is about the average diff between the IDT estimates and our tsic2 estimates. This is the diagnostic delay for the aptima assay according to Facente/Grebe. It appears that when the assay is not aptima, it is relative to aptima, but aptima is relative to 1 copy / ml.
  all_assays[[tolower('iSCAv2')]] <- list(
    class = 'RNA',
    full_assayname = 'iSCA v2.0',
    short_assayname = 'iSCAv2',
    form = 'weib3',
    source = "tosiano_2019",
    fun = 'weib3_assay_dynamics',
    params = list(location = 0, shape = 1, scale = 0.187)
  )
 all_assays[['taqman']] <- list(
   class = 'RNA',
   full_assayname = 'Roche Taqman v2.0',
   short_assayname = 'Taqman',
   form = 'weib3',
   source = "manufacturer's details and IDT",
   fun = 'weib3_assay_dynamics',
   params = list(location =0, shape =2.371, scale =2.234)
  )
 all_assays[['abbott_real_time']] <- list(
   class = 'RNA',
   full_assayname = 'Abbott Real Time HIV01 v1.0 m2000sp/m2000rt',
   short_assayname = 'Abbott Real Time',
   form = 'weib3',
   source = "manufacturer's details and IDT",
   fun = 'weib3_assay_dynamics',
   params = list(location = 5.1252, shape = 1.350, scale = 9.660)
 )
  all_assays[['architect']] <- list(
    class = 'Ag/Ab',
    full_assayname = 'Abbott Architect HIV Ag/Ab Combo',
    short_assayname = 'Architect',
    form = 'weib3',
    source = 'IDT',
    fun = 'weib3_assay_dynamics',
    params = list(location = 6.48, shape =1, scale = 1/0.553)
  )
  all_assays[['gs_combo']] <- list(
    class = 'Ag/Ab',
    full_assayname = 'BioRad GS HIV Combo Ag/Ab EIA',
    short_assayname = 'GS_Combo',
    form = 'weib3',
    source = 'IDT',
    fun = 'weib3_assay_dynamics',
    params = list(location = 6.06, shape = 1, scale =1/0.592)
  )
  all_assays[['geenius_indet']] <- list(
    class = 'IgG_Rapid',
    full_assayname = 'BioRad Geenius Indeterminate',
    short_assayname = 'Geenius_Indet',
    form = 'weib3',
    source = 'IDT',
    fun = 'weib3_assay_dynamics',
    params = list(location = 14.88, shape = 1, scale = 1/0.241 )
  )
  all_assays[['geenius_fr']] <- list(
    class = 'IgG_Supp',
    full_assayname = 'BioRad Geenius Fully Reactive',
    short_assayname = 'Geenius_FR',
    form = 'weib3',
    source = 'IDT',
    fun = 'weib3_assay_dynamics',
    params = list(location = 17.28, shape = 1, scale = 1/0.208 )
  )

all_assay_dynamics <- all_assays;
