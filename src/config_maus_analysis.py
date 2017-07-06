import os
#########################################################################################
  # Contains all the constants and predefined variables need for all
  #  verb_math scripts
check_config = {
# Limits and Cuts
"event_cut":5000,             # Maximum event to process to, negative for max
"min_event":0,              # Minimum event to process, testing purposes
"event_out":1000,           # How often to output event number
"SP_Limit":12,              # Number of PE needed to recognize a space point


#"upstream_Tmin":26.65,      # Minimum timing for TOF0 to TOF1
#"upstream_Tmax":29.83,      # Maximum timing for TOF0 to TOF1
"upstream_Tmin":20.0,      # Minimum timing for TOF0 to TOF1
"upstream_Tmax":60.0,      # Maximum timing for TOF0 to TOF1
#"downstream_Tmin":26.04,    # Minimum timing for TOF1 to TOF2
#"downstream_Tmax":34.84,    # Maximum timing for TOF1 to TOF2
"downstream_Tmin":20.0,    # Minimum timing for TOF1 to TOF2
"downstream_Tmax":60.0,    # Maximum timing for TOF1 to TOF2

# Analysis functions, turns functions On(False)/Off(True)
"ignore_SP_to_Virt":             True,
"ignore_SP_Fill_ROOT":           False,
"ignore_Virt_Fill_ROOT":         True,
"ignore_TOF_Timing_Info":        False,
"ignore_Generate_Virtual_Map":   True,
"ignore_Station_Alignment":      True,
"ignore_TOF_Tkr_Res":            True,
"ignore_MC_Study":               True,
"ignore_Check_Channel_Overlap":  True,
"ignore_emittance":              False,
"ignore_Cut_Collection":         False,


# Counters
"counter":{"upstream_residual_count":0, \
           "downstream_residual_count":0, \
           "no_upstream_tofs":0, \
           "no_downstream_tof":0, \
           "total_events":0},


# Maus data
#"file_name":["08873_recon.root", "08874_recon.root", "08875_recon.root", "08877_recon.root", \
#             "08878_recon.root", "08879_recon.root", "08880_recon.root", "08881_recon.root", \
#             "08882_recon.root"],
#"data_directory": '/media/chris/Research/data/6-240/',
#"data_directory": '/home/chris/work/simulation/simple_beam/',
"data_directory": '/media/chris/Research/data/MAUS-2.7.0/6-240/',
"file_name": ["08873_recon.root", "08874_recon.root", "08875_recon.root"], #, \
#              "08877_recon.root", "08878_recon.root", "08879_recon.root", \
#              "08880_recon.root", "08881_recon.root", "08882_recon.root"],
#"data_directory":"/vols/fets2/heidt/offline/simulation/Test/alpha_rotation_flip/",
"log_file":os.environ['ANALYSIS_DIR']+"/tmp/log.log",
}

fill_config = {
"Virtual_to_Tracker_Map":{ 31:[0,5,2], 32:[0,5,1], 33:[0,5,0], \
                           35:[0,4,2], 36:[0,4,1], 37:[0,4,0], \
                           39:[0,3,2], 40:[0,3,1], 41:[0,3,0], \
                           42:[0,2,2], 43:[0,2,1], 44:[0,2,0], \
                           46:[0,1,2], 47:[0,1,1], 48:[0,1,0], \
                           56:[1,1,0], 57:[1,1,1], 58:[1,1,2], \
                           60:[1,2,0], 61:[1,2,1], 62:[1,2,2], \
                           63:[1,3,0], 64:[1,3,1], 65:[1,3,2], \
                           67:[1,4,0], 68:[1,4,1], 69:[1,4,2], \
                           70:[1,5,0], 71:[1,5,1], 72:[1,5,2] }

}

analysis_config = {
"tk_data_req":14,          # Number of reconstructed track points with data
"min_data":9,              # Bad data criteria for events with two tracks.
"p_cut":0.05,              # Lower limit for reconstructed track P-values
"min_trip":3,
"req_trip":5,
}

out_config = {
"output_dir":os.environ['ANALYSIS_DIR']+"/output/",
"output_file":"_output_140_Diff0_lattice1_5_LiH.root",
"log_file":os.environ['ANALYSIS_DIR']+"/tmp/log.log",
"verbose":False,
}

emittance_config = {
"cut_list":{"Diff_Scrap":True}
}

cut_config = {
"c_mu": 51240,
"c_pi": 60150,
"mu_mu": 27.24,
"mu_pi": 28.56,
"sigma_mu": 0.1547,
"sigma_pi": 0.2045
}