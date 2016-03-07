import os
#########################################################################################
  # Contains all the constants and predefined variables need for all 
  #  verb_math scripts
check_config = {
# Limits and Cuts
"event_cut":-5,             # Number of events to process, -1 for all events
"event_out":1000,           # How often to output event number
"SP_Limit":12,              # Number of PE needed to recognize a space point
"T_Data":11,                # Number of reconstructed track points with data
"P_Cut":0.05,               # Lower limit for reconstructed track P-values
"upstream_Tmin":26.65,      # Minimum timing for TOF0 to TOF1 
"upstream_Tmax":29.83,      # Maximum timing for TOF0 to TOF1
"downstream_Tmin":28.04,    # Minimum timing for TOF1 to TOF2 
"downstream_Tmax":34.84,    # Maximum timing for TOF1 to TOF2


# Lists
"Tker_List":["upstream", "downstream"],          # Labels for tracking detectors
"TOF_List":["tof0", "tof1", "tof2"],             # Labels for TOF detectors


# Analysis functions, turns functions On(False)/Off(True)
"ignore_SP_to_Virt":             True,
"ignore_SP_Fill_ROOT":           False,
"ignore_Virt_Fill_ROOT":         True,
"ignore_TOF_Timing_Info":        False,
"ignore_Generate_Virtual_Map":   True,
"ignore_Station_Alignment":      True,
"ignore_TOF_to_TOF_Tkr_Res":     False,


# Counters
"counter":{"upstream_residual_count":0, \
           "downstream_residual_count":0, \
           "no_upstream_tofs":0, \
           "no_downstream_tof":0, \
           "total_events":0},


# Maus data
"file_name":"maus_output.root",
"data_directory":"/vols/fets2/heidt/offline/analysis/7417/",
#"file_name":"maus_output_geo_123_geo_inverse.root",
#"data_directory":"/vols/fets2/heidt/offline/simulation/Test/alpha_rotation_flip/",
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
"TtT_trip_req":5,           # Required triplets PR for TOF_to_TOF_Tkr_Res
"Number_of_Points": 5,
}

out_config = {
"output_dir":os.environ['ANALYSIS_DIR']+"/output/",
"output_file":"_output_geo_123_geo_7417.root"
}
