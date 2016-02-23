#########################################################################################
  # Contains all the constants and predefined variables need for all 
  #  verb_math scripts
check_config = {

# Limits
"event_cut":5,          # Number of events to process, -1 for all events
"event_out":10,             # How often to output event number
"SP_Limit":12,              # Number of PE needed to recognize a space point
"T_Data":11,                # Number of reconstructed track points with data
"P_Cut":0.05,               # Lower limit for reconstructed track P-values


# Lists
"Tker_List":["upstream", "downstream"],          # Labels for tracking detectors
"TOF_List":["tof0", "tof1", "tof2"],             # Labels for TOF detectors
"Station_Location_List":{0:{5:[13963.04,13962.39,13961.73], \
                            4:[14313.11,14312.16,14311.81], \
                            3:[14613.04,14612.39,14611.74], \
                            2:[14862.04,14861.39,14860.74], \
                            1:[15062.94,15062.28,15061.63]},\
                         1:{1:[18848.41,18849.07,18849.72], \
                            2:[19048.38,19049.03,19049.68], \
                            3:[19298.23,19298.88,19299.54], \
                            4:[19598.23,19598.88,19599.53], \
                            5:[19948.14,19948.79,19949.45]} },


"ignore_list":{"ignore_SP_to_Virt": False,         \
              "ignore_SP_Fill_ROOT": False,        \
              "ignore_Virt_Fill_ROOT": False,      \
              "ignore_TOF_Timing_Info": True,      \
              "ignore_Generate_Virtual_Map": True, \
              "ignore_Station_Alignment": False},


# Counters
"counter":{"upstream_residual_count":0, \
           "downstream_residual_count":0, \
           "no_upstream_tofs":0, \
           "no_downstream_tof":0, \
           "total_events":0},


# Maus data
#"data_identifier":"7417",
#"data_directory":"/vols/fets2/heidt/offline/data/7417/",
"data_identifier":"maus_output_1.4.0_736_preprod_RyRx_correction.root",
"data_directory":"/vols/fets2/heidt/offline/simulation/Test/alpha_rotation_flip/",
"output_file":"test.root"
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

"Number_of_Points": 5,

}