import numpy as np
import labber_script_extension.script_object
# export LABBERPATH = /usr/share/Labber/Script/

globals_path = r'/home/antti/CloudStation/codes/labber-script-extension/examples/'
globals_file_name = 'globals_test'
template_path = '/home/antti/CloudStation/codes/labber-script-extension/examples/test_master.labber'
output_directory_root = '/home/antti/CloudStation/codes/labber-script-extension/examples/'
output_file_name = 'test_3'
tags = []
comment = 'no comments'
local_step_channels = {}
looped_variables = {}
local_parameters = {}
local_instrument_values = {}
log_channels = []

labber_script_extension.script_object.updateAndPerformMeasurement(template_path, output_directory_root, output_file_name, tags, comment, globals_path,
                                                                  globals_file_name, local_step_channels, looped_variables, local_parameters, local_instrument_values, log_channels=log_channels)
