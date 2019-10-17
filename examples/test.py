import numpy as np
#import labber_script_extension.script_object
import labber_script_extension.script_object
# export LABBERPATH = /usr/share/Labber/Script/

globals_path = r'./'
globals_file_name = 'globals_test'
template_path = './test_master2.labber'
output_directory_root = os.getcwd()
output_file_name = 'test_9'
tags = []
comment = 'no comments'
local_step_channels = {}
looped_variables = {}
local_parameters = {}
local_instrument_values = {}
log_channels = ['Manual - Value 3']

so = labber_script_extension.script_object.ScriptObject(
    template_path, os.path.join(output_directory_root, '2019', output_file_name))

# so.printStepChannels()
# so.getStepChannels()
print(so.getChannelNames())
print(so.getLogChannels())
so.setLogChannels('Manual - Value 3')
print(so.getLogChannels())
so.performMeasurement()
# labber_script_extension.script_object.updateAndPerformMeasurement(template_path, output_directory_root, output_file_name, tags, comment, globals_path,
#                                                                  globals_file_name, local_step_channels, looped_variables, local_parameters, local_instrument_values, log_channels=log_channels)
