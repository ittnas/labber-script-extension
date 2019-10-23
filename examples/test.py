import numpy as np
#import labber_script_extension.script_object
import labber_script_extension.script_object
# export LABBERPATH = /usr/share/Labber/Script/

globals_path = r'./'
globals_file_name = 'globals_test'
template_path = './test_master3.labber'
output_directory_root = os.getcwd()
output_file_name = 'test5'
tags = []
comment = 'no comments'
local_step_channels = {}
looped_variables = {}
local_parameters = {}
local_instrument_values = {}
log_channels = ['Manual - Value 3']

so = labber_script_extension.script_object.ScriptObject(
    template_path, os.path.join(output_directory_root, '2019/10/Data_1018/', output_file_name))

# so.printStepChannels()
# so.getStepChannels()
# print(so.getChannelNames())
# print(so.getLogChannels())
so.setLogChannels('Manual - Value 3')
# so.add
# print(so.getLogChannels())
#so.setSignalConnection('FE - State vector, QB1', 'multi - Trace - I1')

so.setSignalConnectionsByDict({'FE - State vector, QB1': 'multi - Trace - I1'})
so.getStepChannelNames()
so.addStepChannel('Manual - Value 1')
so.getStepChannelNames()
so.getChannelNames()
so.getInstrumentValueNames('tau', None, is_case_sensitive=False)

# so.performMeasurement()
# labber_script_extension.script_object.updateAndPerformMeasurement(template_path, output_directory_root, output_file_name, tags, comment, globals_path,
#                                                                  globals_file_name, local_step_channels, looped_variables, local_parameters, local_instrument_values, log_channels=log_channels)
