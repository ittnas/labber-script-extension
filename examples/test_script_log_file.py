import numpy as np
# import labber_script_extension.script_object
from labber_script_extension import script_object
from labber_script_extension import script_log_file
# export LABBERPATH = /usr/share/Labber/Script/

template_path = './test_master3.labber'
output_directory_root = os.getcwd()
output_file_name = 'test_log_file_6'
tags = ['LogFile test']
comment = 'no comments'
local_step_channels = {}
looped_variables = {}
local_parameters = {}
local_instrument_values = {}
log_channels = ['Manual - Value 3']

so = script_object.ScriptObject(
    template_path, os.path.join(output_directory_root, '2019/10/Data_1018/', output_file_name))

so.setLogChannels(['FE - Frequency estimate, QB1', 'multi - Trace - I1'])
so.setSignalConnection('FE - State vector, QB1', 'multi - Trace - I1')
so.updateStepChannelsByDict(
    {'multi - Plateau': np.linspace(0, 10e-9, 3), 'multi - Width': np.linspace(0, 10e-9, 4), 'Manual - Value 1': np.linspace(0, 10e-9, 5)},)

so.performMeasurement()

lf = script_log_file.ScriptLogFile(so.getOutputPathOfPreviousMeasurement())
data = lf.getDataMatrix('FE - Frequency estimate, QB1')
print(data[0])
