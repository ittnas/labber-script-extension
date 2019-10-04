import script_object
import time
import numpy as np
import sys
sys.path.append(r'C:\Program Files (x86)\Labber\Script')
import Labber as Labber
import script_log_file
import os

class SimpleMeasurementObject(script_object.ScriptObject):
	def __init__(self,template_path,output_directory_root,globals_path,globals_file_name,comment = '',output_name=None):
		suffix = '_{:s}'.format(time.strftime('%H%M%S'))
		output_file_name = '_script scan resonator'
		if output_name is not None:
			output_file_name = output_name + output_file_name + suffix
		self.output_directory_root = output_directory_root
		self.output_file_name = output_file_name
		self.template_path = template_path
		self.globals_path = globals_path
		self.globals_file_name = globals_file_name
		self.comment = comment

	def updateAndPerformMeasurement(self,tags = [],local_step_channels={},looped_variables={},local_parameters={},local_instrument_values={},log_channels={}):
		(x,y) = script_object.updateAndPerformMeasurement(self.template_path,self.output_directory_root,self.output_file_name,tags,self.comment,self.globals_path,self.globals_file_name,local_step_channels,looped_variables,local_parameters,local_instrument_values,log_channels)
		return(x,y)

	def get_full_output_path(self):
		return script_object.get_full_output_path(self.output_directory_root,self.output_file_name)
		

	def scan_resonator(self,resonator_frequencies=None,bias_voltages=None,list_of_qubits=[1],digitizer_instrument_name='FPGA Digitizer',pulse_generator_instrument_name = 'Multi-Qubit Pulse Generator',qubit_control_rf_instrument_name_list = ['QB 1 RF'],trigger_period_channel_name='Trigger - Trig period',readout_rf_instrument_name='ReadOut RF', number_of_averages=100,integration_time=20e-6,readout_power = -18,extra_local_step_channels = {}):
		"""
		Scans the resonator as a function of probe tone and bias voltage.

		Arguments:
			resonator_frequencies (dict) -- A dictionary of channel name and sweep parameter pairs. Each pair corresponds to a single qubit. The first channel name acts as the master channel. Other channels are lookup tables to that channel.
		"""
		tags = ['Resonator Scan']
		# TODO allow equations to be set
		local_step_channels = {
			digitizer_instrument_name + ' - Number of averages':100,
			#'Multi-Qubit Pulse Generator - # of pulses':1,
			#'Multi-Qubit Pulse Generator - Alternate pulse direction':False,
			#'Multi-Qubit Pulse Generator - Plateau #1':50e-6,
			#'QB 1 RF - Power':-30
			pulse_generator_instrument_name + ' - # of pi pulses':-1,
			readout_rf_instrument_name + ' - Power':readout_power,
			#'Multi-Qubit Pulse Generator - Plateau #1':50E-6,
			#'FPGA Digitizer - Number of averages':1e5,
			#'QB 1 AWG - Ch4 - Offset':0.0,
			#'QB 1 RF - Power':-40,
			#'Discriminator - Training, input state':{'VARS':{'p1':'Multi-Qubit Pulse Generator - # of pulses'},'EQ':'p1','PARAM':{'use_relations':True}}
			#'Discriminator - Perform training':True
			#'Discriminator - Perform training':True,
			#'Discriminator - Use median value':True,
			digitizer_instrument_name +' - Integration time':integration_time
			   #'Trigger - Trig period':1.6e-6+500e-9
		}

		local_step_channels.update(extra_local_step_channels)
		buffer_time = 10e-6
		local_step_channels[trigger_period_channel_name] = integration_time + buffer_time

		for inst_name in qubit_control_rf_instrument_name_list:
			local_step_channels[inst_name + ' - Output'] = False

		looped_variables = {
#			'R1 freq offset':np.linspace(-10e6,10e6,21),
#			'QB1 bias offset':np.linspace(-100e-3,100e-3,21),
		}
		counter = 0
		master_channel_name = ''
		master_channel_values = None

		if resonator_frequencies is not None:
			for freq_channel_name,freq_channel_values in resonator_frequencies.items():
				if counter == 0:
					master_channel_name = freq_channel_name
					master_channel_values = freq_channel_values
					looped_variables[freq_channel_name] = freq_channel_values
				else:
					local_step_channels[freq_channel_name] = {'VARS':{'p1':{'channel_name':master_channel_name,'lookup_x':master_channel_values,'lookup_y':freq_channel_values}},'EQ':'p1'}

				counter = counter+1

		master_channel_name = ''
		master_channel_values = None

		if bias_voltages is not None:
			counter = 0
			for bias_channel_name,bias_channel_values in bias_voltages.items():
				if counter == 0:
					master_channel_name = bias_channel_name
					master_channel_values = bias_channel_values
					looped_variables[bias_channel_name] = bias_channel_values
				else:
					local_step_channels[bias_channel_name] = {'VARS':{'p1':{'channel_name':master_channel_name,'lookup_x':master_channel_values,'lookup_y':bias_channel_values}},'EQ':'p1'}

				counter = counter+1

		local_parameters = {'hardware_loop':False}

		local_instrument_values = {
			pulse_generator_instrument_name:{'Sequence':1},
			#'Discriminator':{'Training type':0},
			#'Discriminator':{'Perform training':False},
			#'Discriminator':{'Pointer, QB1-S0':29.1E-6-9.19E-3j},
			#'Discriminator':{'Pointer, QB1-S1':3.84E-3-4.08E-3j},
			}

		log_channels = [
			]

		for QB in list_of_qubits:
			log_channels.append(digitizer_instrument_name + ' - FPGA Voltage, QB'+str(QB))

		# print(log_channels)
		# print(local_instrument_values)
		# print(local_step_channels)
		
		self.updateAndPerformMeasurement(tags = tags,local_step_channels=local_step_channels,looped_variables=looped_variables,local_parameters=local_parameters,local_instrument_values=local_instrument_values,log_channels=log_channels)
		#if len(resonator_frequencies) > 1:
		return SimpleMeasurementObject.get_resonator_frequencies(self.get_full_output_path(),list_of_qubits=list_of_qubits,log_channel_name = 'FPGA Digitizer - FPGA Voltage, QB',analysis_method='min')
		#else:
	#	return None

	def get_resonator_frequencies(log_file_path,list_of_qubits=[1],log_channel_name = 'FPGA Digitizer - FPGA Voltage, QB',analysis_method='min'):
		output_file = script_log_file.ScriptLogFile(log_file_path)
		resonator_frequencies_all_qb = []
		channel_values_all_qb = []
		for qb in list_of_qubits:
			(data,channel_values,channel_names) = output_file.getDataMatrix(log_channel_name + str(qb))
			
	
			resonator_frequencies = []
			if analysis_method is 'min':
				if len(channel_values)>1:
					for ii in range(len(channel_values[1])):
						resonator_frequencies.append(channel_values[0][np.argmin(np.abs(data[:,ii]))])
					resonator_frequencies_all_qb.append(resonator_frequencies)
					channel_values_all_qb.append(channel_values[1])
						
				else:
					resonator_frequencies = [channel_values[0][np.argmin(np.abs(data))]]
					resonator_frequencies_all_qb.append(resonator_frequencies)
					channel_values_all_qb.append(None)
					#return  (resonator_frequencies,None)
		return resonator_frequencies_all_qb,channel_values_all_qb



