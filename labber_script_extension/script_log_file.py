import sys
sys.path.append(r'C:\Program Files (x86)\Labber\Script')
import Labber as Labber
from Labber import ScriptTools
import numpy as np

class ScriptLogFile(Labber.LogFile):
	def getDataMatrix(self,log_channel=None):
		step_channels = self.getStepChannels()

		channel_names = []
		channel_values = []
		dimensions = []
		for channel in step_channels:
			if isinstance(channel['values'],np.ndarray) and len(channel['values']) > 1:
				#print(len(channel['values']))
				channel_names.append(channel['name'])
				channel_values.append(channel['values'])
				dimensions.append(len(channel['values']))

		data = np.transpose(self.getData(log_channel))
		#print(data.shape)
		#print(dimensions)
		data = np.reshape(data,dimensions)

		return (data,channel_values,channel_names)
		#print(data.shape)
