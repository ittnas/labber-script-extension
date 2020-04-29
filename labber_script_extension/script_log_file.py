import numpy as np
import sys
import os
sys.path.append(os.environ.get(
    'LABBERPATH', r'C:\Program Files (x86)\Labber\Script'))  # NOQA: E402
import Labber as Labber


class ScriptLogFile(Labber.LogFile):
    """
        Extension of Labber LogFile that provides some additional features for reading the data.
    """

    def getDataMatrix(self, log_channel=None):
        """
        Returns the data as a numpy array.

        Parameters:
            log_channel (str) -- name of the log channel which data is returned.
        """
        step_channels = self.getStepChannels()
        log_channels = self.getLogChannels()
        if log_channel is None:
            log_channel = log_channels[0]['name']

        log_channel_dict = None
        for log_channel_dict_iter in log_channels:
            if log_channel_dict_iter['name'] == log_channel:
                log_channel_dict = log_channel_dict_iter
                break
        channel_names = []
        channel_values = []
        dimensions = []
        for channel in step_channels:
            if isinstance(channel['values'], np.ndarray) and len(channel['values']) > 1:
                channel_names.append(channel['name'])
                channel_values.append(channel['values'])
                dimensions.append(len(channel['values']))

        data = self.getData(log_channel)
        if log_channel_dict['vector']:
            dimensions.insert(0, data.shape[-1])
            channel_names.insert(0, 'vector_data')
            channel_values.insert(0, np.linspace(0, 1, data.shape[-1]))

        first_element = dimensions.pop(0)
        dimensions = dimensions[::-1]
        dimensions.append(first_element)

        data = np.reshape(data, dimensions)
        data = np.transpose(data)

        return (data, channel_values, channel_names)
