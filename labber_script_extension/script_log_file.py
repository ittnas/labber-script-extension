import numpy as np
import sys
import os
sys.path.append(os.environ.get(
    'LABBERPATH', r'C:\Program Files (x86)\Labber\Script'))  # NOQA: E402
import Labber as Labber
from Labber import ScriptTools


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

        channel_names = []
        channel_values = []
        dimensions = []
        for channel in step_channels:
            if isinstance(channel['values'], np.ndarray) and len(channel['values']) > 1:
                # print(len(channel['values']))
                channel_names.append(channel['name'])
                channel_values.append(channel['values'])
                dimensions.append(len(channel['values']))

        data = np.transpose(self.getData(log_channel))
        # print(data.shape)
        # print(dimensions)
        data = np.reshape(data, dimensions)

        return (data, channel_values, channel_names)
        # print(data.shape)
